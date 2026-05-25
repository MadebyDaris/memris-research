"""
    MemristorModels

Comprehensive memristor models for hybrid Julia-Xyce simulation framework.
Includes threshold-based and VTEAM models, crossbar array simulation with
non-idealities (IR drop, variability, noise), and sparse MNA solver.
"""
module MemristorModels

using DifferentialEquations
using LinearAlgebra
using SparseArrays
using Random

export ThresholdParams, VTEAMParams
export memristor_ode!, vteam_dynamics!
export simulate_memristor, simulate_vteam
export resistance, current, triangular_wave, sinusoidal_wave
export window_joglekar, window_biolek
export CrossbarArray, simulate_crossbar_mvm, simulate_crossbar_mvm_with_ir
export solve_crossbar_mna, build_crossbar_mna_matrix
export ideal_mvm, compute_mvm_error
export get_ir_drop_map

# ============================================================================
# Window Functions
# ============================================================================

"""
    window_joglekar(w, w_on, w_off, p=1)

Joglekar window function to bound state variable.
f(w) = 1 - ((2w - w_on - w_off)/(w_off - w_on))^(2p)
"""
function window_joglekar(w::Real, w_on::Real, w_off::Real, p::Int=1)
    w_norm = (2w - w_on - w_off) / (w_off - w_on)
    return 1 - w_norm^(2p)
end

"""
    window_biolek(w, w_on, w_off, dw, p=1)

Biolek window function (direction-dependent).
"""
function window_biolek(w::Real, w_on::Real, w_off::Real, dw::Real, p::Int=1)
    if dw >= 0
        return 1 - ((w - w_on) / (w_off - w_on))^(2p)
    else
        return 1 - ((w_off - w) / (w_off - w_on))^(2p)
    end
end

# ============================================================================
# Threshold-Based Model
# ============================================================================

"""
    ThresholdParams

Parameters for threshold-based memristor model.
"""
Base.@kwdef struct ThresholdParams
    R_on::Float64 = 1e3       # ON resistance (Ω)
    R_off::Float64 = 1e5      # OFF resistance (Ω)
    Vth_p::Float64 = 0.5      # Positive threshold voltage (V)
    Vth_n::Float64 = -0.5     # Negative threshold voltage (V)
    k_on::Float64 = 1e4       # ON switching rate (1/s)
    k_off::Float64 = 1e4      # OFF switching rate (1/s)
    w_init::Float64 = 0.5     # Initial state [0, 1]
    p_window::Int = 1         # Window function parameter
end

"""
    resistance(w, p::ThresholdParams)
Calculate memristor resistance for normalized state w ∈ [0, 1].
"""
function resistance(w::Real, p::ThresholdParams)
    return p.R_on * w + p.R_off * (1 - w)
end

"""
    current(V, w, p::ThresholdParams)
Calculate memristor current for voltage V and state w.
"""
function current(V::Real, w::Real, p::ThresholdParams)
    R = resistance(w, p)
    return V / R
end
"""
    memristor_ode!(dw, w, params, t; voltage_func)
ODE for threshold-based memristor state evolution using Biolek window.
"""
function memristor_ode!(dw, w, params::ThresholdParams, t; voltage_func)
    p = params
    V = voltage_func(t)

    # Calculate raw derivative (using linear overdrive for continuity at threshold)
    dw_raw = 0.0
    
    if V > p.Vth_p
        # SET operation (V - Vth_p is positive -> dw_raw is positive)
        dw_raw = p.k_on * (V - p.Vth_p)
        
    elseif V < p.Vth_n
        # RESET operation (V - Vth_n is negative -> dw_raw is negative)
        dw_raw = p.k_off * (V - p.Vth_n)
        
    else
        # Sub-threshold region
        dw[1] = 0.0
        return
    end

    # Apply direction-dependent Biolek window function
    fw = window_biolek(w[1], 0.0, 1.0, dw_raw, p.p_window)
    dw[1] = dw_raw * fw

    # Hard clamp derivative at boundaries to assist the solver
    if w[1] <= 0.0 && dw[1] < 0.0
        dw[1] = 0.0
    elseif w[1] >= 1.0 && dw[1] > 0.0
        dw[1] = 0.0
    end
end

"""
    simulate_memristor(params::ThresholdParams, tspan, voltage_func; kwargs...)
Simulate single memristor device using threshold model.
"""
function simulate_memristor(params::ThresholdParams, tspan::Tuple, voltage_func;
                           solver=Rodas5(), abstol=1e-8, reltol=1e-6,
                           saveat=nothing)
    w0 = [params.w_init]

    ode_func!(dw, w, p, t) = memristor_ode!(dw, w, p, t; voltage_func=voltage_func)

    prob = ODEProblem(ode_func!, w0, tspan, params)

    if saveat === nothing
        saveat = range(tspan[1], tspan[2], length=1000)
    end

    sol = solve(prob, solver; abstol=abstol, reltol=reltol, saveat=saveat)

    t = sol.t
    w = [s[1] for s in sol.u]
    V = voltage_func.(t)
    
    # CRITICAL FIX: Clamp the state variable array before resistance calculation
    # This guarantees no division-by-zero or negative resistance spikes
    w_safe = clamp.(w, 0.0, 1.0)
    
    R = resistance.(w_safe, Ref(params))
    I = V ./ R

    return (sol=sol, t=t, w=w, V=V, I=I, R=R)
end

# ============================================================================
# VTEAM Model (Industry Standard)
# ============================================================================

"""
    VTEAMParams

Parameters for VTEAM (Voltage ThrEshold Adaptive Memristor) model.
This is the industry-standard model for system-level simulation.
"""
Base.@kwdef struct VTEAMParams
    v_on::Float64 = 0.5       # Threshold voltage ON (V)
    v_off::Float64 = -0.5     # Threshold voltage OFF (V)
    k_on::Float64 = -200.0    # Settling rate ON (m/s) - negative for decreasing w
    k_off::Float64 = 200.0    # Settling rate OFF (m/s) - positive for increasing w
    alpha_on::Float64 = 3.0   # Nonlinearity exponent ON
    alpha_off::Float64 = 3.0  # Nonlinearity exponent OFF
    w_on::Float64 = 0.0       # Min state (m) - conducting
    w_off::Float64 = 10e-9    # Max state (m) - insulating
    R_on::Float64 = 1e3       # LRS resistance (Ω)
    R_off::Float64 = 1e6      # HRS resistance (Ω)
    w_init::Float64 = 5e-9    # Initial state (m)
end

"""
    resistance(w, p::VTEAMParams)

Calculate VTEAM memristor resistance from physical state w.
R(w) = R_on + (w/w_off) * (R_off - R_on)
"""
function resistance(w::Real, p::VTEAMParams)
    w_clamped = clamp(w, p.w_on, p.w_off)
    w_norm = (w_clamped - p.w_on) / (p.w_off - p.w_on)
    return p.R_on + w_norm * (p.R_off - p.R_on)
end

"""
    conductance(w, p::VTEAMParams)

Calculate VTEAM memristor conductance from physical state w.
"""
function conductance(w::Real, p::VTEAMParams)
    return 1.0 / resistance(w, p)
end

"""
    vteam_dynamics!(dw, w, params, t; voltage_func)

VTEAM state derivative function dw/dt.
"""
function vteam_dynamics!(dw, w, params::VTEAMParams, t; voltage_func)
    p = params
    V = voltage_func(t)

    # Clamp state to physical boundaries
    w_current = clamp(w[1], p.w_on, p.w_off)

    # Window function (Joglekar)
    fw = window_joglekar(w_current, p.w_on, p.w_off, 1)

    if V > p.v_on && w_current > p.w_on
        # SET: decreasing w (toward conducting state)
        dw[1] = p.k_on * ((V / p.v_on) - 1)^p.alpha_on * fw
    elseif V < p.v_off && w_current < p.w_off
        # RESET: increasing w (toward insulating state)
        dw[1] = p.k_off * ((V / p.v_off) - 1)^p.alpha_off * fw
    else
        # Sub-threshold
        dw[1] = 0.0
    end

    # Enforce boundaries
    if w[1] <= p.w_on && dw[1] < 0.0
        dw[1] = 0.0
    elseif w[1] >= p.w_off && dw[1] > 0.0
        dw[1] = 0.0
    end
end

"""
    simulate_vteam(params::VTEAMParams, tspan, voltage_func; kwargs...)

Simulate single memristor device using VTEAM model.
"""
function simulate_vteam(params::VTEAMParams, tspan::Tuple, voltage_func;
                        solver=Rodas5P(), abstol=1e-10, reltol=1e-8,
                        saveat=nothing)
    w0 = [params.w_init]

    ode_func!(dw, w, p, t) = vteam_dynamics!(dw, w, p, t; voltage_func=voltage_func)

    prob = ODEProblem(ode_func!, w0, tspan, params)

    if saveat === nothing
        saveat = range(tspan[1], tspan[2], length=1000)
    end

    sol = solve(prob, solver; abstol=abstol, reltol=reltol, saveat=saveat)

    t = sol.t
    w = [s[1] for s in sol.u]
    V = voltage_func.(t)
    R = resistance.(w, Ref(params))
    I = V ./ R

    return (sol=sol, t=t, w=w, V=V, I=I, R=R)
end

# ============================================================================
# Voltage Waveforms
# ============================================================================

"""
    triangular_wave(t; amplitude, period)

Generate triangular wave voltage signal.
"""
function triangular_wave(t; amplitude=1.0, period=1e-3)
    phase = mod(t, period) / period
    if phase < 0.25
        return 4 * amplitude * phase
    elseif phase < 0.75
        return amplitude * (2 - 4 * phase)
    else
        return 4 * amplitude * (phase - 1)
    end
end

"""
    sinusoidal_wave(t; amplitude, frequency)

Generate sinusoidal voltage signal.
"""
function sinusoidal_wave(t; amplitude=1.0, frequency=1e3)
    return amplitude * sin(2π * frequency * t)
end

"""
    pulse_train(t; V_high, V_low, t_high, t_low, t_rise=0.0)

Generate pulse train for SET/RESET programming.
"""
function pulse_train(t; V_high=1.0, V_low=-1.0, t_high=1e-6, t_low=1e-6, t_rise=0.0)
    T = t_high + t_low + 2*t_rise
    phase = mod(t, T)

    if phase < t_rise
        return V_low + (V_high - V_low) * (phase / t_rise)
    elseif phase < t_rise + t_high
        return V_high
    elseif phase < 2*t_rise + t_high
        return V_high + (V_low - V_high) * ((phase - t_rise - t_high) / t_rise)
    else
        return V_low
    end
end

# ============================================================================
# Crossbar Array Modeling
# ============================================================================

"""
    CrossbarArray

Memristive crossbar array structure with non-ideality parameters.
"""
struct CrossbarArray
    rows::Int
    cols::Int
    conductances::Matrix{Float64}    # G = 1/R matrix (S)
    R_wire::Float64                  # Wire resistance per segment (Ω)
    variability::Float64             # Device-to-device variation (σ/μ)
    noise_level::Float64             # Read noise level (σ/μ)
end

"""
    CrossbarArray(rows, cols; R_on, R_off, R_wire, variability, noise_level)

Create crossbar array with random conductances and specified non-idealities.
"""
function CrossbarArray(rows::Int, cols::Int;
                      R_on=1e3, R_off=1e6,
                      R_wire=0.0,
                      variability=0.0,
                      noise_level=0.0)
    G_on = 1/R_on
    G_off = 1/R_off

    # Random base conductances (uniform distribution)
    G_base = G_off .+ (G_on - G_off) .* rand(rows, cols)

    # Apply device-to-device variability
    if variability > 0
        G_var = G_base .* (1 .+ variability .* randn(rows, cols))
        G_array = clamp.(G_var, G_off, G_on)
    else
        G_array = G_base
    end

    return CrossbarArray(rows, cols, G_array, R_wire, variability, noise_level)
end

"""
    CrossbarArray(rows, cols, G::Matrix; R_wire, variability, noise_level)

Create crossbar array with specified conductance matrix.
"""
function CrossbarArray(rows::Int, cols::Int, G::Matrix{Float64};
                      R_wire=0.0, variability=0.0, noise_level=0.0)
    @assert size(G) == (rows, cols) "Conductance matrix size mismatch"
    return CrossbarArray(rows, cols, copy(G), R_wire, variability, noise_level)
end

"""
    ideal_mvm(G::Matrix, V_in::Vector)

Compute ideal matrix-vector multiplication: I = G' * V
"""
function ideal_mvm(G::Matrix{Float64}, V_in::Vector{Float64})
    return G' * V_in
end

"""
    simulate_crossbar_mvm(xbar::CrossbarArray, V_in::Vector; apply_noise=true)

Simulate MVM with device variability and read noise (no IR drop).
"""
function simulate_crossbar_mvm(xbar::CrossbarArray, V_in::Vector{Float64};
                               apply_noise=true)
    @assert length(V_in) == xbar.rows "Input vector length must match crossbar rows"

    # Ideal MVM
    I_out = xbar.conductances' * V_in

    # Add read noise
    if apply_noise && xbar.noise_level > 0
        noise = randn(length(I_out)) .* (abs.(I_out) .* xbar.noise_level)
        I_out = I_out .+ noise
    end

    return I_out
end

# ============================================================================
# IR Drop Analysis using Sparse MNA
# ============================================================================

"""
    build_crossbar_mna_matrix(rows, cols, G, R_wire)

Build the Modified Nodal Analysis (MNA) matrix for crossbar with wire resistance.

The crossbar is modeled as a resistor network where:
- Each crosspoint has a memristor connecting wordline node to bitline node
- Wordlines and bitlines have series resistance R_wire between adjacent nodes

Returns: (A, b_template) where A*x = b gives node voltages
"""
function build_crossbar_mna_matrix(rows::Int, cols::Int,
                                   G::Matrix{Float64}, R_wire::Float64)
    # Node indexing:
    # - Wordline nodes: 1 to rows*cols (node (i,j) = (i-1)*cols + j)
    # - Bitline nodes: rows*cols+1 to 2*rows*cols

    n_nodes = 2 * rows * cols
    G_wire = 1.0 / R_wire

    # Build conductance matrix using sparse arrays
    I_idx = Int[]
    J_idx = Int[]
    V_vals = Float64[]

    function wl_node(i, j)
        return (i-1)*cols + j
    end

    function bl_node(i, j)
        return rows*cols + (j-1)*rows + i
    end

    # Add memristor conductances
    for i in 1:rows
        for j in 1:cols
            wl = wl_node(i, j)
            bl = bl_node(i, j)
            g = G[i, j]

            # Self conductance
            push!(I_idx, wl); push!(J_idx, wl); push!(V_vals, g)
            push!(I_idx, bl); push!(J_idx, bl); push!(V_vals, g)

            # Cross conductance
            push!(I_idx, wl); push!(J_idx, bl); push!(V_vals, -g)
            push!(I_idx, bl); push!(J_idx, wl); push!(V_vals, -g)
        end
    end

    # Add wordline wire resistances (horizontal)
    for i in 1:rows
        for j in 1:(cols-1)
            n1 = wl_node(i, j)
            n2 = wl_node(i, j+1)

            push!(I_idx, n1); push!(J_idx, n1); push!(V_vals, G_wire)
            push!(I_idx, n2); push!(J_idx, n2); push!(V_vals, G_wire)
            push!(I_idx, n1); push!(J_idx, n2); push!(V_vals, -G_wire)
            push!(I_idx, n2); push!(J_idx, n1); push!(V_vals, -G_wire)
        end
    end

    # Add bitline wire resistances (vertical)
    for j in 1:cols
        for i in 1:(rows-1)
            n1 = bl_node(i, j)
            n2 = bl_node(i+1, j)

            push!(I_idx, n1); push!(J_idx, n1); push!(V_vals, G_wire)
            push!(I_idx, n2); push!(J_idx, n2); push!(V_vals, G_wire)
            push!(I_idx, n1); push!(J_idx, n2); push!(V_vals, -G_wire)
            push!(I_idx, n2); push!(J_idx, n1); push!(V_vals, -G_wire)
        end
    end

    # Add ground connection for bitline outputs (large conductance)
    G_ground = 1e6  # Virtual ground conductance
    for j in 1:cols
        bl = bl_node(rows, j)  # Bottom of bitline
        push!(I_idx, bl); push!(J_idx, bl); push!(V_vals, G_ground)
    end

    A = sparse(I_idx, J_idx, V_vals, n_nodes, n_nodes)

    return A
end

"""
    solve_crossbar_mna(xbar::CrossbarArray, V_in::Vector)

Solve crossbar MNA equations to get output currents with IR drop.
"""
function solve_crossbar_mna(xbar::CrossbarArray, V_in::Vector{Float64})
    @assert length(V_in) == xbar.rows "Input vector length must match rows"

    rows, cols = xbar.rows, xbar.cols

    if xbar.R_wire ≈ 0.0
        # No wire resistance - use ideal MVM
        return simulate_crossbar_mvm(xbar, V_in; apply_noise=true)
    end

    # Build MNA matrix
    A = build_crossbar_mna_matrix(rows, cols, xbar.conductances, xbar.R_wire)

    # Build RHS vector (current injection at wordline inputs)
    n_nodes = 2 * rows * cols
    b = zeros(n_nodes)

    G_source = 1e6  # Source conductance (voltage source approximation)

    # Inject current at first column of each wordline (V_in * G_source)
    for i in 1:rows
        wl_node = (i-1)*cols + 1
        b[wl_node] = V_in[i] * G_source
        # Also add source conductance to diagonal
        A[wl_node, wl_node] += G_source
    end

    # Solve sparse system
    x = A \ b

    # Extract output currents from bitline bottom nodes
    I_out = zeros(cols)
    G_ground = 1e6
    for j in 1:cols
        bl_node = rows*cols + (j-1)*rows + rows
        I_out[j] = x[bl_node] * G_ground  # Current to ground
    end

    # Add read noise
    if xbar.noise_level > 0
        noise = randn(cols) .* (abs.(I_out) .* xbar.noise_level)
        I_out = I_out .+ noise
    end

    return I_out
end

"""
    simulate_crossbar_mvm_with_ir(xbar::CrossbarArray, V_in::Vector)

Simulate MVM with full IR drop analysis using sparse MNA solver.
"""
function simulate_crossbar_mvm_with_ir(xbar::CrossbarArray, V_in::Vector{Float64})
    return solve_crossbar_mna(xbar, V_in)
end

"""
    compute_mvm_error(I_actual::Vector, I_ideal::Vector)

Compute relative MVM error metrics.
"""
function compute_mvm_error(I_actual::Vector{Float64}, I_ideal::Vector{Float64})
    abs_error = norm(I_actual - I_ideal)
    rel_error = abs_error / norm(I_ideal)

    max_abs_error = maximum(abs.(I_actual - I_ideal))
    max_rel_error = maximum(abs.(I_actual - I_ideal) ./ (abs.(I_ideal) .+ 1e-12))

    return (
        absolute = abs_error,
        relative = rel_error,
        max_absolute = max_abs_error,
        max_relative = max_rel_error,
        rmse = sqrt(mean((I_actual - I_ideal).^2))
    )
end

"""
    get_ir_drop_map(xbar::CrossbarArray, V_in::Vector)

Compute IR drop across the crossbar array.
Returns voltage drop at each crosspoint relative to ideal.
"""
function get_ir_drop_map(xbar::CrossbarArray, V_in::Vector{Float64})
    @assert length(V_in) == xbar.rows

    if xbar.R_wire ≈ 0.0
        return zeros(xbar.rows, xbar.cols)
    end

    rows, cols = xbar.rows, xbar.cols

    # Build and solve MNA
    A = build_crossbar_mna_matrix(rows, cols, xbar.conductances, xbar.R_wire)
    n_nodes = 2 * rows * cols
    b = zeros(n_nodes)
    G_source = 1e6

    for i in 1:rows
        wl_node = (i-1)*cols + 1
        b[wl_node] = V_in[i] * G_source
        A[wl_node, wl_node] += G_source
    end

    x = A \ b

    # Extract wordline voltages and compute drops
    V_drop = zeros(rows, cols)
    for i in 1:rows
        V_ideal = V_in[i]
        for j in 1:cols
            wl_node = (i-1)*cols + j
            V_drop[i, j] = V_ideal - x[wl_node]
        end
    end

    return V_drop
end

end # module
