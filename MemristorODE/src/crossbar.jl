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