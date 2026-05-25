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
