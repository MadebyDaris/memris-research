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
