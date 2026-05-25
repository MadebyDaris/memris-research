using Pkg; Pkg.activate("memres/memres/julia_poc")
include("memres/memres/julia_poc/src/MemristorModels.jl")
using .MemristorModels

vp = VTEAMParams(
    v_on = 0.5, v_off = -0.5,
    k_on = -2e-7, k_off = 2e-7,
    alpha_on = 3.0, alpha_off = 3.0,
    w_on = 0.0, w_off = 10e-9,
    R_on = 1e3, R_off = 100e3, w_init = 8e-9
)
V_amplitude = 1.2
T_period = 1e-2
voltage_func(t) = sinusoidal_wave(t; amplitude=V_amplitude, frequency=1 / T_period)
tspan = (0.0, 6 * T_period)
res = simulate_vteam(vp, tspan, voltage_func)
println("Min w: ", minimum(res.w), " Max w: ", maximum(res.w))
