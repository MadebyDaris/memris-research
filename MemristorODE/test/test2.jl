using Pkg; Pkg.activate("memres/memres/julia_poc")
include("memres/memres/julia_poc/src/MemristorModels.jl")
using .MemristorModels
using Printf

# Try scanning k_on
for exp_val in [-9, -8, -7, -6]
    vp = VTEAMParams(v_on=0.5, v_off=-0.5, k_on=-10.0^exp_val, k_off=10.0^exp_val, 
                     w_off=10e-9, R_on=1e3, R_off=100e3, w_init=5e-9)
    res = simulate_vteam(vp, (0.0, 1e-2), t -> sinusoidal_wave(t; amplitude=1.0, frequency=100.0))
    @printf "VTEAM 10^%d -> w range: %g to %g\n" exp_val minimum(res.w) maximum(res.w)
end
