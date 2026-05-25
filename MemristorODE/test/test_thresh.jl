using Pkg; Pkg.activate("memres/memres/julia_poc")
include("memres/memres/julia_poc/src/MemristorModels.jl")
using .MemristorModels
using Printf

for exp_val in [0, 1, 2, 3, 4]
    tp = ThresholdParams(
        R_on = 1e3, R_off = 100e3,
        Vth_p = 0.5, Vth_n = -0.5,
        k_on = 10.0^exp_val, k_off = 10.0^exp_val,
        w_init = 0.5, p_window = 1
    )
    res = simulate_memristor(tp, (0.0, 6e-2), t -> sinusoidal_wave(t; amplitude=1.2, frequency=100.0))
    @printf "Threshold 10^%d -> w range: %g to %g\n" exp_val minimum(res.w) maximum(res.w)
end
