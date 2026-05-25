using Pkg; Pkg.activate("memres/memres/julia_poc")
include("memres/memres/julia_poc/src/MemristorModels.jl")
using .MemristorModels

tp = ThresholdParams(
    R_on = 1e3, R_off = 100e3,
    Vth_p = 0.5, Vth_n = -0.5,
    k_on = 25.0, k_off = 25.0,
    w_init = 0.5, p_window = 1
)

voltage_func(t) = sinusoidal_wave(t; amplitude=1.5, frequency=1.0)
res = simulate_memristor(tp, (0.0, 2.0), voltage_func; saveat=range(0, 2, length=1000))

# Extract last cycle
last_cycle_start = 1.0
idx = findall(t -> t >= last_cycle_start, res.t)

println("Last cycle w range: ", minimum(res.w[idx]), " to ", maximum(res.w[idx]))
println("Last cycle I range: ", minimum(res.I[idx]) * 1e3, " to ", maximum(res.I[idx]) * 1e3, " mA")

# Check negative half
idx_neg = findall(t -> (t >= last_cycle_start) && (res.V[findfirst(x->x==t, res.t)] < 0), res.t)
if length(idx_neg) > 0
    println("Negative half w range: ", minimum(res.w[idx_neg]), " to ", maximum(res.w[idx_neg]))
    println("Negative half I range: ", minimum(res.I[idx_neg]) * 1e3, " to ", maximum(res.I[idx_neg]) * 1e3, " mA")
end
