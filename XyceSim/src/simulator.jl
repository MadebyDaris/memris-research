"""
    XyceSim.Simulator

Functions for creating and configuring Xyce simulator instances via Jyce.
"""
module Simulator

using Jyce
using Printf

"""
    create_simulator(; sub_path::String = "./memristor_model.sub", plugin_path::String = "")

Create and configure a simulator instance with an optional subsckt library.
"""
function create_simulator(; sub_path::String = "./memristor_model.sub", plugin_path::String = "")
    # Ensure native available before pure-Jyce run
    if !Jyce.native_available()
        error("Jyce native backend not available. Run Cell 1 diagnostics and set JYCE_XYCESOLVER_JULIA_LIB or JYCE_XYCESOLVER_ROOT if needed.")
    end

    sub_candidates = [joinpath(@__DIR__, "biolek_mmrstor_model.sub"), "biolek_mmrstor_model.sub", sub_path]
    sub_idx = findfirst(isfile, sub_candidates)
    if sub_idx === nothing
        error("biolek_mmrstor_model.sub not found. Place it next to this script or pass sub_path=...")
    end
    sub = abspath(sub_candidates[sub_idx])
    println("Using subcircuit: ", sub)

    sim = try Jyce.XyceSimulator(false) catch e
        @error "Could not create XyceSimulator: " e
        nothing
    end
    @assert sim !== nothing "Failed to construct Jyce.XyceSimulator; see diagnostics"

    # Optionally load memristor plugin if available
    plugin_candidates = filter(!isempty, [
        plugin_path,
        joinpath(@__DIR__, "..", "Jyce", "plugins", "memristor_plugin.so"),
        joinpath(@__DIR__, "..", "Jyce", "plugins", "build", "libmemristor_plugin.so"),
    ])
    plugin_idx = findfirst(isfile, plugin_candidates)
    if plugin_idx !== nothing
        plugin = plugin_candidates[plugin_idx]
        @info "Loading memristor plugin" plugin
        Jyce.add_plugin_library(sim, plugin)
    else
        @warn "Memristor plugin not found; crossbar simulation may fail" plugin_candidates
    end

    return sim
end

end # module