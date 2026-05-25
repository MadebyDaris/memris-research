"""
    XyceSim

Julia interface for Xyce simulator via Jyce wrapper for memristor crossbar simulation.
"""
module XyceSim

using Jyce
using Printf

# Export types and functions from submodules
export MemristorParams, CrossbarArray
export create_simulator
export create_crossbar_array, create_crossbar_array_subckt
export write_netlist
export run_crossbar, run_and_plot_crossbar

# Include sub-modules
include("simulator.jl")
include("netlist.jl")
include("postprocess.jl")

end # module