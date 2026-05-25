"""
    MemristorODE

Comprehensive memristor models for pure Julia ODE simulation.
Includes threshold-based and VTEAM models, crossbar array simulation with
non-idealities (IR drop, variability, noise), and sparse MNA solver.
"""
module MemristorODE

using DifferentialEquations
using LinearAlgebra
using SparseArrays
using Random
using Statistics
using Printf

# Export types and parameters
export ThresholdParams, VTEAMParams, CrossbarArray

# Export ODE and dynamics functions
export memristor_ode!, vteam_dynamics!

# Export simulation functions
export simulate_memristor, simulate_vteam

# Export model helpers
export resistance, conductance, current

# Export window functions
export window_joglekar, window_biolek

# Export signals
export triangular_wave, sinusoidal_wave, pulse_train

# Export crossbar functions
export simulate_crossbar_mvm, simulate_crossbar_mvm_with_ir
export solve_crossbar_mna, build_crossbar_mna_matrix
export ideal_mvm, compute_mvm_error
export get_ir_drop_map

# Include sub-modules
include("models/windows.jl")
include("models/threshold.jl")
include("models/vteam.jl")
include("signals.jl")
include("crossbar.jl")

end # module