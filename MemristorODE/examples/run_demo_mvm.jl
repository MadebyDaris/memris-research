# Demo script: program a known weight pattern and validate MVM using Xyce (via Jyce)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Jyce
using Random
using Printf

include("crossbar.jl")

function build_resistor_crossbar_netlist(G::Matrix{Float64}; Vdrive=0.1, rload=10e3, tstep=1e-6, tstop=1e-4)
    rows, cols = size(G)
    io = IOBuffer()
    println(io, "* Resistor crossbar netlist for MVM validation")
    println(io)
    # Row sources
    for i in 1:rows
        println(io, @sprintf("VROW_%d in_%d 0 PULSE(0 %0.6f 0 1n 1n %0.8f %0.8f)", i, i, Vdrive, tstep, tstop))
    end
    # Column loads
    for j in 1:cols
        println(io, @sprintf("RLOAD_%d out_%d 0 %0.8f", j, j, rload))
    end
    println(io)
    # Crosspoint resistors: R = 1/G (handle zero conductance as large resistor)
    for i in 1:rows
        for j in 1:cols
            Rval = G[i,j] > 0 ? 1.0 / G[i,j] : 1e12
            println(io, @sprintf("R_%d_%d in_%d out_%d %0.8f", i, j, i, j, Rval))
        end
    end
    println(io)
    println(io, ".TRAN $(tstep) $(tstop)")
    # Print all voltages and all load currents
    vnames = ["V(in_$(i))" for i in 1:rows]
    vouts = ["V(out_$(j))" for j in 1:cols]
    inames = ["I(VROW_$(i))" for i in 1:rows]
    rnames = ["I(RLOAD_$(j))" for j in 1:cols]
    all_prints = join(vcat(vnames, vouts, inames, rnames), " ")
    println(io, ".PRINT TRAN " * all_prints)
    println(io, ".END")
    return String(take!(io))
end

# Program a simple weight matrix: column 1 matches input pattern [1,0,1], columns 2-3 random
function demo()
    rows, cols = 3, 3
    # Define w in [0,1] (1 => Ron, 0 => Roff)
    w = zeros(rows, cols)
    # Make column 1 a perfect match: rows 1 and 3 at w=1 (LRS), row2 at 0 (HRS)
    w[:,1] = [1.0, 0.0, 1.0]
    # Other columns random mid-states for demonstration
    w[:,2] = [0.2, 0.8, 0.1]
    w[:,3] = [0.0, 1.0, 0.0]

    # Choose Ron/Roff to produce a clear on/off ratio
    Ron = 1e3
    Roff = 100e3
    G = zeros(rows, cols)
    for i in 1:rows, j in 1:cols
        Rval = Ron * w[i,j] + Roff * (1 - w[i,j])
        G[i,j] = 1.0 / Rval
    end

    netlist = build_resistor_crossbar_netlist(G; Vdrive=0.1, rload=10e3, tstep=1e-6, tstop=1e-4)
    sim = create_simulator()
    # Load and run
    loaded = Jyce.loadNetlistString(sim, netlist)
    @assert loaded "Failed to load netlist string"
    prn = joinpath(@__DIR__, "crossbar_resistor_demo.prn")
    result = Jyce.runSimulationOutput(sim, prn)
    @assert Jyce.simulation_success(result) Jyce.simulation_error_message(result)
    prn_file = isfile(prn) ? prn : String(Jyce.simulation_prn_file_path(result))

    # Construct a CrossbarArray object with same G represented as w
    cb = CrossbarArray(rows, cols, MemristorParams(Ron, Roff, 3e-9, 1e-15), w, netlist, prn_file)

    # Read data and generate validation plots
    data = Jyce.read_simulation_data(prn_file)
    res = generate_crossbar_validation_plots(data, cb, prn_file; outdir=joinpath(@__DIR__, "figures"))
    println("Generated plots:", res)
end

# Run demo when executed
if abspath(PROGRAM_FILE) == @__FILE__
    demo()
end
