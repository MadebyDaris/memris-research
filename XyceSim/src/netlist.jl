"""
    XyceSim.Netlist

Functions for creating and writing Xyce netlists for memristor crossbar arrays.
"""
module Netlist

using Printf
using Random
using LinearAlgebra

export MemristorParams, CrossbarArray
export create_crossbar_array, create_crossbar_array_subckt
export write_netlist

"""
    MemristorParams

Parameters for memristor device model.
"""
Base.@kwdef struct MemristorParams
    Ron::Float64 = 1000.0      # ON resistance (Ω)
    Roff::Float64 = 10000.0    # OFF resistance (Ω)
    D::Float64 = 3e-9          # Diffusion coefficient (cm²/s)
    uv::Float64 = 1e-15        # Ion mobility (cm²/V/s)
end

"""
    CrossbarArray

Crossbar array structure for Xyce simulation.
"""
struct CrossbarArray
    rows::Int
    cols::Int
    memristor_params::MemristorParams
    memristor_values::Matrix{Float64}
    netlist::String
    prn_file::String
end

"""
    create_crossbar_array(
        rows::Int,
        cols::Int;
        file_name::String = "",
        params::MemristorParams = MemristorParams(),
        analysis::Symbol = :tran,
        tstep::Float64 = 1e-3,
        tstop::Float64 = 1e-1,
        rng::AbstractRNG = Random.default_rng(),
    )

Build a crossbar netlist for the ADMS memristor plugin.

The generated netlist includes:
- one DC source per row,
- one load resistor per column (to avoid floating nodes),
- one memristor per cross-point, with a dedicated internal state node.
"""
function create_crossbar_array(
    rows::Int,
    cols::Int;
    file_name::String = "",
    params::MemristorParams = MemristorParams(),
    analysis::Symbol = :tran,
    tstep::Float64 = 1e-3,
    tstop::Float64 = 1e-1,
    rng::AbstractRNG = Random.default_rng(),
)
    @assert rows > 0 "rows must be > 0"
    @assert cols > 0 "cols must be > 0"

    memristor_values = zeros(rows, cols)
    for i in 1:rows
        for j in 1:cols
            memristor_values[i, j] = rand(rng)
        end
    end

    io = IOBuffer()
    println(io, "* Crossbar array netlist")
    println(io, ".OPTIONS DEVICE level=2")
    println(io)

    for i in 1:rows
        println(io, "VROW_$(i) in_$(i) 0 DC 0.1")
    end
    for j in 1:cols
        println(io, "RLOAD_$(j) out_$(j) 0 10k")
    end
    println(io)

    println(io, ".MODEL memristor_model Memristor level=1")
    println(io, "+ model=4 window_type=5 dt=1ms")
    println(io, "+ ron=$(params.Ron) roff=$(params.Roff) D=$(params.D) uv=$(params.uv)")
    println(io)

    for i in 1:rows
        for j in 1:cols
            # Xyce ADMS plugin syntax: Y<device_type> <instance_name> <nodes...> <model_name>
            println(io, "YMEMRISTOR M_$(i)_$(j) in_$(i) out_$(j) w_$(i)_$(j) memristor_model")
        end
    end
    println(io)

    if analysis == :tran
        println(io, ".TRAN $(tstep) $(tstop)")
        println(io, ".PRINT TRAN V(in_1) V(out_1) I(VROW_1)")
    else
        println(io, ".OP")
        println(io, ".PRINT OP V(in_1) V(out_1) I(VROW_1)")
    end
    println(io, ".END")

    netlist = String(take!(io))
    
    prn_file = isempty(file_name) ? "" : file_name
    return CrossbarArray(rows, cols, params, memristor_values, netlist, prn_file)
end

"""
    create_crossbar_array_subckt(
        rows::Int,
        cols::Int;
        file_name::String = "",
        subckt_path::String = "biolek_mmrstor_model.sub",
        row_source::String = "SIN(0 0.1 1000)",
        analysis::Symbol = :tran,
        tstep::Float64 = 1e-3,
        tstop::Float64 = 1e-1,
        rng::AbstractRNG = Random.default_rng(),
        params::MemristorParams = MemristorParams(1e3, 100e3, 3e-9, 1e-15),
    )

Build a crossbar netlist using a .subckt memristor model (no plugin required).

Uses XMEM instances with dedicated state nodes and supports transient or OP analysis.
"""
function create_crossbar_array_subckt(
    rows::Int,
    cols::Int;
    file_name::String = "",
    subckt_path::String = "biolek_mmrstor_model.sub",
    row_source::String = "SIN(0 0.1 1000)",
    analysis::Symbol = :tran,
    tstep::Float64 = 1e-3,
    tstop::Float64 = 1e-1,
    rng::AbstractRNG = Random.default_rng(),
    params::MemristorParams = MemristorParams(1e3, 100e3, 3e-9, 1e-15),
)
    @assert rows > 0 "rows must be > 0"
    @assert cols > 0 "cols must be > 0"

    sub_candidates = [subckt_path, joinpath(@__DIR__, subckt_path), joinpath(@__DIR__, "biolek_mmrstor_model.sub")]
    sub_idx = findfirst(isfile, sub_candidates)
    sub_idx === nothing && error("Subckt not found: " * subckt_path)
    sub = abspath(sub_candidates[sub_idx])

    memristor_values = zeros(rows, cols)
    for i in 1:rows
        for j in 1:cols
            memristor_values[i, j] = rand(rng)
        end
    end

    io = IOBuffer()
    println(io, "* Crossbar array netlist (subckt)")
    println(io, ".INCLUDE \"$(sub)\"")
    println(io)

    # Row voltage sources
    for i in 1:rows
        println(io, "VROW_$(i) in_$(i) 0 $(row_source)")
    end
    # Column load resistors
    for j in 1:cols
        println(io, "RLOAD_$(j) out_$(j) 0 10k")
    end
    println(io)

    # Instantiate memristor subcircuits
    for i in 1:rows
        for j in 1:cols
            println(io, "XMEM_$(i)_$(j) in_$(i) out_$(j) w_$(i)_$(j) MEM_BIOLEK")
        end
    end
    println(io)

    # Build a comprehensive .PRINT line that records all row/col voltages and load currents
    if analysis == :tran
        println(io, ".TRAN $(tstep) $(tstop)")
        # Voltages
        vnames = ["V(IN_$(i))" for i in 1:rows]
        vnames_out = ["V(OUT_$(j))" for j in 1:cols]
        # Currents through row sources and load resistors
        inames = ["I(VROW_$(i))" for i in 1:rows]
        rnames = ["I(RLOAD_$(j))" for j in 1:cols]

        all_prints = join(vcat(vnames, vnames_out, inames, rnames), " ")
        println(io, ".PRINT TRAN " * all_prints)
    else
        println(io, ".OP")
        println(io, ".PRINT OP V(in_1) V(out_1) I(VROW_1)")
    end
    println(io, ".END")

    netlist = String(take!(io))
    prn_file = isempty(file_name) ? "" : file_name
    return CrossbarArray(rows, cols, params, memristor_values, netlist, prn_file)
end

"""
    write_netlist(crossbar::CrossbarArray, output_path::String)

Write the generated netlist to disk.
"""
function write_netlist(crossbar::CrossbarArray, output_path::String)
    open(output_path, "w") do io
        write(io, crossbar.netlist)
    end
    return output_path
end

end # module