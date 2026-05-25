"""
    XyceSim.Postprocess

Functions for running Xyce simulations and postprocessing results.
"""
module Postprocess

using Jyce
using Plots
using Statistics
using Printf

export run_crossbar, run_and_plot_crossbar

"""
    run_crossbar(sim::Jyce.XyceSimulator, crossbar::CrossbarArray; prn_path::String="")

Run a crossbar netlist using pure Jyce (native backend).

Writes a stable PRN file (avoids /tmp cleanup issues) and returns the path.
"""
function run_crossbar(sim::Jyce.XyceSimulator, crossbar::CrossbarArray; prn_path::String="")
    # Load netlist into sim (string-based for pure Jyce)
    loaded = Jyce.loadNetlistString(sim, crossbar.netlist)
    @assert loaded "Failed to load crossbar netlist string"

    # Force a stable PRN path
    forced_prn = isempty(prn_path) ? joinpath(@__DIR__, "crossbar_output.prn") : prn_path
    result = Jyce.runSimulationOutput(sim, forced_prn)
    @assert Jyce.simulation_success(result) Jyce.simulation_error_message(result)

    prn_file = isfile(forced_prn) ? forced_prn : String(Jyce.simulation_prn_file_path(result))
    @assert isfile(prn_file) "Data file not found: $(prn_file)"
    return prn_file
end

"""
    run_and_plot_crossbar(
        sim::Jyce.XyceSimulator,
        crossbar::CrossbarArray;
        prn_path::String = "",
        nodes::Vector{String} = ["V(in_1)", "V(out_1)"]
    )

Run a crossbar simulation and generate plots using Jyce utilities.

Returns the PRN path and a NamedTuple of plots.
"""
function run_and_plot_crossbar(
    sim::Jyce.XyceSimulator,
    crossbar::CrossbarArray;
    prn_path::String = "",
    nodes::Vector{String} = ["V(in_1)", "V(out_1)"]
)
    prn_file = run_crossbar(sim, crossbar; prn_path=prn_path)

    data = Jyce.read_simulation_data(prn_file)
    cols = Jyce.detect_iv_columns(data)

    transient_plot = Jyce.plot_transient_voltages(
        prn_file;
        nodes=nodes,
        title="Crossbar Transient"
    )

    iv_plot = Jyce.plot_iv_characteristic(
        prn_file;
        voltage_col=cols.voltage_col,
        current_col=cols.current_col,
        title="Crossbar I-V",
        cycle_overlay=true
    )

    # Additionally compute validation plots: per-column currents and measured vs expected
    try
        plots = generate_crossbar_validation_plots(data, crossbar, prn_file)
    catch e
        @warn "Crossbar validation plotting failed" e
        plots = (currents_plot=nothing, scatter_plot=nothing)
    end

    return (prn_file=prn_file, transient_plot=transient_plot, iv_plot=iv_plot, validation=plots)
end

"""
    generate_crossbar_validation_plots(data, crossbar::CrossbarArray, prn_file::String; outdir::String = joinpath(@__DIR__, "figures"))

Generate validation plots for crossbar simulation results.
"""
function generate_crossbar_validation_plots(data, crossbar::CrossbarArray, prn_file::String; outdir::String = joinpath(@__DIR__, "figures"))
    # Read data frame if necessary
    df = isa(data, String) ? Jyce.read_simulation_data(data) : data

    rows, cols = crossbar.rows, crossbar.cols

    # Build expected column currents from memristor_values and params
    # Conductance model: G = 1 / (Ron * w + Roff * (1-w))
    Ron = crossbar.memristor_params.Ron > 0 ? crossbar.memristor_params.Ron : 1e3
    Roff = crossbar.memristor_params.Roff > 0 ? crossbar.memristor_params.Roff : 100e3
    wmat = crossbar.memristor_values
    G = zeros(rows, cols)
    for i in 1:rows, j in 1:cols
        R = Ron * wmat[i, j] + Roff * (1.0 - wmat[i, j])
        G[i, j] = 1.0 / R
    end

    # Extract time vector
    t = Float64.(df[:, Symbol("TIME")])

    # Extract V_in and V_out arrays over time
    V_in_time = zeros(length(t), rows)
    for i in 1:rows
        name1 = Symbol("V(IN_$(i))")
        name2 = Symbol("V(in_$(i))")
        if hasproperty(df, name1)
            V_in_time[:, i] = Float64.(df[:, name1])
        elseif hasproperty(df, name2)
            V_in_time[:, i] = Float64.(df[:, name2])
        else
            V_in_time[:, i] .= 0.0
        end
    end

    # Choose a representative timestep where the drive amplitude is maximal
    drive_amplitude = vec(sum(abs.(V_in_time), dims=2))
    max_idx = findmax(drive_amplitude)[2]

    # V_in at representative timestep (ensure it's a proper vector)
    V_in = Float64.(vec(Matrix(V_in_time)[max_idx, :]))

    # Extract measured load currents per column over time
    I_loads = zeros(length(t), cols)
    for j in 1:cols
        colname = Symbol("I(RLOAD_$(j))")
        colname2 = Symbol("I(rload_$(j))")
        if hasproperty(df, colname)
            I_loads[:, j] = Float64.(df[:, colname])
        elseif hasproperty(df, colname2)
            I_loads[:, j] = Float64.(df[:, colname2])
        else
            I_loads[:, j] .= 0.0
        end
    end

    # Measured currents at the representative timestep
    I_measured = Float64.(vec(I_loads[max_idx, :]))

    # Compute expected column node voltages accounting for column load resistor (RLOAD)
    # Netlist uses RLOAD_# = 10k; model expected measured load current accordingly.
    Rload = 10e3
    # Sum conductance per column
    col_G_sum = vec(sum(G, dims=1))
    # Currents injected into each column from rows (A)
    raw_col_currents = transpose(G) * V_in
    # Column node voltage (V) = injected_current / (col_G_sum + 1/Rload)
    V_col = raw_col_currents ./ (col_G_sum .+ 1.0 / Rload)
    # Expected measured current through the load resistor = V_col / Rload
    I_expected = Float64.(vec(V_col ./ Rload))

    # Compute read margin: ratio between largest and second-largest expected column current (with epsilon to avoid Inf)
    idxs = sortperm(abs.(I_expected), rev=true)
    max1 = abs(I_expected[idxs[1]])
    max2 = length(idxs) > 1 ? abs(I_expected[idxs[2]]) : 0.0
    eps_margin = 1e-12
    read_margin = max1 / (max2 + eps_margin)

    # Create output directory
    isdir(outdir) || mkpath(outdir)

    # Plot transient currents per column with better visibility
    I_loads_uA = I_loads .* 1e6
    p1 = plot(t .* 1e6, I_loads_uA[:, 1], label = "Col 1", xlabel = "Time (μs)", ylabel = "Current (μA)", 
              title = "Bit-line Currents (Per Column)", linewidth = 2.5, size = (800, 500), legend = :topleft)
    for j in 2:cols
        plot!(p1, t .* 1e6, I_loads_uA[:, j], label = "Col $(j)", linewidth = 2.5)
    end
    max_I = maximum(abs.(I_loads_uA))
    annotate!(p1, maximum(t .* 1e6) * 0.5, max_I * 0.85, 
              text(@sprintf("Read margin = %.2f", read_margin), 10, :black, :left))
    currents_plot = joinpath(outdir, "crossbar_bitline_currents.png")
    savefig(p1, currents_plot)
    println("Saved: ", currents_plot)

    # Scatter measured vs expected with annotations
    # Convert to microamps for readability
    I_expected_uA = I_expected .* 1e6
    I_measured_uA = I_measured .* 1e6
    # Determine plotting limits (tight around data with margin)
    vmin = min(minimum(I_expected_uA), minimum(I_measured_uA))
    vmax = max(maximum(I_expected_uA), maximum(I_measured_uA))
    margin = max(1e-3, (vmax - vmin) * 0.25)
    xmin = vmin - margin
    xmax = vmax + margin

    p2 = scatter(I_expected_uA, I_measured_uA, markersize = 12, label = "Columns",
              xlabel = "Expected I (μA)", ylabel = "Measured I (μA)",
              title = "Measured vs Expected Column Currents", size = (1200, 800),
              legend = :topleft, linewidth = 2, markercolor = :blue, markerstrokecolor = :black, markerstrokewidth = 0.6,
              titlefont = font(18), guidefont = font(14), tickfont = font(12), legendfont = font(12))

    # Identity line across plotting limits
    plot!(p2, [xmin, xmax], [xmin, xmax], lc = :black, linestyle = :dash,
          label = "Ideal (Measured = Expected)", linewidth = 2)

    # Annotate each point with column number and measured value
    for j in 1:cols
      lbl = @sprintf("C%d: %.2f μA", j, I_measured_uA[j])
      annotate!(p2, I_expected_uA[j], I_measured_uA[j] + (vmax - vmin) * 0.02,
            text(lbl, 12, :black, :center))
    end

    # Add read margin annotation in upper-left corner
    annotate!(p2, xmin + (xmax - xmin) * 0.05, xmax - (xmax - xmin) * 0.05,
          text(@sprintf("Read margin = %.2f", read_margin), 14, :black, :left))

    xlims!(p2, xmin, xmax)
    ylims!(p2, xmin, xmax)
    
    scatter_plot = joinpath(outdir, "crossbar_measured_vs_expected.png")
    savefig(p2, scatter_plot)
    println("Saved: ", scatter_plot)

    return (currents_plot=currents_plot, scatter_plot=scatter_plot, read_margin=read_margin)
end

end # module