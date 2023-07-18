# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

"""
    plot_pacco(experiment; vars2plot, plot_MPT, plot_MIS, plot_PSD, times, time_anth, plot_proxies, proxy_files)
calculates spectrum and plots results from PACCO (given or not the variables to plot)

### Arguments
* `experiment` experiment/experiments name/names (`string` or `vector of strings`)

### Optional arguments
* `vars2plot::Vector = ["I", "H", "T", "C", "Vol"]` vector of variables to plot
* `plot_MPT::Bool = false` plot Mid-Pleistocene Transition?
* `plot_MIS::Bool = false` plot Marine Isotope Stages?
* `plot_PSD::Bool = true` calculate and plot Power Spectrum Density?
* `times::Tuple = ()` start and end years (in years)
* `time_anth::Real = 2000.0` when does Anthropocene start?
* `plot_proxies::Bool = true` include proxy curves?
* `proxy_files::Dict = Dict("T" => "barker-etal_2011.nc", "C" => "luthi-etal_2008.nc", "Vol" => "spratt-lisiecki_2016.nc")` dictionary with the names of the proxy files to use in T/, C/ and Vol/
"""
function plot_pacco(experiment; vars2plot::Vector=["I", "H", "T", "C", "Vol"], plot_MPT::Bool=false, plot_MIS::Bool=false, plot_PSD::Bool=true, times::Tuple=(), time_anth::Real=2000.0, plot_proxies::Bool=true, proxy_files::Dict=Dict("T" => "barker-etal_2011.nc", "C" => "luthi-etal_2008.nc", "Vol" => "spratt-lisiecki_2016.nc"))
    # 1. Define some local variables and check if ensemble
    proxy_path = pwd() .* "/data/"

    MIS = Dict(
        "2" => (-12, -115),
        "6" => (-191, -130),
        "8" => (-300, -243),
        "10" => (-374, -337),
        "12" => (-478, -424),
        "14" => (-563, -533),
        "16" => (-676, -621),
        "18" => (-761, -712),
        "20" => (-814, -790),
        "22" => (-900, -866))

    if ~isdir(proxy_path)
        (plot_proxies) && (printstyled("I can't find proxy data directory \n", color=:red))
        plot_proxies = false
    end

    data_to_load, data_labels, experiment = get_runs(experiment)
    multi_colormap = :darkrainbow

    if length(data_to_load) > 1
        pacco_colors = collect(cgrad(multi_colormap, length(data_to_load), categorical=true))
    else
        pacco_colors = [:grey20]
    end

    if plot_proxies == true
        proxy_colors = [:red4, :olive, :royalblue4, :darkorange]
    end

    # 2. Load data
    time, data, units = [], [], []          # data[variable][element][values]
    for v in vars2plot
        datav, timev, unitsv = [], [], []
        for e in data_to_load
            units_ve = NCDataset(e, "r") do ds
                ds[v].attrib["units"]
            end # ds is closed

            t = NCDataset(e, "r") do ds
                ds["time"][:]
            end # ds is closed

            data_ve = NCDataset(e, "r") do ds
                ds[v][:]
            end # ds is closed

            if times != ()
                idx = get_index_from_time(t, [times[1], times[2]])
                push!(timev, t[idx[1]:idx[2]])
                push!(datav, data_ve[idx[1]:idx[2]])
            else
                push!(timev, t)
                push!(datav, data_ve)
            end
            push!(unitsv, units_ve)
        end
        push!(units, unitsv)
        push!(time, timev)
        push!(data, datav)
    end

    # 3. Plot
    # -- 3.1 Create figure
    nrows, ncols = length(vars2plot), Int(plot_PSD) + 1
    fgsz = (600 * nrows * 1.5, 600 * nrows)
    fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)
    fig = Figure(resolution=fgsz)
    (length(data_to_load) > 1) ? (linewidth = 4) : (linewidth = 8)
    k = 0
    for v in eachindex(vars2plot)
        data_v = copy(data[v])
        time_v = copy(time[v])

        # -- 3.2 Create axis
        ax = Axis(fig[v, 1], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz,
            xlabel="Time (kyr)", ylabel=vars2plot[v] * " ($(units[v][1]))",
            xgridcolor=:darkgrey, ygridcolor=:darkgrey,
            xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz)

        # -- 3.3 Create PSD axis if desired
        if plot_PSD
            ax_PSD = Axis(fig[v, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz,
                xlabel="Period (kyr)",
                xgridcolor=:darkgrey, ygridcolor=:darkgrey,
                xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz, xgridvisible=false)
        end

        # -- 3.4 Remove some decorations
        if v < nrows
            hidexdecorations!(ax, grid=false)
            (plot_PSD) && (hidexdecorations!(ax_PSD, grid=false))
        end

        # -- 3.5 Plot MPT and main Milankovitch periods
        (plot_MPT) && (vspan!(ax, -1.25e6, -0.7e6, color=(:red, 0.2)))    # plot Mid-Pleistocene Transition
        (plot_PSD) && (vlines!(ax_PSD, [19e3, 23e3, 41e3, 100e3], linewidth=3, color=:red, linestyle=:dash))  # main Milankovitch periods
        vlines!(ax, [time_anth], linewidth=3, color=(:red, 0.5))


        # -- 3.6 If desired, plot MIS and proxy files
        if plot_MIS
            for m in keys(MIS)
                vspan!(ax, MIS[m][1] * 1e3, MIS[m][2] * 1e3, color=(:lightblue, 0.2))
            end
        end

        if plot_proxies
            if vars2plot[v] in collect(keys(proxy_files))
                k += 1
                proxy_file_v = proxy_path * "/" * vars2plot[v] * "/" * proxy_files[vars2plot[v]]
                proxy_ds = NCDataset(proxy_file_v, "r")
                if ~ismissing(proxy_ds[vars2plot[v]*"_lo"][1])
                    band!(ax, proxy_ds["time"], proxy_ds[vars2plot[v]*"_lo"], proxy_ds[vars2plot[v]*"_up"], color=(proxy_colors[k], 0.4))
                end
                scatter!(ax, proxy_ds["time"], proxy_ds[vars2plot[v]], label=proxy_files[vars2plot[v]][1:end-3], color=proxy_colors[k], width=6)

                if plot_PSD
                    d_nomiss = collect(skipmissing(proxy_ds[vars2plot[v]]))
                    t_nomiss = proxy_ds["time"][broadcast(!, ismissing.(proxy_ds[vars2plot[v]]))]
                    new_d = Vector{Float64}(undef, length(d_nomiss))
                    new_d[:] = d_nomiss
                    G, f = calc_spectrum(d_nomiss, 1 / (t_nomiss[2] - t_nomiss[1]))
                    G, f = G[f.>1/150e3], f[f.>1/150e3]   # filtering
                    G = G ./ sum(G)
                    periods = 1 ./ f
                    lines!(ax_PSD, periods, G, color=proxy_colors[k], markersize=24, linewidth=linewidth)
                end
            end
        end

        # -- 3.7 Plot each simulation (and PSD if desired)
        for e in eachindex(data_to_load)
            label = "PACCO"
            if length(data_to_load) > 1
                label = label * ", $(data_labels[e])"
            end

            if vars2plot[v] in ["T", "Tsurf"]
                lines!(ax, time_v[e][:], data_v[e][:] .- 273.15, markersize=4 * linewidth, linewidth=linewidth, color=pacco_colors[e], label=label)
            else
                lines!(ax, time_v[e][:], data_v[e][:], markersize=4 * linewidth, linewidth=linewidth, color=pacco_colors[e], label=label)
            end

            if plot_PSD
                G, f = calc_spectrum(data_v[e][:], 1 / (time_v[e][2] - time_v[e][1]))
                G, f = G[f.>1/150e3], f[f.>1/150e3]   # filtering
                G = G ./ sum(G)
                periods = 1 ./ f
                lines!(ax_PSD, periods, G, color=pacco_colors[e], linewidth=linewidth)
            end
        end

        if typeof(experiment) == String
            if length(data_to_load) == 1
                fig[v, ncols+1] = Legend(fig, ax, framevisible=false, labelsize=0.6 * fntsz)
            end
        else
            fig[v, ncols+1] = Legend(fig, ax, framevisible=false, labelsize=0.6 * fntsz)
        end

        # -- 3.8 Formatting ...
        if times == ()
            min_time = minimum([time_v[e][1] for e in eachindex(data_to_load)])
            max_time = maximum([time_v[e][end] for e in eachindex(data_to_load)])
            xlims!(ax, (min_time, max_time))
        else
            min_time, max_time = times
            xlims!(ax, times)
        end

        if (max_time - min_time) > 3.0e6
            tstep = 300e3
        elseif (max_time - min_time) > 2.0e6
            tstep = 250e3
        elseif (max_time - min_time) > 1.0e6
            tstep = 200e3
        else
            tstep = 100e3
        end
        ax.xticks = -5e6:tstep:5e6
        ax.xtickformat = k -> string.(Int.(k / 1000))

        if plot_PSD
            ax_PSD.xticks = 0:20e3:150e3
            ax_PSD.xtickformat = k -> string.(Int.(ceil.(k / 1000)))

            hideydecorations!(ax_PSD)
            hidespines!(ax_PSD, :l, :r, :t)
        end
    end

    (plot_PSD) && (colsize!(fig.layout, 2, Relative(1 / 3)))
    resize_to_layout!(fig)

    if typeof(experiment) == String
        if length(data_to_load) > 1
            step = length(data_to_load) / 20
            if step < 1
                step = 1
            else
                step = Int(ceil(step))
            end

            ticks_selected = 1:step:length(data_to_load)
            new_labels = data_labels[ticks_selected]

            if length(data_to_load) > 1
                Colorbar(fig[1:end, end+1], height=Relative(1), colormap=multi_colormap,
                    limits=(1, length(data_to_load)),
                    ticks=(ticks_selected, new_labels),
                    ticklabelsize=0.7 * fntsz
                )
            end

            save(pacco_path * "/output/" * experiment * "/results/" * "pacco_results.png", fig)
        else
            save(pacco_path * "/output/" * experiment * "/" * "pacco_results.png", fig)
        end
    else
        save(pacco_path * "/output/" * experiment[1] * "/" * "pacco-nruns_results.png", fig)
    end
end

"""
    plot_wavelet(experiment; var2plot, fs, sigma, MPT, time_anth)
plots the wavelet map of `experiment`

### Arguments
* `experiment::String` experiment name

### Optional arguments
* `var2plot::String = "H"` variable to plot
* `fs::Real = 1 / 1000` sampling frequency
* `sigma::Real = π` value of central frequency
* `plot_MPT::Bool = false` plot Mid-Pleistocene Transition?
* `time_anth::Real = 2000.0` when does Anthropocene start?
"""
function plot_wavelet(experiment::String; var2plot::String="H", fs::Real=1 / 1000, sigma::Real=π, plot_MPT::Bool=false, time_anth::Real=2000.0)
    ## 1. Load output data 
    path2data = pwd() * "/output/$(experiment)/pacco.nc"
    t = NCDataset(path2data, "r") do ds
        ds["time"][:]
    end # ds is closed

    data = NCDataset(path2data, "r") do ds
        ds[var2plot][:]
    end # ds is closed

    units = NCDataset(path2data, "r") do ds
        ds[var2plot].attrib["units"]
    end # ds is closed

    ## 2. Compute Wavelets
    Wnorm, freqs = calc_wavelet(data, fs; sigma=sigma)
    periods = 1 ./ freqs
    t_coi, p_coi = calc_coi(t, freqs, 2 * pi / sigma)  # convert rad/s to Hz (for cf)
    p_coi = log10.(p_coi ./ 1e3)

    ## 3. Plot
    # 3.1 Define figure and axes
    fig, fntsz = Figure(resolution=(800, 400)), 20
    ax = Axis(fig[1, 1], ylabelsize=fntsz, ylabel="$(var2plot) ($(units))")
    ax2 = Axis(fig[2, 1], titlesize=0.8 * fntsz,
        xlabelsize=fntsz, ylabelsize=fntsz, xlabel="Time (kyr)", ylabel="Period (kyr)")

    # 3.2 Define colors and limits
    cmap = :cividis
    minW, maxW = minimum(Wnorm), maximum(Wnorm)
    stepW = 0.1 * max(minW, maxW)

    # 3.3 Plot variable time series and wavelet map
    lines!(ax, t, data, color=:black, linewidth=1)
    c = contourf!(ax2, t, log10.(periods ./ 1e3), Wnorm, colormap=cmap, levels=minW:stepW:maxW)

    # 3.4 Plot if desired MPT and Anthropocene
    (plot_MPT) && (vlines!(ax, [-1.25e6, -0.7e6], linewidth=1, color=:black, linestyle=:dash))
    (t[end] > time_anth) && (vlines!(ax, [time_anth], linewidth=1, color=:black, linestyle=:dash))
    (plot_MPT) && (vlines!(ax2, [-1.25e6, -0.7e6], linewidth=1, color=:black, linestyle=:dash))
    (t[end] > time_anth) && (vlines!(ax2, [time_anth], linewidth=1, color=:black, linestyle=:dash))

    # 3.5 Plot main Milankovitch periodicities
    hlines!(ax2, log10.([19, 23, 41, 100]), linewidth=1, color=:black, linestyle=:dash)

    # 3.6 Create wavelet colorbar
    c.extendlow = :auto
    c.extendhigh = :auto
    Colorbar(fig[2, 2], c, height=Relative(2 / 3), width=20, label="Normalized PSD", ticklabelsize=fntsz)

    # 3.7 Plot COI
    t1, p1 = t_coi[1:Int(length(p_coi) / 2)], p_coi[1:Int(length(p_coi) / 2)]
    t2, p2 = t_coi[Int(length(p_coi) / 2)+1:end], p_coi[Int(length(p_coi) / 2)+1:end]

    lines!(ax2, t1, p1, color=:red, linestyle=:dash)
    lines!(ax2, t2, p2, color=:red, linestyle=:dash)
    band!(ax2, t1, p1, 50, color=("red", 0.3))
    vspan!(ax2, t[1], t1[end], color=("red", 0.3))
    band!(ax2, t2, p2, 50, color=("red", 0.3))
    vspan!(ax2, t2[end], t[end], color=("red", 0.3))

    # 3.8 Adjust figure
    if (t[end] - t[1]) > 2.0e6
        tstep = 300e3
    elseif (t[end] - t[1]) > 1e6
        tstep = 200e3
    else
        tstep = 100e3
    end
    ax.xticks = -5e6:tstep:5e6
    ax.xtickformat = k -> string.(Int.(k / 1000))
    ax2.xticks = -5e6:tstep:5e6
    ax2.xtickformat = k -> string.(Int.(k / 1000))

    xlims!(ax, (minimum(t), maximum(t)))
    xlims!(ax2, (minimum(t), maximum(t)))

    ax2.ytickformat = k -> string.(Int.(ceil.(10 .^ k)))
    ax2.yticks = log10.([19, 23, 41, 100])

    ylims!(ax2, (1, maximum(p_coi)))

    hidexdecorations!(ax, grid=false)
    rowsize!(fig.layout, 1, Relative(1 / 3))
    rowgap!(fig.layout, 5.0)

    save(pacco_path * "/output/" * experiment * "/" * var2plot * "_wavelet.png", fig)
end

"""
    plot_states(experiment)
plots all the states of PACCO `experiment`

### Arguments
* `experiment::String` experiment name
"""
function plot_pacco_states(experiment)
    prognostic = ["T", "C", "iceage", "albedo", "H", "Hsed", "B", "Tice", "fstr"]
    diagnostic = ["I", "Tsl", "Tref",  # radiative forcing and climate response
    "z", "Surf", "Vol",  # ice geometry
    "albedo_ref", "Tsurf",   # climate parameters
    "s", "a",   # ice-sheet mass balance
    "taud", "taub", "vd", "vb", "fstr_ref", "v", # ice dynamics
    "Qdif", "Qdrag"]   # ice thermodynamics
    pacco_states = vcat(prognostic, diagnostic)

    path2data = pwd() * "/output/$(experiment)/pacco.nc"
    data_frame = NCDataset(path2data, "r")
    fig, k = Figure(resolution=(2000, 800)), 1
    for i in 1:4, j in 1:6
        ax = Axis(fig[i, j], ylabel=pacco_states[k])
        lines!(ax, data_frame["time"] ./ 1e3, data_frame[pacco_states[k]])
        k += 1
    end

    save(pwd() * "/output/" * experiment * "/pacco_states.png", fig)
end

"""
    plot_states(experiment)
makes some composite plots of PACCO `experiment`

### Arguments
* `experiment::String` experiment name
"""
function plot_pacco_comp_states(experiment::String)
    path2data = pwd() * "/output/$(experiment)/pacco.nc"
    parameters2use = JLD2.load_object(pwd() * "/output/$(experiment)/params.jld2")
    data_frame = NCDataset(path2data, "r")

    selected_albedo = Vector{Float32}(undef, length(data_frame["time"]))
    for t in eachindex(data_frame["time"])
        if data_frame["MB"][t] <= 0
            selected_albedo[t] = data_frame["albedo"][t]
        else
            selected_albedo[t] = parameters2use.albedo_newice
        end
    end

    # Temperature evolution
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="Temperature evolution")
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* parameters2use.cZ .* data_frame["Z"] ./ parameters2use.tau_T, label="-cz⋅Z/tau", color=:royalblue4)
    lines!(ax, data_frame["time"] ./ 1e3, (data_frame["Tref"] .- data_frame["T"]) ./ parameters2use.tau_T, color=:olive, label="(Tref-T)/tau")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.cI .* data_frame["Ianom"] ./ parameters2use.tau_T, label="cᵢ⋅Ianom", color=:orange)
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.cC .* 5.35 .* NaNMath.log.(data_frame["C"] ./ 280.0) ./ parameters2use.tau_T, label="cc⋅5.35⋅log(C/280)", color=:maroon)
    lines!(ax, data_frame["time"] ./ 1e3, (data_frame["Tref"] .+ parameters2use.cI .* data_frame["Ianom"] .+ parameters2use.cC .* 5.35 .* NaNMath.log.(data_frame["C"] ./ 280.0) - parameters2use.cZ .* data_frame["Z"] - data_frame["T"]) ./ parameters2use.tau_T, label="dT/dt", color=:grey20)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[2, 1], ylabel="T")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["T"] .- data_frame["Tref"], color=:grey20, label="T")
    save(pwd() * "/output/" * experiment * "/pacco_comp_states_T.png", fig)

    # Ice evolution
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="Ice evolution")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["A"], label="A")
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* data_frame["M"], label="M")
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* data_frame["U"] .* data_frame["H"] ./ parameters2use.L, label="U⋅H/L", color=:green)
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["A"] .- data_frame["M"] .- data_frame["U"] .* data_frame["H"] ./ parameters2use.L, label="dH/dt", color=:grey20)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[2, 1], ylabel="Ice velocity", yaxisposition=:right, yticklabelcolor=:green, ygridvisible=false, ylabelcolor=:green)
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["U"], color=:green)
    ax = Axis(fig[2, 1], ylabel="Size (m)")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["H"], color=:grey20, label="H")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["Z"], color=:royalblue4, label="Z")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["Beq"], color=:violet, label="Beq")
    fig[2, 2] = Legend(fig, ax, framevisible=false)
    save(pwd() * "/output/" * experiment * "/pacco_comp_states_H.png", fig)

    # Melting terms
    fig = Figure(resolution=(1000, 500))
    ax = Axis(fig[1, 1], ylabel="Melting terms")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["MB"], color=:grey, label="MB")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["A"], label="A", color="blue")
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* data_frame["M"], label="M", color="red")
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* max.(parameters2use.kI .* (1 .- selected_albedo) .* data_frame["Ianom"], 0.0), label="-kᵢ⋅(1-albedo)⋅Ianom", color=:orange)
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* max.(parameters2use.lambda .* (data_frame["T"] .- data_frame["Tref"]), 0.0), label="-lambda⋅(T-Tref)", color=:maroon)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    save(pwd() * "/output/" * experiment * "/pacco_comp_states_MB.png", fig)

    # Albedo
    fig = Figure(resolution=(1000, 500))
    ax = Axis(fig[1, 1], ylabel="Albedo", yticklabelcolor=:green, ylabelcolor=:green)
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["albedo"], label="albedo", color=:olive, markersize=5)
    scatter!(ax, data_frame["time"] ./ 1e3, selected_albedo, label="albedo selected", color=:green, markersize=5)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[1, 1], ylabel="iceage", yaxisposition=:right, ylabelcolor=:grey20)
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["iceage"], color=:grey20)
    save(pwd() * "/output/" * experiment * "/pacco_comp_states_albedo.png", fig)
end

"""
    fastplot(experiment, y; x="time", use_colormap=false, plot_function=lines!)
makes some composite plots of PACCO `experiment`

### Arguments
* `experiment` experiment/experiments name/names (`string` or `vector of strings`)
* `y::String` variable to plot in y axis
* `x::String = "time"` variable to plot in x axis
* `use_colormap::Bool = false` use colormap? (use only if one simulation)
* `plot_function::Function = lines!` CairoMakie function to use when plotting (lines!, scatter!) 
"""
function fastplot(experiment, y::String; x::String="time", use_colormap::Bool=false, plot_function::Function=lines!)
    data_to_load, data_labels, experiment = get_runs(experiment)

    fig = Figure()
    ax = Axis(fig[1, 1])
    times = []
    for i in eachindex(data_to_load)
        df = NCDataset(data_to_load[i])
        if use_colormap
            times = df["time"]
            plot_function(ax, df[x], df[y], label=data_labels[i], color=df["time"], colormap=:blues)
        else
            plot_function(ax, df[x], df[y], label=data_labels[i])
        end
        if i == 1
            if x == "time"
                ax.xlabel = "Time (yr)"
            else
                ax.xlabel = "$(x) ($(df[x].attrib["units"]))"
            end
            ax.ylabel = "$(y) ($(df[y].attrib["units"]))"
        end
    end

    fig[1, 2] = Legend(fig, ax, framevisible=false)

    if use_colormap
        Colorbar(fig[2, 1], colormap=:blues, limits=(times[1], times[end]), vertical=false, label="Time (yr)")
    end

    if typeof(experiment) == String
        if length(data_to_load) > 1
            save(pacco_path * "/output/" * experiment * "/results/" * "pacco_fastplot_$(y)vs$(x).png", fig)
        else
            save(pacco_path * "/output/" * experiment * "/" * "pacco_fastplot_$(y)vs$(x).png", fig)
        end
    else
        save(pacco_path * "/output/" * experiment[1] * "/" * "pacco_nruns_fastplot_$(y)vs$(x).png", fig)
    end
end

function plot_recurrence(experiment; var2plot::String="H")

end