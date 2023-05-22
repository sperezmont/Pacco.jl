# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

"""
    plot_pacco(experiment; vars2plot, plot_MPT, plot_MIS, plot_PSD, times, time_anth, plot_proxies, proxy_files)
calculates spectrum and plots results from PACCO (given or not the variables to plot)

## Attributes
* `experiment`      --> experiment/experiments name/names
* `vars2plot`       --> vector of variables to plot
* `plot_MPT`        --> plot MPT?   (Mid-Pleistocene Transition)
* `plot_PSD`        --> calculate and plot PSD? (Power Spectrum Density)
* `plot_proxies`    --> include proxy curves?
* `proxy_files`     --> dictionary with the names of the proxy files to use in T/, co2/ and V/
"""
function plot_pacco(experiment; vars2plot=["I", "H", "T", "co2", "V"], plot_MPT=false, plot_MIS=false, plot_PSD=true, times=(), time_anth=2000.0, plot_proxies=true, proxy_files=Dict("T" => "barker-etal_2011.nc", "co2" => "luthi-etal_2008.nc", "V" => "spratt-lisiecki_2016.nc"))
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

    data_to_load, data_labels = get_runs(experiment)
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
    fgsz = (2000 * ncols, 600 * nrows)
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
                          xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz, xgridvisible = false)
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
                vspan!(ax, MIS[m][1] * 1e3, MIS[m][2] * 1e3, color = (:lightblue, 0.2))
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

            if vars2plot[v] == "T"
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

        if (max_time - min_time) > 2.0e6
            tstep = 300e3
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
"""
function plot_wavelet(experiment; var2plot="H", fs=1 / 1000, sigma=π, MPT=false, time_anth=2000.0)
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
    (MPT) && (vlines!(ax, [-1.25e6, -0.7e6], linewidth=1, color=:black, linestyle=:dash)) 
    (t[end] > time_anth) && (vlines!(ax, [time_anth], linewidth=1, color=:black, linestyle=:dash))
    (MPT) && (vlines!(ax2, [-1.25e6, -0.7e6], linewidth=1, color=:black, linestyle=:dash)) 
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

function plot_pacco_states(experiment)
    prognostic = ["T", "co2", "iceage", "alpha", "H", "Hsed", "B", "Tice", "fstream"]
    diagnostic = ["I", "R", "Tsl", "Tref",  # radiative forcing and climate response
        "Z", "E", "V",  # ice geometry
        "alpha_ref", "Tsurf",   # climate parameters
        "A", "M",   # ice-sheet mass balance
        "taud", "taub", "Ud", "Ub", "fstream_ref", "U", # ice dynamics
        "Qdif", "Qdrag"]   # ice thermodynamics
    pacco_states = vcat(prognostic, diagnostic)

    path2data = pwd() * "/output/$(experiment)/pacco.nc"
    data_frame = NCDataset(path2data, "r")
    fig, k = Figure(resolution=(2000, 800)), 1
    for i in 1:4, j in 1:7
        ax = Axis(fig[i, j], ylabel=pacco_states[k])
        lines!(ax, data_frame["time"] ./ 1e3, data_frame[pacco_states[k]])
        k += 1
    end

    save(pwd()*"/output/"*experiment*"/pacco_states.png", fig)
end

function plot_pacco_comp_states(experiment::String)
    path2data = pwd() * "/output/$(experiment)/pacco.nc"
    parameters2use = JLD2.load_object(pwd() * "/output/$(experiment)/params.jld2")
    data_frame = NCDataset(path2data, "r")

    selected_alpha = Vector{Float32}(undef, length(data_frame["time"]))
    for t in eachindex(data_frame["time"])
        if data_frame["MB"][t] <= 0
            selected_alpha[t] = data_frame["alpha"][t]
        else
            selected_alpha[t] = parameters2use.alpha_newice
        end
    end
    
    # Temperature evolution
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="Temperature evolution")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.ci .* data_frame["Ianom"] ./ parameters2use.tau_T, label="Insolation forcing")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.cc .* 5.35 .* NaNMath.log.(data_frame["co2"] ./ 280.0) ./ parameters2use.tau_T, label="CO2 feedback")
    lines!(ax, data_frame["time"] ./ 1e3, -1 .* parameters2use.cz .* data_frame["Z"] ./ parameters2use.tau_T, label="Cooling effect")
    lines!(ax, data_frame["time"] ./ 1e3, (data_frame["Tref"] .+ parameters2use.ci .* data_frame["Ianom"] .+ parameters2use.cc .* 5.35 .* NaNMath.log.(data_frame["co2"] ./ 280.0) - parameters2use.cz .* data_frame["Z"] - data_frame["T"]) ./ parameters2use.tau_T, label="dT/dt", color=:grey20)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[2, 1], ylabel="Temperature")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["T"], color=:grey20, label="T")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["Tref"], color=:grey, label="Tref")
    fig[2, 2] = Legend(fig, ax, framevisible=false)
    save(pwd()*"/output/"*experiment*"/pacco_comp_states_T.png", fig)

    # Ice evolution
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="Ice evolution")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["A"], label="A")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["M"], label="M")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["U"] .* data_frame["H"] ./ parameters2use.L, label="U⋅H/L", color=:green)
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["A"] .- data_frame["M"] .- data_frame["U"] .* data_frame["H"] ./ parameters2use.L, label="dH/dt", color=:grey20)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[2, 1], ylabel="Size (m)")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["H"], color=:grey20, label="H")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["Z"], color=:royalblue4, label="Z")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["B"], color=:violet, label="B")
    fig[2, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[2, 1], ylabel="Ice velocity", yaxisposition=:right, yticklabelcolor=:green, ygridvisible=false, ylabelcolor=:green)
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["U"], color=:green)
    save(pwd()*"/output/"*experiment*"/pacco_comp_states_H.png", fig)

    # Melting terms
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="Melting terms")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.ki .* (1 .- data_frame["alpha"]) .* data_frame["Ianom"], label="Insolation term, MB ≤ 0")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.ki .* (1 .- parameters2use.alpha_newice) .* data_frame["Ianom"], label="Insolation term, MB > 0")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.ki .* (1 .- selected_alpha) .* data_frame["Ianom"], label="Insolation term, α selected")
    lines!(ax, data_frame["time"] ./ 1e3, parameters2use.lambda .* (data_frame["T"] .- data_frame["Tref"]), label="Temperature term")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["A"], label="A", color="navy")
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["M"], label="M", color="red")
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[2, 1], ylabel="Mass balance")
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["MB"], color=:grey20)
    save(pwd()*"/output/"*experiment*"/pacco_comp_states_MB.png", fig)

    # Albedo
    fig = Figure(resolution=(1000, 500))
    ax = Axis(fig[1, 1], ylabel="Albedo", yticklabelcolor=:green, ylabelcolor=:green)
    scatter!(ax, data_frame["time"] ./ 1e3, data_frame["alpha"], label="α", color=:olive, markersize=5)
    scatter!(ax, data_frame["time"] ./ 1e3, selected_alpha, label="α selected", color=:green, markersize=5)
    fig[1, 2] = Legend(fig, ax, framevisible=false)
    ax = Axis(fig[1, 1], ylabel="MB", yaxisposition=:right)
    barplot!(ax, data_frame["time"] ./ 1e3, data_frame["MB"], color=:grey)
    if minimum(data_frame["MB"]) != -minimum(data_frame["MB"])
        ylims!(ax, minimum(data_frame["MB"]), -minimum(data_frame["MB"]))
    end
    ax = Axis(fig[1, 1], ylabel="H", yaxisposition=:right, ylabelcolor=:royalblue4, yticklabelcolor=:royalblue4)
    lines!(ax, data_frame["time"] ./ 1e3, data_frame["H"], color=:royalblue4)
    ylims!(ax, -maximum(data_frame["H"]), maximum(data_frame["H"]))
    save(pwd()*"/output/"*experiment*"/pacco_comp_states_alpha.png", fig)
end


# OLD CODE, do not use these functions!

function collect_variable(str::String)
    if str[1:end-2] in ["T_surf", "T_ice", "T_rf", "T_ref"]
        return str[1:end-2]
    else
        idx = 1
        for i in eachindex(str)
            if str[i] == '_'
                idx = i
                break
            end
        end
        return str[1:idx-1]
    end
end

"""
    plot_std(experiment, filename)

"""
function plot_std(experiment; var2plot="H_n", cmap=:darktest)
    path_to_results = get_path_to_results_or_runs(experiment, "runs")
    runs = path_to_results .* "/" .* readdir(path_to_results) .* "/pacco.nc"
    labels = get_exp_labels(runs)
    # if length(runs) > 1000
    #     error("Too much runs to plot ... (> 1000)")
    # end

    fig = Figure()
    ax = Axis(fig[1, 1], ylabel="std $(var2plot)")
    color_list = collect(cgrad(cmap, length(runs), categorical=true))
    palettes, k = (color=color_list,)[1], 1
    for run in runs
        df = NCDataset(run, "r")
        x = df[var2plot][:]

        scatter!(ax, k, std(x), color=palettes[k], markersize=30)

        k += 1
    end

    step = length(runs) / 20
    if step < 1
        step = 1
    else
        step = Int(ceil(step))
    end

    ticks_selected = 1:step:length(runs)
    new_labels = labels[ticks_selected]

    ax.xticks = (ticks_selected, new_labels)
    Colorbar(fig[1, 2], height=Relative(1), colormap=cmap,
        limits=(1, length(runs)),
        ticks=(ticks_selected, new_labels),
        label="runs"
    )
    save(get_path_to_results_or_runs(experiment, "results") * "/ensemble_std.png", fig)
end


"""
    fast_histogram(experiment, filename)
this function assumes we want a no so fancy plot so it just plots the histogram of an ensembleof good runs stored in `good_runs.txt`

"""
function fast_histogram(experiment, filename; all_runs=true, plot_name="fast_histogram.png")
    path_to_results = get_path_to_results_or_runs(experiment, "results")

    if all_runs
        runs = readlines(path_to_results * "/$(filename)") .* "/pacco.nc"
    else
        runs = get_good_runs_from_file(path_to_results * "/$(filename)") .* "/pacco.nc"
    end
    labels = get_exp_labels(runs)

    aux_runs_vector = readlines(path_to_results * "/$(filename)")
    par_names, parameters = get_parameters_from_runs(experiment, filename, all_runs=all_runs)

    fig = Figure(resolution=(500 * length(par_names), 500))

    for i in eachindex(par_names)
        if i == 1
            ax = Axis(fig[1, i], title=par_names[i], xlabel="Values", ylabel="Amount of runs", xticklabelrotation=pi / 2)
        else
            ax = Axis(fig[1, i], title=par_names[i], xlabel="Values", xticklabelrotation=pi / 2)
        end
        ax2 = Axis(fig[1, i], yaxisposition=:right, ygridstyle=:dash, ygridcolor=(:red, 0.4), rightspinecolor=:red, yticklabelcolor=:red, xticklabelrotation=pi / 2)

        data = [parameters[j][i] for j in 1:length(parameters)]

        color = 1:length(data)
        if all_runs == true
            color = 1:length(aux_runs_vector)
            for d in eachindex(data)
                if aux_runs_vector[d][1:3] == "NS_"
                    data[d] = NaN64
                end
            end
        end

        ocurrences = counter(data)
        display(ocurrences)
        kys, vls = collect(keys(ocurrences)), collect(values(ocurrences))

        barplot!(ax, kys, vls, color=(:black, 0.5))
        scatter!(ax2, data, 1:length(data), colormap=:darktest, color=color)

        step = maximum(vls) / 20
        if step < 1
            step = 1
        else
            step = Int(ceil(step))
        end

        ticks_selected = 0:step:maximum(vls)
        ax2.yticks = ticks_selected

        step = length(runs) / 20
        if step < 1
            step = 1
        else
            step = Int(ceil(step))
        end

        ticks_selected = 1:step:length(runs)
        new_labels = labels[ticks_selected]

        ax.xticks, ax2.xticks = kys, kys
        ax2.yticks = (ticks_selected, new_labels)
        ylims!(ax, (0, length(runs) + 5))
        ylims!(ax2, (0, length(runs) + 5))
        linkxaxes!(ax, ax2)

        if i != length(par_names)
            hideydecorations!(ax2, grid=false)
        end
        if i != 1
            hideydecorations!(ax, grid=false)
        end

    end

    save(path_to_results * "/$(plot_name)", fig)
end


"""
    fast_plot(experiment, filename; var2plot="H_n", cmap=:heat)
this function assumes we want a no so fancy plot so it just plots variable `var2plot` of good runs stored in `good_runs.txt` with a gradation of colors

"""
function fast_plot(experiment, filename; var2plot="H_n", cmap=:darkrainbow, all_runs=true, plot_PSD=false, plot_name="fast_plot.png")
    path_to_results = get_path_to_results_or_runs(experiment, "results")

    if all_runs
        runs = readlines(path_to_results * "/$(filename)") .* "/pacco.nc"
    else
        runs = get_good_runs_from_file(path_to_results * "/$(filename)") .* "/pacco.nc"
    end
    labels = get_exp_labels(runs)

    # if length(runs) > 1000
    #     error("Too much runs to plot ... (> 1000)")
    # end

    fig = Figure(resolution=(2000, 600))
    (plot_PSD) ? (xlabel = "Period (kyr)") : (xlabel = "Time (kyr)")
    ax = Axis(fig[1, 1], ylabel=var2plot, xlabel=xlabel)
    if length(runs) > 1
        color_list = collect(cgrad(cmap, length(runs), categorical=true))
        palettes, k = (color=color_list,)[1], 1
    else
        palettes, k = :red, 1
    end
    t, periods = [], []
    for run in runs
        if run[1:3] == "NS_"
            k += 1
            continue
        end

        df = NCDataset(run, "r")
        x, t = df[var2plot][:], df["time"][:]

        if length(runs) > 1
            color = palettes[k]
        else
            color = palettes
        end

        if plot_PSD == true
            G, f = calc_spectrum(x, 1 / (t[2] - t[1]))
            time2cut = t[Int(ceil(length(t) / 2))]
            G, f = G[f.>-1/time2cut], f[f.>-1/time2cut]   # filtering
            G = G ./ sum(G)
            periods = 1 ./ f
            scatterlines!(ax, periods ./ 1000, G, color=color, markersize=3, linewidth=0.5)
        else
            scatterlines!(ax, t ./ 1000, x, color=color, markersize=3, linewidth=0.5)
        end
        k += 1
    end

    if plot_PSD == false
        xlims!(ax, (t[1] / 1000, t[end] / 1000))
        if (t[end] - t[1]) > 2.0e3
            tstep = 300
        elseif (t[end] - t[1]) > 800
            tstep = 200
        else
            tstep = 100
        end
        ax.xticks = -5e3:tstep:5e3
        ax.xtickformat = k -> string.(Int.(k))
    else
        ax.xticks = 0:20:500
    end

    step = length(runs) / 20
    if step < 1
        step = 1
    else
        step = Int(ceil(step))
    end

    ticks_selected = 1:step:length(runs)
    new_labels = labels[ticks_selected]

    if length(runs) > 1
        Colorbar(fig[1, 2], height=Relative(1), colormap=cmap,
            limits=(1, length(runs)),
            ticks=(ticks_selected, new_labels),
            label="runs"
        )
    end
    if plot_PSD
        save(path_to_results * "/$(var2plot)_PSD_$(plot_name)", fig)
    else
        save(path_to_results * "/$(var2plot)_$(plot_name)", fig)
    end
end

