# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Functions

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
"""
    time        --> time vector
    times2find  --> times to find in time vector (ONLY two elements)
"""
function get_index_from_time(time::Vector, times2find::Vector; tol = 1e-5)
    idx1, idx2 = 1, 2

    for i in eachindex(time)
        error = abs.(time[i] .- times2find) ./ abs.(times2find .+ tol)  # to avoid 1/0

        if error[1] < tol
            idx1 = i
        end
        if error[2] < tol
            idx2 = i
        end
    end
    return (idx1, idx2)
end

@doc """
    plot_pacco: calculates spectrum and plots results from PACCO (given or not the variables to plot)
        experiment      --> experiment's name
        vars2plot       --> vector of variables to plot
        plot_MPT        --> plot MPT?   (Mid-Pleistocene Transition)
        plot_PSD        --> calculate and plot PSD? (Power Spectrum Density)
        plot_proxies    --> include proxy curves?
        proxy_files     --> dictionary with the names of the proxy files to use in T/, co2/ and V/
"""
function plot_pacco(; experiment="test_default", experiments=[], vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_MIS=false, plot_PSD=true, times=(), time_anth=2000.0, plot_proxies=true, proxy_files=Dict("T" => "barker-etal_2011.nc", "co2" => "luthi-etal_2008.nc", "V" => "spratt-lisiecki_2016.nc"))
    # 1. Define some local variables and check if ensemble
    out_path = pwd() * "/output/" * experiment * "/"
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

    isensemble, data_to_load = is_experiment_or_experiments(experiment, experiments)

    colors = [:black, :royalblue4, :red4, :olive, :steelblue4]
    if length(vars2plot) == 1
        palettes = [[colors[1]]]
    else
        color_list = collect(cgrad(colors, length(vars2plot), categorical=true))
        palettes = (color=color_list,)
    end

    if length(data_to_load) == 1
        colors_2 = [:black]
    else
        colors_2 = cgrad(:darktest, length(data_to_load), categorical=true)
    end

    # 2. Load data
    time, data, units = [], [], []          # data[variable][element][values]
    for v in vars2plot
        datav, timev = [], []
        for e in data_to_load
            if e == data_to_load[1]
                units_v = NCDataset(e, "r") do ds
                    ds[v].attrib["units"]
                end # ds is closed
                push!(units, units_v)
            end

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
        end
        push!(time, timev)
        push!(data, datav)
    end

    # 3. Plot
    # -- 3.1 Create figure
    nrows, ncols = length(vars2plot), Int(plot_PSD) + 1
    fgsz = (2000 * ncols, 600 * nrows)
    fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)
    fig = Figure(resolution=fgsz)
    (isensemble) ? (linewidth = 2) : (linewidth = 4)
    for v in eachindex(vars2plot)
        data_v = data[v]
        time_v = time[v]
        var_v = collect_variable(vars2plot[v])

        # -- 3.2 Create axis
        if v == 1
            ax = Axis(fig[v, 1], title="PACCO variables", titlesize=0.8 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[v] * " ($(units[v]))", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz)
        else
            ax = Axis(fig[v, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[v] * " ($(units[v]))", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz)
        end

        # -- 3.3 Create PSD axis if desired
        if plot_PSD
            if v == 1
                ax_PSD = Axis(fig[v, 2], title="Normalized PSD", titlesize=0.8 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period (kyr)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz, xgridvisible = false)
            else
                ax_PSD = Axis(fig[v, 2], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period (kyr)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz, xgridvisible = false)
            end
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
            if var_v in collect(keys(proxy_files))
                proxy_file_v = proxy_path * "/" * var_v * "/" * proxy_files[var_v]
                proxy_ds = NCDataset(proxy_file_v, "r")
                if ~ismissing(proxy_ds[var_v*"_lo"][1])
                    band!(ax, proxy_ds["time"], proxy_ds[var_v*"_lo"], proxy_ds[var_v*"_up"], color=(:grey, 0.4))
                end
                scatter!(ax, proxy_ds["time"], proxy_ds[var_v], label=proxy_files[var_v][1:end-3], color=:black, width=4)
                if plot_PSD
                    d_nomiss = collect(skipmissing(proxy_ds[var_v]))
                    t_nomiss = proxy_ds["time"][broadcast(!, ismissing.(proxy_ds[var_v]))]
                    new_d = Vector{Float64}(undef, length(d_nomiss))
                    new_d[:] = d_nomiss
                    G, f = calc_spectrum(d_nomiss, 1 / (t_nomiss[2] - t_nomiss[1]))
                    G, f = G[f.>1/150e3], f[f.>1/150e3]   # filtering
                    G = G ./ sum(G)
                    periods = 1 ./ f
                    scatterlines!(ax_PSD, periods, G, color=:black, markersize=24, linewidth=4)
                end
            end
        end

        # -- 3.7 Plot each simulation (and PSD if desired)
        lins = []
        for e in eachindex(data_to_load)
            if experiments == []
                if length(data_to_load) == 1
                    color_to_use = palettes[1][v]
                else
                    color_to_use = colors_2[e]
                end
                lbl = "PACCO"
            else
                lbl = experiments[e]
                color_list = collect(cgrad(:darkrainbow, length(experiments), categorical=true))
                color_to_use = (color=color_list,)[1][e]
            end

            if var_v == "T" # CHECK THIS, CUTRE
                l = scatterlines!(ax, time_v[e][:], data_v[e][:] .- 273.15, markersize=4 * linewidth, linewidth=linewidth, color=color_to_use)
            else
                l = scatterlines!(ax, time_v[e][:], data_v[e][:], markersize=4 * linewidth, linewidth=linewidth, color=color_to_use)
            end
            push!(lins, l)

            if plot_PSD
                G, f = calc_spectrum(data_v[e][:], 1 / (time_v[e][2] - time_v[e][1]))
                G, f = G[f.>1/150e3], f[f.>1/150e3]   # filtering
                G = G ./ sum(G)
                periods = 1 ./ f
                scatterlines!(ax_PSD, periods, G, color=color_to_use, markersize=6 * linewidth, linewidth=linewidth)
            end
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

        if (experiments != []) && (v == 1)
            fig[1, ncols+1] = Legend(fig, lins, experiments, framevisible=false, labelsize=0.6 * fntsz)
        end

        if plot_proxies
            if var_v in collect(keys(proxy_files))
                fig[v, ncols+1] = Legend(fig, ax, framevisible=false, labelsize=0.6 * fntsz)
            end
        end

    end

    (plot_PSD) && (colsize!(fig.layout, 2, Relative(1 / 3)))
    resize_to_layout!(fig)

    if experiments == []
        if isensemble
            step = length(data_to_load) / 20
            if step < 1
                step = 1
            else
                step = Int(ceil(step))
            end

            labels = get_exp_labels(data_to_load)
            ticks_selected = 1:step:length(data_to_load)
            new_labels = labels[ticks_selected]

            if length(data_to_load) > 1
                Colorbar(fig[1:end, end+1], height=Relative(1), colormap=:darktest,
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
        save(pacco_path * "/output/" * experiments[1] * "/" * "pacco-nruns_results.png", fig)
    end
end

@doc """
    plot_wavelet: Plots a map of the wavelet scalogram
"""
function plot_wavelet(; experiment="test_default", var2plot="H_n", fs=1 / 1000, sigma=π, MPT=false, time_anth=2000.0)
    ## Load output data 
    df = NCDataset(pacco_path * "/output/" * experiment * "/pacco.nc")
    data, time = df[var2plot], df["time"]

    ## Compute Wavelets
    Wnorm, freqs = calc_wavelet(data, fs; sigma=sigma)
    periods = 1 ./ freqs
    t_coi, p_coi = calc_coi(time, freqs, 2 * pi / sigma)  # convert rad/s to Hz (for cf)
    p_coi = log10.(p_coi ./ 1e3)

    ## Plot
    fig, fntsz = Figure(resolution=(800, 400)), 20

    ax = Axis(fig[1, 1], titlesize=0.8 * fntsz,
        xlabelsize=fntsz, ylabelsize=fntsz, xlabel="Time (kyr)", ylabel="Period (kyr)")
    ax2 = Axis(fig[1, 1], ylabelsize=fntsz, yaxisposition=:right)

    cmap = :cividis
    minW, maxW = minimum(Wnorm), maximum(Wnorm)
    stepW = 0.1 * max(minW, maxW)
    c = contourf!(ax, time, log10.(periods ./ 1e3), Wnorm, colormap=cmap, levels=minW:stepW:maxW)
    (MPT) && (vlines!(ax, [-1.25e6, -0.7e6], linewidth=1, color=:black, linestyle=:dash))    # plot MPT
    (time[end] > time_anth) && (vlines!(ax, [time_anth], linewidth=1, color=:black, linestyle=:dash))
    hlines!(ax, log10.([19, 23, 41, 100]), linewidth=1, color=:black, linestyle=:dash)
    c.extendlow = :auto
    c.extendhigh = :auto
    Colorbar(fig[1, 2], c, height=Relative(1 / 3), width=20, label="Normalized PSD", ticklabelsize=fntsz)

    t1, p1 = t_coi[1:Int(length(p_coi) / 2)], p_coi[1:Int(length(p_coi) / 2)]
    t2, p2 = t_coi[Int(length(p_coi) / 2)+1:end], p_coi[Int(length(p_coi) / 2)+1:end]

    lines!(ax, t1, p1, color=:red, linestyle=:dash)
    lines!(ax, t2, p2, color=:red, linestyle=:dash)
    lines!(ax2, time, data, color=:white, linewidth=1)

    band!(ax, t1, p1, 50, color=("red", 0.3))
    band!(ax, t2, p2, 50, color=("red", 0.3))

    if (time[end] - time[1]) > 1.5e6
        tstep = 300e3
    elseif (time[end] - time[1]) > 1e6
        tstep = 200e3
    else
        tstep = 100e3
    end
    ax.xticks = -5e6:tstep:5e6
    ax.xtickformat = k -> string.(Int.(k / 1000))
    xlims!(ax, (minimum(time), maximum(time)))
    xlims!(ax2, (minimum(time), maximum(time)))

    ax.ytickformat = k -> string.(Int.(ceil.(10 .^ k)))
    ax.yticks = log10.([19, 23, 41, 100])

    ax2.ytickformat = k -> string.(Int.(ceil.(k))) .* " $(df[var2plot].attrib["units"])"
    ax2.yticks = range(minimum(data), stop=maximum(data), length=3)
    ylims!(ax, (1, maximum(p_coi)))
    ylims!(ax2, (minimum(data), 3 * maximum(data)))

    hidespines!(ax2)
    hidexdecorations!(ax2)

    save(pacco_path * "/output/" * experiment * "/" * var2plot * "_wavelet.png", fig)
end

@doc """
    plot_ensemble:
        plots ensemble experiment and variables in vars2plot
"""
function plot_ensemble(; experiment="test_default_ens", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], sims2highlight=[])
    locdir = pwd() * "/output/" * experiment * "/"
    runs = []
    for element in readdir(locdir)
        if element[end-3:end] != ".png"
            push!(runs, element)
        end
    end

    fig = Figure()
    for v in eachindex(vars2plot)
        if v == length(vars2plot)
            ax = Axis(fig[v, 1], ylabel=vars2plot[v], xlabel="Time (kyr)")
        else
            ax = Axis(fig[v, 1], ylabel=vars2plot[v])
        end
        for r in eachindex(runs)
            run = NCDataset(locdir * "/" * runs[r] * "/pacco.nc")
            t, d = run["time"], run[vars2plot[v]]
            if runs[r] in sims2highlight
                lines!(ax, t, d, label=runs[r], linewidth=2, overdraw=true)
            else
                lines!(ax, t, d, color="grey")
            end
            close(run)
        end
        if (sims2highlight != []) && (v == length(vars2plot))
            axislegend(ax, framevisible=false)
        end

        if v != length(vars2plot)
            ax.xticklabelsvisible = false
        end
    end
    save(locdir * "pacco_results.png", fig)
end

# Shortcuts
@doc """
    plot_all: plots H wavelet and pacco main results
"""
function plot_all(; experiment="test_default", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], MPT=false, fs=1 / 1000, sigma=π)
    plot_pacco(experiment=experiment, vars2plot=vars2plot, MPT=MPT)
    plot_wavelet(experiment=experiment, var2plot="H_n", fs=fs, sigma=sigma, MPT=MPT)
    plot_wavelet(experiment=experiment, var2plot="H_n", fs=fs, sigma=sigma, MPT=MPT) # cutre, arreglar -- spm
end

@doc """
"""
function runplot_pacco(; experiment="test_default", par_file="pacco_default.jl", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], MPT=false, fs=1 / 1000, sigma=π)
    run_pacco(experiment=experiment, par_file=par_file)
    plot_all(experiment=experiment, vars2plot=vars2plot, MPT=MPT, fs=fs, sigma=sigma)
end