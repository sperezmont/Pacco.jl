# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Functions

function is_ensemble(exp)
    elems = readdir(pwd() * "/output/" * exp * "/")
    if "pacco.nc" in elems
        is_ens = false
        new_elems = elems
    else
        is_ens = true

        new_elems = []
        for e in elems # drop "/results/" directory
            if e != "results"
                push!(new_elems, e)
            end
        end
    end
    return new_elems, is_ens
end

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

@doc """
    plot_pacco: calculates spectrum and plots results from PACCO (given or not the variables to plot)
        experiment      --> experiment's name
        vars2plot       --> vector of variables to plot
        plot_MPT        --> plot MPT?   (Mid-Pleistocene Transition)
        plot_PSD        --> calculate and plot PSD? (Power Spectrum Density)
        plot_proxies    --> include proxy curves?
        proxy_files     --> dictionary with the names of the proxy files to use in T/, co2/ and V/
"""
function plot_pacco(; experiment="test_default", experiments=[], vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true, time_anth=2000.0, plot_proxies=true, proxy_files=Dict("T" => "barker-etal_2011.nc", "co2" => "luthi-etal_2008.nc", "V" => "spratt-lisiecki_2016.nc"))
    # 1. Define some local variables and check if ensemble
    out_path = pwd() * "/output/" * experiment * "/"
    proxy_path = pwd() .* "/data/"

    if ~isdir(proxy_path)
        (plot_proxies) && (printstyled("can't find proxy data directory \n", color=:red))
        plot_proxies = false
    end

    if experiments == []
        elements, isensemble = is_ensemble(experiment)
        if isensemble
            data_to_load = out_path .* elements .* "/pacco.nc"
        else
            data_to_load = out_path .* ["/pacco.nc"]
        end
    else
        isensemble = false
        data_to_load = pwd() .* "/output/" .* experiments .* "/pacco.nc"
    end

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
        colors_2 = cgrad(:copper, length(data_to_load), categorical = true, rev = true)
    end

    # 2. Load data
    time, data, units = [], [], []          # data[variable][element][values]
    for v in vars2plot
        datav = []
        for e in data_to_load
            if e == data_to_load[1]
                if v == vars2plot[1]
                    time = NCDataset(e, "r") do ds
                        ds["time"][:]
                    end # ds is closed
                end

                units_v = NCDataset(e, "r") do ds
                    ds[v].attrib["units"]
                end # ds is closed
                push!(units, units_v)
            end

            data_ve = NCDataset(e, "r") do ds
                ds[v][:]
            end # ds is closed
            push!(datav, data_ve)
        end
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
                ax_PSD = Axis(fig[v, 2], title="Normalized PSD", titlesize=0.8 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period (kyr)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz)
            else
                ax_PSD = Axis(fig[v, 2], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period (kyr)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6 * fntsz, yticklabelsize=0.7 * fntsz)
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
        (time[end] > time_anth) && (vlines!(ax, [time_anth], linewidth=3, color=(:red, 0.5)))


        # -- 3.6 If desired, plot proxy files
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
                    lines!(ax_PSD, periods, G, color=:black, linewidth=4)
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
                colors_experiments_mode = [:royalblue4, :red4, :olive, :purple]
                color_to_use = colors_experiments_mode[e]
            end

            if var_v == "T" # CHECK THIS, CUTRE
                l = lines!(ax, time, data_v[e][:] .- 273.15, linewidth=linewidth, color=color_to_use)
            else
                l = lines!(ax, time, data_v[e][:], linewidth=linewidth, color=color_to_use)
            end
            push!(lins, l)

            if plot_PSD
                G, f = calc_spectrum(data_v[e][:], 1 / (time[2] - time[1]))
                G, f = G[f.>1/150e3], f[f.>1/150e3]   # filtering
                G = G ./ sum(G)
                periods = 1 ./ f
                lines!(ax_PSD, periods, G, color=color_to_use, linewidth=linewidth)
            end
        end

        # -- 3.8 Formatting ...
        xlims!(ax, (time[1], time[end]))

        if (time[end] - time[1]) > 2.0e6
            tstep = 300e3
        elseif (time[end] - time[1]) > 8e5
            tstep = 200e3
        else
            tstep = 100e3
        end
        ax.xticks = -5e6:tstep:5e6
        ax.xtickformat = k -> string.(Int.(k / 1000))

        if plot_PSD
            ax_PSD.xticks = 0:20e3:150e3
            ax_PSD.xtickformat = k -> string.(Int.(ceil.(k / 1000)))
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
    t_coi, p_coi = calc_coi(time, freqs, 2*pi / sigma)  # convert rad/s to Hz (for cf)
    p_coi = log10.(p_coi ./ 1e3)

    ## Plot
    fig, fntsz = Figure(resolution=(800, 400)), 20

    ax = Axis(fig[1, 1], titlesize=0.8 * fntsz,
        xlabelsize=fntsz, ylabelsize=fntsz, xlabel="Time (kyr)", ylabel="Period (kyr)")
    ax2 = Axis(fig[1, 1], ylabelsize=fntsz, yaxisposition = :right)

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

    t1, p1 = t_coi[1:Int(length(p_coi)/2)], p_coi[1:Int(length(p_coi)/2)]
    t2, p2 = t_coi[Int(length(p_coi)/2)+1:end], p_coi[Int(length(p_coi)/2)+1:end]

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

    ax2.ytickformat = k -> string.(Int.(ceil.(k ./ 100) .* 100)) .* " $(df[var2plot].attrib["units"])"
    ax2.yticks = range(minimum(data), stop=maximum(data), length=3)
    ylims!(ax, (1, maximum(p_coi)))
    ylims!(ax2, (minimum(data), 8*maximum(data)))

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