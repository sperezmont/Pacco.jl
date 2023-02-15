# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Functions

function is_ensemble(exp)
    elems = readdir(pwd() * "/output/" * exp * "/")
    if "amod.nc" in elems
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

@doc """
    plot_amod: calculates spectrum and plots results from AMOD (given or not the variables to plot)
        experiment      --> experiment's name
        vars2plot       --> vector of variables to plot
        plot_MPT        --> plot MPT?   (Mid-Pleistocene Transition)
        plot_PSD        --> calculate and plot PSD? (Power Spectrum Density)
        plot_proxies    --> include proxy curves?
        proxy_files     --> dictionary with the names of the proxy files to use in T/, co2/ and V/
"""
function plot_amod(; experiment="test_default", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true, plot_proxies=true, proxy_files=Dict("T" => "barker-etal_2011.nc", "co2" => "luthi-etal_2008.nc", "V" => "spratt-lisiecki_2016.nc"))
    # 1. Define some local variables and check if ensemble
    out_path = pwd() * "/output/" * experiment * "/"
    proxy_path = pwd() .* "/data/"

    color_list = collect(cgrad([:black, :royalblue4, :red4, :olive, :steelblue4], length(vars2plot), categorical=true))
    palettes=(color=color_list,)

    elements, isensemble = is_ensemble(experiment)
    if isensemble
        data_to_load = out_path .* elements .* "/amod.nc"
    else
        data_to_load = out_path .* ["/amod.nc"]
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
    (isensemble) ? (linewidth=2) : (linewidth=4)
    for v in eachindex(vars2plot)
        data_v = data[v]
        var_v = string(vars2plot[v][1])
        # -- 3.2 Create axis
        if v == 1
            ax = Axis(fig[v, 1], title="AMOD variables", titlesize= 0.8 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[v] * " ($(units[v]))", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
        else
            ax = Axis(fig[v, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[v] * " ($(units[v]))", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
        end

        # -- 3.3 Create PSD axis if desired
        if plot_PSD
            if v == 1
                ax_PSD = Axis(fig[v, 2], title="Normalized PSD", titlesize= 0.8 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period (kyr)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
            else
                ax_PSD = Axis(fig[v, 2], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period (kyr)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
            end
        end

        # -- 3.4 Remove some decorations
        if v < nrows
            hidexdecorations!(ax, grid=false)
            (plot_PSD) && (hidexdecorations!(ax_PSD, grid=false))
        end

        # -- 3.5 Plot MPT and main Milankovitch periods
        (plot_MPT) && (vspan!(ax, -1.25e6, -0.7e6, color=(:red, 0.2)))    # plot Mid-Pleistocene Transition
        (plot_PSD) && (vlines!(ax_PSD, [21e3, 41e3, 100e3], linewidth=3, color=:red, linestyle=:dash))  # main Milankovitch periods

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
                    G, f = calc_spectrum(d_nomiss, 1/(t_nomiss[2]-t_nomiss[1]))
                    G, f = G[f .> 1/150e3], f[f .> 1/150e3]   # filtering
                    G = G ./ sum(G)
                    periods = 1 ./ f
                    lines!(ax_PSD, periods, G, color=:black, linewidth=4) 
                end
            end
        end

        # -- 3.7 Plot each simulation (and PSD if desired)
        for e in eachindex(data_to_load)
            if var_v == "T"
                lines!(ax, time, data_v[e][:] .- data_v[e][end], linewidth=linewidth, color=palettes[1][v])
            else
                lines!(ax, time, data_v[e][:], linewidth=linewidth, color=palettes[1][v])
            end
            if plot_PSD
                G, f = calc_spectrum(data_v[e][:], 1/(time[2]-time[1]))
                G, f = G[f .> 1/150e3], f[f .> 1/150e3]   # filtering
                G = G ./ sum(G)
                periods = 1 ./ f
                lines!(ax_PSD, periods, G, color=color=palettes[1][v], linewidth=linewidth) 
            end
        end

        # -- 3.8 Formatting ...
        xlims!(ax, (time[1], time[end]))
        
        if (time[end] - time[1]) > 1.5e6
            tstep = 300e3
        elseif (time[end] - time[1]) > 1e6
            tstep = 200e3
        else
            tstep = 100e3
        end
        ax.xticks = -5e6:tstep:5e6
        ax.xtickformat = k -> string.(k / 1000)

        if plot_PSD
            ax_PSD.xticks = 0:20e3:150e3
            ax_PSD.xtickformat = k -> string.(Int.(ceil.(k / 1000)))
        end

        if plot_proxies
            if var_v in collect(keys(proxy_files))
                fig[v, ncols+1] = Legend(fig, ax, framevisible=false, labelsize=0.6*fntsz)
            end
        end

    end

    (plot_PSD) && (colsize!(fig.layout, 2, Relative(1/3)))
    resize_to_layout!(fig)

    if isensemble
        save(amod_path * "/output/" * experiment * "/results/" * "amod_results.png", fig)
    else
        save(amod_path * "/output/" * experiment * "/" * "amod_results.png", fig)
    end

end

@doc """
    plot_wavelet: Plots a map of the wavelet scalogram
"""
function plot_wavelet(; experiment="test_default", var2plot="H_n", fs=1 / 1000, sigma=π, MPT=false)
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", [var2plot])

    ## Compute Wavelets
    Wnorm, freqs = calc_wavelet(data, fs; sigma=sigma)
    periods = 1 ./ freqs

    ## Plot
    fig, fntsz = Figure(resolution=(800, 400)), 20

    ax = Axis(fig[1, 1], title=var2plot * " (" * data.attrib["units"] * ")", titlesize= 0.8 * fntsz,
        xlabelsize=fntsz, ylabelsize=fntsz, xlabel="Time (kyr)", ylabel="Period (kyr)")

    cmap = :cividis
    minW, maxW = minimum(Wnorm), maximum(Wnorm)
    stepW = 0.1 * max(minW, maxW)
    c = contourf!(ax, time, periods, Wnorm, colormap=cmap, levels=minW:stepW:maxW)
    (MPT) && (vlines!(ax, [-1.25e6, -0.7e6], linewidth=3, color=:red, linestyle=:dash))    # plot MPT
    hlines!(ax, [21e3, 41e3, 100e3], linewidth=3, color=:red, linestyle=:dash)
    c.extendlow = :auto
    c.extendhigh = :auto
    Colorbar(fig[1, 2], c, height=Relative(1 / 3), width=20, label="Normalized PSD", ticklabelsize=fntsz)

    text!(ax, -1.24e6, 5e3, text="MPT starts", fontsize=0.9*fntsz, color=:red)
    text!(ax, -0.69e6, 5e3, text="MPT ends", fontsize=0.9fntsz, color=:red)

    if (time[end] - time[1]) > 1.5e6
        tstep = 300e3
    elseif (time[end] - time[1]) > 1e6
        tstep = 200e3
    else
        tstep = 100e3
    end
    ax.xticks = -5e6:tstep:5e6
    ax.xtickformat = k -> string.(k / 1000)
    xlims!(ax, (time[1], time[end]))

    ax.ytickformat = k -> string.(Int.(ceil.(k / 1000)))
    ylims!(ax, (0, 120e3))

    save(amod_path * "/output/" * experiment * "/" * var2plot * "_wavelet.png", fig)
end

@doc """
    plot_ensemble:
        plots ensemble experiment and variables in vars2plot
"""
function plot_ensemble(;experiment="test_default_ens", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], sims2highlight=[])
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
            run = NCDataset(locdir * "/" * runs[r] * "/amod.nc")
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
    save(locdir * "amod_results.png", fig)
end

# Shortcuts
@doc """
    plot_all: plots H wavelet and amod main results
"""
function plot_all(;experiment="test_default", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], MPT=false, fs=1 / 1000, sigma=π)
    plot_amod(experiment=experiment, vars2plot=vars2plot, MPT=MPT)
    plot_wavelet(experiment=experiment, var2plot="H_n", fs=fs, sigma=sigma, MPT=MPT)
    plot_wavelet(experiment=experiment, var2plot="H_n", fs=fs, sigma=sigma, MPT=MPT) # cutre, arreglar -- spm
end

@doc """
"""
function runplot_amod(; experiment="test_default", par_file="amod_default.jl", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], MPT=false, fs=1/1000, sigma=π)
    run_amod(experiment=experiment, par_file=par_file)
    plot_all(experiment=experiment, vars2plot=vars2plot, MPT=MPT, fs=fs, sigma=sigma)
end