# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Functions
@doc """
    plot_amod: calculates spectrum and plots results from AMOD (given or not the variables to plot)
"""
function plot_amod(; experiment="test_default", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], MPT=false)
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", vars2plot)
    include(amod_path * "/output/" * experiment * "/namelist.jl")
    proxies_data = Dict("T_proxy" => T_proxy, "co2_proxy" => co2_proxy, "V_proxy" => V_proxy)
    
    ## Calculate spectra
    G_data, freqs_data, periods_data = [], [], []
    for v in 1:length(vars2plot)
        if length(vars2plot) == 1
            new_data = copy(data)
        else
            new_data = copy(data[v])
        end
        Gv, freqv = calc_spectrum(new_data, 1/(time[2]-time[1]))
        Gv, freqv = Gv[freqv .> 1/120e3], freqv[freqv .> 1/120e3]   # filtering
        Gv = Gv ./ sum(Gv)
        periodsv = 1 ./ freqv    # kyr
        push!(G_data, Gv)
        push!(freqs_data, freqv)  # save frequencies
        push!(periods_data, periodsv)
    end

    ## Plot
    ## First determine plot parameters
    # -- number of rows and columns
    nrows, ncols = length(vars2plot), 2

    # -- figure size
    fgsz = (2000 * ncols, 600 * nrows)

    # -- check fontsize
    fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Now, define the figure and static elements
    fig = Figure(resolution=fgsz)

    ## Plotting
    # -- Plot variables
    color_list = collect(cgrad([:black, :royalblue4, :red4, :olive, :steelblue4], nrows, categorical=true))
    palettes=(color=color_list,)
    k = 1
    for i in 1:nrows
        if length(vars2plot) == 1
            di = data
        else
            di = data[i]
        end

        if i == 1
            ax = Axis(fig[i, 1], title="AMOD variables", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[i] * " ( " * di.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
            ax_spect = Axis(fig[i, 2], title="Normalized PSD", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (" * time.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
        elseif vars2plot[i] == "T_n"
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel="ΔT_n (K)", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (k" * time.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
        else
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[i] * " (" * di.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (k" * time.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.6*fntsz, yticklabelsize=0.7*fntsz)
        end
        if i < nrows
            hidexdecorations!(ax, grid=false)
            hidexdecorations!(ax_spect, grid=false)
        end

        update_theme!()

        proxy_colors = [:limegreen, :indigo, :aquamarine4, :orangered, :hotpink4]
        if vars2plot[i] in ["V_n", "co2_n", "T_n"]
            proxy_i = proxies_data[vars2plot[i][1:end-2]*"_proxy"]
            for j in eachindex(proxy_i)
                df = NCDataset(amod_path * "/data/" * proxy_i[j])
                proxy_data, proxy_time = df[vars2plot[i][1:end-2]], df["time"]
                if haskey(df, vars2plot[i][1:end-2]*"_lo")
                    proxy_err_lo = df[vars2plot[i][1:end-2]*"_lo"]
                    proxy_err_up = df[vars2plot[i][1:end-2]*"_up"]
                    band!(ax, proxy_time, proxy_err_lo, proxy_err_up, linewidth=2, color=(proxy_colors[k], 0.6/j))
                end
                idx = 7
                for n in eachindex(proxy_i[j])
                    if proxy_i[j][n] == '/'
                        idx = n
                    end
                end
                label_ij = proxy_i[j][1:idx-1]
                scatter!(ax, proxy_time, proxy_data, color=proxy_colors[k], width=5, label=label_ij)
                k += 1
            end
            if vars2plot[i] in ["T_n"]
                lines!(ax, time, di .- (t_ref_n + 273.15), linewidth=5, color=palettes[1][i], label="AMOD")
            else
                lines!(ax, time, di, linewidth=5, color=palettes[1][i], label="AMOD")
            end
            axislegend(ax, position=:lb, labelsize=0.5*fntsz)
        else
            lines!(ax, time, di, linewidth=5, color=palettes[1][i])
        end        

        
        (MPT) && (vlines!(ax, [-1.25e6, -0.7e6], linewidth=3, color=:red, linestyle=:dash))    # plot MPT
        vlines!(ax_spect, [21e3, 41e3, 100e3], linewidth=3, color=:red, linestyle=:dash)
        barplot!(ax_spect, periods_data[i], G_data[i], width=fgsz[1], color=palettes[1][i])
        lines!(ax_spect, periods_data[i], G_data[i], linestyle=:dash, linewidth=5, color=palettes[1][i])

        xlims!(ax, (time[1], time[end]))
        xlims!(ax_spect, (0, 120e3))

        ylims!(ax_spect, (0.0, 0.55))

        xlen = length(time)
        if mod(xlen, 2) == 0
            xstep = Int(xlen / 10)
        else
            xstep = Int((xlen - 1) / 10)
        end
        ax.xticks = time[1:xstep:end]
        ax.xtickformat = k -> string.(k / 1000)

        ax_spect.xtickformat = k -> string.(Int.(ceil.(k / 1000)))
    end

    # Resizing
    colsize!(fig.layout, 2, Relative(1/3))
    resize_to_layout!(fig)

    # Saving
    save(amod_path * "/output/" * experiment * "/" * "amod_results.png", fig)
end

@doc """
    plot_wavelet: Plots a map of the wavelet scalogram
"""
function plot_wavelet(; experiment="test_default", var2plot="H_n", fs=1 / 1000, sigma=π)
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", [var2plot])

    ## Compute Wavelets
    Wnorm, freqs = calc_wavelet(data, fs; sigma=sigma)
    periods = 1 ./ freqs

    ## Plot
    fig, fntsz = Figure(resolution=(1000, 600)), 20
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    ax = Axis(fig[1, 1], title=var2plot * " (" * data.attrib["units"] * ")",
        xlabelsize=fntsz, ylabelsize=fntsz, xlabel="Time (kyr)", ylabel="Period (kyr)")
    update_theme!()

    cmap = :cividis
    minW, maxW = minimum(Wnorm), maximum(Wnorm)
    stepW = 0.1 * max(minW, maxW)
    c = contourf!(ax, time, periods, Wnorm, colormap=cmap, levels=minW:stepW:maxW)
    hlines!(ax, [21e3, 41e3, 100e3], color=:red, linestyle=:dash)
    c.extendlow = :auto
    c.extendhigh = :auto
    Colorbar(fig[1, 2], c, height=Relative(1 / 3), width=20, label="Normalized PSD", ticklabelsize=fntsz)

    xlen = length(time)
    if mod(xlen, 2) == 0
        xstep = Int(xlen / 10)
    else
        xstep = Int((xlen - 1) / 10)
    end
    ax.xticks = time[1:xstep:end]
    ax.xtickformat = k -> string.(k / 1000)
    xlims!(ax, (time[1], time[end]))

    ax.ytickformat = k -> string.(Int.(ceil.(k / 1000)))
    ylims!(ax, (0, 120e3))

    save(amod_path * "/output/" * experiment * "/" * var2plot * "_wavelet.png", fig)
end