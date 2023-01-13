# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Functions
@doc """
    plot_amod: calculates spectrum and plots results from AMOD (given or not the variables to plot)
"""
function plot_amod(; experiment="test_default", vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"])
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", vars2plot)

    if ("T_n" in vars2plot)
        data_snyder2016, time_snyder2016 = load_nc(amod_path * "/data/Snyder_2016/snyder_2016.nc", ["GAST(2.5%)", "GAST(5%)", "GAST(25%)", "GAST(50%)",	"GAST(75%)", "GAST(95%)", "GAST(97.5%)"]; time_name="Time")
    end

    if ("co2_n" in vars2plot)
        data_luthi2008, time_luthi2008 = load_nc(amod_path * "/data/Luthi-etal_2008/luthi-etal_2008.nc", ["CO2"]; time_name="time")
    end

    if ("V_n" in vars2plot)
        # Waelbroeck 2002
        data_waelbroeck2002, time_waelbroeck2002 = load_nc(amod_path * "/data/Waelbroeck-etal_2002/waelbroeck-etal_2002.nc", ["RSL-", "RSL", "RSL+"]; time_name="Age")
        # Spratt and Lisiecki 2016
        if time[1] > -5e5
            spratt_vars = ["SeaLev_shortPC1_err_lo", "SeaLev_shortPC1", "SeaLev_shortPC1_err_up"]
        else
            spratt_vars = ["SeaLev_longPC1_err_lo", "SeaLev_longPC1", "SeaLev_longPC1_err_up"]
        end
        data_spratt2016, time_spratt2016 = load_nc(amod_path * "/data/Spratt-Lisiecki_2016/spratt-lisiecki_2016.nc", spratt_vars; time_name="age_calkaBP")
    end

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
    fgsz = (1600 * ncols, 600 * nrows)

    # -- check fontsize
    fntsz = 0.02 * sqrt(fgsz[1]^2 + fgsz[2]^2)
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Now, define the figure and static elements
    fig = Figure(resolution=fgsz)

    ## Plotting
    # -- Plot variables
    color_list = collect(cgrad(:darkrainbow, nrows, categorical=true, rev=true))
    palettes=(color=color_list,)
    for i in 1:nrows
        if length(vars2plot) == 1
            di = data
        else
            di = data[i]
        end

        if i == 1
            ax = Axis(fig[i, 1], title="AMOD variables", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[i] * " ( " * di.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.5*fntsz, yticklabelsize=0.5*fntsz)
            ax_spect = Axis(fig[i, 2], title="Normalized PSD", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (" * time.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.5*fntsz, yticklabelsize=0.5*fntsz)
        else
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vars2plot[i] * " (" * di.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.5*fntsz, yticklabelsize=0.5*fntsz)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (k" * time.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey, xticklabelsize=0.5*fntsz, yticklabelsize=0.5*fntsz)
        end
        if i < nrows
            hidexdecorations!(ax, grid=false)
            hidexdecorations!(ax_spect, grid=false)
        end

        update_theme!()

        if vars2plot[i] == "V_n"
            band!(ax, time_waelbroeck2002 .* 1000, data_waelbroeck2002[1], data_waelbroeck2002[3], linewidth=2, color=(:black, 0.2))
            lines!(ax, time_waelbroeck2002 .* 1000, data_waelbroeck2002[2], linewidth=3, color=:black, label="Waelbroeck et al. (2002)")
            band!(ax, time_spratt2016 .* 1000, data_spratt2016[1], data_spratt2016[3], linewidth=3, color=(:red, 0.2))
            lines!(ax, time_spratt2016 .* 1000, data_spratt2016[2], linewidth=3, color=:red, label="Spratt and Lisiecki (2016)")
            lines!(ax, time, di, linewidth=5, color=palettes[1][i], label="AMOD")
            axislegend(ax, position=:lb, labelsize=0.5*fntsz)
        elseif vars2plot[i] == "co2_n"
            lines!(ax, time_luthi2008 .* 1000, data_luthi2008, linewidth=3, color=:black, label="Lüthi et al. (2008)")
            lines!(ax, time, di, linewidth=5, color=palettes[1][i], label="AMOD")
            axislegend(ax, position=:lb, labelsize=0.5*fntsz)
        elseif vars2plot[i] == "T_n"
            snyder = Array([i[:] for i in data_snyder2016])
            snyder_up = maximum(snyder, dims=1)[1]
            snyder_lo = minimum(snyder, dims=1)[1]
            lines!(ax, time_snyder2016 .* 1000, data_snyder2016[4], linewidth=3, color=:black, label="Snyder (2016)")
            band!(ax, time_snyder2016 .* 1000, snyder_lo, snyder_up, color=(:black, 0.2))
            lines!(ax, time, di .- 273.15, linewidth=5, color=palettes[1][i], label="AMOD")  # cutre, tengo que ver cómo hacer lo de la referencia
            axislegend(ax, position=:lb, labelsize=0.5*fntsz)
        else
            lines!(ax, time, di, linewidth=5, color=palettes[1][i])
        end        

        vlines!(ax_spect, [21e3, 41e3, 100e3], linewidth=3, color=:black, linestyle=:dash)
        barplot!(ax_spect, periods_data[i], G_data[i], width=fgsz[1], color=palettes[1][i])

        xlims!(ax, (time[1], time[end]))
        xlims!(ax_spect, (0, 120e3))

        ylims!(ax_spect, (0.0, 0.5))

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
    colsize!(fig.layout, 2, Relative(0.33))
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