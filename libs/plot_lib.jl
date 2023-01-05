# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Functions
@doc """
    plot_spectrum: plots time series and their spectrums 
"""
function plot_spectrum(x, d::Any, f::Any, G::Any, vrs::Any, plotpath::String; fntsz=nothing)
    ## First determine plot parameters
    # -- number of rows and columns
    nrows, ncols = length(vrs), 2

    # -- figure size
    fgsz = (1500 * ncols, 500 * nrows)

    # -- check fontsize
    isnothing(fntsz) && (fntsz = 0.01 * sqrt(fgsz[1]^2 + fgsz[2]^2))
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Now, define the figure and static elements
    fig = Figure(resolution=fgsz)

    ## Plotting
    # -- Plot variables
    for i in 1:nrows
        if length(vrs) == 1
            di = d
        else
            di = d[i]
        end

        if i == 1
            ax = Axis(fig[i, 1], title="AMOD variables", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vrs[i] * " ( " * di.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], title="Normalized Power density", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (" * x.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        else
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vrs[i] * " (" * di.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period" * " (k" * x.attrib["units"] * ")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        end
        if i < nrows
            hidexdecorations!(ax, grid=false)
            hidexdecorations!(ax_spect, grid=false)
        end

        update_theme!()

        #maxdi, mindi = maximum(di[:]), minimum(di[:])
        #maxi = max(abs(maxdi), abs(mindi))
        if var(di) < 1e-1
            lines!(ax, x, di, linewidth=3)
        else
            lines!(ax, x, di, linewidth=3)
        end

        lines!(ax_spect, f[i], G[i], color=:black, linewidth=3)
        band!(ax_spect, f[i], 0.0, G[i], linewidth=2)

        xlims!(ax, (x[1], x[end]))
        xlims!(ax_spect, (1 / 500e3, 1 / 21e3))

        xlen = length(x)
        if mod(xlen, 2) == 0
            xstep = Int(xlen / 10)
        else
            xstep = Int((xlen - 1) / 10)
        end
        ax.xticks = x[1:xstep:end]
        ax.xtickformat = k -> string.(k / 1000)

        ax_spect.xticks = [1 / 100e3, 1 / 41e3, 1 / 23e3]
        ax_spect.xtickformat = k -> string.(Int.(ceil.(1 ./ k) / 1000))
    end

    # Resizing
    colsize!(fig.layout, 2, Relative(0.33))
    resize_to_layout!(fig)

    # Saving
    save(plotpath, fig)
end


@doc """
    plot_amod: calculates spectrum and plots results from AMOD (given or not the variables to plot)
"""
function plot_amod(; experiment="test_default", vars2plot=["ins_n", "SMB_n", "H_n", "Hsed_n"])
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", vars2plot)

    ## Calculate spectra
    G_data, freqs_data = [], []
    for v in 1:length(vars2plot)
        if length(vars2plot) == 1
            new_data = copy(data)
        else
            new_data = copy(data[v])
        end

        # -- spectrum blackman tuckey
        N = length(new_data)
        Nmax = Int(ceil(N / 2))
        P = periodogram(new_data, onesided=false, fs=1 / (time[2] - time[1]), window=blackman(N))
        G, freq = P.power, P.freq
        G, freq = G[1:Nmax] .* 2, freq[1:Nmax]
        G, freq = G[freq.>=1/150.5e3], freq[freq.>=1/150e3] # we eliminate values above and below Milankovitch cycles
        G, freq = G[freq.<=1/22e3], freq[freq.<=1/22e3]
        G = sqrt.(G) / (N / 2)
        G = G ./ sum(G)
        if var(new_data) < 1e-1
            G = zeros(length(G))
        end
        push!(G_data, G)
        push!(freqs_data, freq)  # save frequencies
    end

    ## Plot
    plot_spectrum(time, data, freqs_data, G_data, vars2plot, amod_path * "/output/" * experiment * "/" * "amod_results.png")
end

@doc """
    plot_wavelet: Plots a map of the wavelet scalogram
"""
function plot_wavelet(; experiment="test_default", var2plot="H_n", fs=1 / 1000)
    ## Load output data 
    data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", [var2plot])

    ## Check if not dyadic size, else, interpolate to dyadic size 
    if ~isdyadic(data[:])
        data_dyadic = vector2dyadic(copy(data[:]))
    else
        data_dyadic = copy(data[:])
    end

    ## Compute Wavelets
    Wnorm, freqs = calc_wavelet(data_dyadic, fs)

    ## Plot
    fig, fntsz = Figure(resolution=(800, 600)), 0.01 * sqrt(800^2 + 600^2)
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    ax = Axis(fig[1, 1], title=var2plot * " (" * data.attrib["units"] * ")",
        xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel="Period (kyr)")
    update_theme!()

    cmap = :vik
    minW, maxW = minimum(Wnorm), maximum(Wnorm)
    stepW = 0.1 * max(minW, maxW)
    c = contourf!(ax, vector2dyadic(time), freqs, Wnorm, colormap=cmap, levels=minW:stepW:maxW)
    c.extendlow = :auto
    c.extendhigh = :auto
    Colorbar(fig[1, 2], c, height=Relative(1 / 3), width=20, label="Normalized power density", ticklabelsize=fntsz)

    xlen = length(time)
    if mod(xlen, 2) == 0
        xstep = Int(xlen / 10)
    else
        xstep = Int((xlen - 1) / 10)
    end
    ax.xticks = time[1:xstep:end]
    ax.xtickformat = k -> string.(k / 1000)
    xlims!(ax, (time[1], time[end]))

    ax.yticks = [1 / 100e3, 1 / 41e3, 1 / 23e3]
    ax.ytickformat = k -> string.(Int.(ceil.(1 ./ k) / 1000))
    ylims!(ax, (1 / 500e3, 1 / 21e3))

    save(amod_path * "/output/" * experiment * "/" * var2plot * "_wavelet.png", fig)
end