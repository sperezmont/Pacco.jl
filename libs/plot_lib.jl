# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

using NCDatasets
using CairoMakie

function plot_spectrum(x, d::Any, f::Any, G::Any, vrs::Any, clrs, plotpath::String; fntsz=nothing)
    ## First determine plot parameters
    # -- number of rows and columns
    nrows, ncols = length(vrs), 2

    # -- figure size
    fgsz = (1500*ncols, 500*nrows)

    # -- check fontsize
    isnothing(fntsz) && (fntsz = 0.01 * sqrt(fgsz[1]^2 + fgsz[2]^2))
    fontsize_theme = Theme(font="Dejavu Serif", fontsize=fntsz)
    set_theme!(fontsize_theme)

    # Now, define the figure and static elements
    fig = Figure(resolution=fgsz)
    ax_spect = Axis(fig[1:nrows, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" ("*x.attrib["units"]*")", ylabel="Normalized Power")
    #ax_spect.xtickformat = k -> string.(1 ./ k)

    ## Plotting
    # -- Plot variables
    for i in 1:nrows
        ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time", ylabel=vrs[i]*" ("*d[i].attrib["units"]*")")
        if i < nrows
            hidexdecorations!(ax, ticks=false)
        end
        
        update_theme!()
        
        lines!(ax, x, d[i], linewidth=3, color=clrs[i])
        lines!(ax_spect, f[i], G[i], linewidth=3, color=clrs[i])

        xlims!(ax, (x[1], x[end]))
    end

    xlims!(ax_spect, (0, 1/20e3))
    ax_spect.xticks = [1/100e3, 1/41e3, 1/23e3]
    ax_spect.xtickformat = k -> string.(ceil.(1 ./ k))

    ax_spect.yticks = 0:0.1:1.0

    # Resizing
    resize_to_layout!(fig)

    # Saving
    save(plotpath, fig)
end