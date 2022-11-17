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
    

    ## Plotting
    # -- Plot variables
    for i in 1:nrows
        if i == 1
            ax = Axis(fig[i, 1], title = "AMOD variables", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time", ylabel=vrs[i]*" ("*d[i].attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], title="Normalized Power density", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" ("*x.attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        else
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time", ylabel=vrs[i]*" ("*d[i].attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" ("*x.attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        end
        if i < nrows
            hidexdecorations!(ax, grid=false)
            hidexdecorations!(ax_spect, grid=false)
        end
        
        update_theme!()
        
        lines!(ax, x, d[i], linewidth=3, color=clrs[i])
        barplot!(ax_spect, f[i], G[i], linewidth=3, color=clrs[i])

        xlims!(ax, (x[1], x[end]))
        xlims!(ax_spect, (0, 1/20e3))
        ylims!(ax_spect, (0.0, 1.0))
        ax_spect.xticks = [1/100e3, 1/41e3, 1/23e3]
        ax_spect.xtickformat = k -> string.(ceil.(1 ./ k))
        ax_spect.yticks = 0:0.2:1.0
    end

    # Resizing
    colsize!(fig.layout, 2, Relative(0.33))
    resize_to_layout!(fig)

    # Saving
    save(plotpath, fig)
end