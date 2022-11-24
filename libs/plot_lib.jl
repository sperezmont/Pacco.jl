# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

using NCDatasets
using CairoMakie
using Statistics

function plot_spectrum(x, d::Any, f::Any, G::Any, vrs::Any, clrmp, plotpath::String; fntsz=nothing, fancy=false)
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
            ax = Axis(fig[i, 1], title = "AMOD variables", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vrs[i]*" ( "*d[i].attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], title="Normalized Power density", xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" ("*x.attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        else
            ax = Axis(fig[i, 1], titlesize=0.7 * fntsz, xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Time (kyr)", ylabel=vrs[i]*" ("*d[i].attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
            ax_spect = Axis(fig[i, 2], xlabelsize=0.8 * fntsz, ylabelsize=0.8 * fntsz, xlabel="Period"*" (k"*x.attrib["units"]*")", xgridcolor=:darkgrey, ygridcolor=:darkgrey)
        end
        if i < nrows
            hidexdecorations!(ax, grid=false)
            hidexdecorations!(ax_spect, grid=false)
        end
        
        update_theme!()

        maxdi, mindi = maximum(d[i][:]), minimum(d[i][:])
        maxi = max(abs(maxdi), abs(mindi))
        if maxdi * mindi < 0
            if var(d[i]) < 1e-1
                lines!(ax, x, d[i], linewidth=3, color = x, colormap = clrmp[i])
            else
                lines!(ax, x, d[i], linewidth=3, color = d[i], colormap = clrmp[i], colorrange = (-maxi, maxi))
            end   
        else
            if var(d[i]) < 1e-1
                lines!(ax, x, d[i], linewidth=3, color = x, colormap = clrmp[i])
            else
                lines!(ax, x, d[i], linewidth=3, color = d[i], colormap = clrmp[i], colorrange = (mindi, maxi))
            end
        end
        lines!(ax_spect, f[i], G[i], color=:black, linewidth=3)
        if fancy
            band!(ax_spect, f[i], 0.0, G[i], linewidth=2, color = G[i], colormap = :berlin)
        else
            band!(ax_spect, f[i], 0.0, G[i], linewidth=2, color = :darkred)
        end

        xlims!(ax, (x[1], x[end]))
        xlims!(ax_spect, (1/500e3, 1/21e3))

        xlen = length(x)
        if mod(xlen, 2) == 0
            xstep = Int(xlen / 10)
        else
            xstep = Int((xlen - 1) / 10)
        end 
        ax.xticks = x[1:xstep:end]
        ax.xtickformat = k -> string.(k/1000)

        ax_spect.xticks = [1/100e3, 1/41e3, 1/23e3]
        ax_spect.xtickformat = k -> string.(Int.(ceil.(1 ./ k)/1000)) 
    end

    # Resizing
    colsize!(fig.layout, 2, Relative(0.33))
    resize_to_layout!(fig)

    # Saving
    save(plotpath, fig)
end