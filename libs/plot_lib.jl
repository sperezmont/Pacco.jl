# =============================
#     Program: plot_lib.jl
#     Aim: functions for plotting
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

using NCDatasets
using CairoMakie

function plot_varspectra(d::Any, vrs::Any, plotpath::String)
    ## First determine plot parameters
    # -- number of rows and columns
    nrows, ncols = length(vrs), 2

    # -- figure size
    fgsz = (1000*ncols, 750*nrows)

    # Now, define the figure and static elements
    fig = Figure(resolution=fgsz)

    save(plotpath, fig)
end