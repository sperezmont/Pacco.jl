# =============================
#     Program: plot_experiment.jl
#     Aim: Plot AMOD output
# =============================
println("#### Preparing Julia to plot AMOD experiment ... ####")
using Pkg     
Pkg.activate("amod_env")                    # -- activate amod virtual environment

global amod_path = pwd()                    # -- determine amod path. This line assumes we are working from the main AMOD directory

using NCDatasets, DataStructures, Statistics, DSP, CairoMakie            # -- import dependencies

include("./libs/nc.jl")                     # -- include libraries and functions
include("./libs/plot_lib.jl")

println("Done!")
println("Now use PLTAMOD.plot_amod(experiment, variables) to plot variables of some experiment")
println("#### -------- ####")