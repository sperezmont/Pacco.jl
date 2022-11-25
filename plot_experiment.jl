# =============================
#     Program: plot_experiment.jl
#     Aim: Plot AMOD output
# =============================
println("#### Preparing Julia to plot AMOD experiment ... ####")
using Pkg
Pkg.activate("amod_env")        # -- activate amod virtual environment
amod_path = pwd()               # -- get amod directory
include("./src/PLTAMOD.jl")   # -- include PLOT_AMOD module

println("Done!")
println("Now use PLTAMOD.plot_amod(experiment, variables) to plot variables of some experiment")
println("#### -------- ####")