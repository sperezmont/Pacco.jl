# =============================
#     Program: amod.jl
#     Aim: prepare Julia for simulations of AMOD
# =============================
println("#### Preparing Julia to run AMOD ... ####")
using Pkg
Pkg.activate("amod_env")                    # -- activate amod virtual environment
global amod_path = pwd()                    # -- determine amod path. This line assumes we are working from the main AMOD directory

# -- import dependencies
using NCDatasets, DataStructures             # amod
using Statistics, DSP, CairoMakie, Wavelets, ContinuousWavelets, Interpolations  # plotting

include("./libs/misc.jl")              # -- include libraries
include("./libs/nc.jl")
include("./libs/plot_lib.jl")

include("./par/earth_const.jl")             # -- include earth constants

include("./src/lib_amod.jl")                # -- include amod functions
include("./src/amod_defs.jl")
include("./src/amod_update.jl")
include("./src/amod_radiative.jl")
include("./src/amod_orbital.jl")
include("./src/amod_dynamics.jl")
include("./src/amod_thermodynamics.jl")

println("Done! Now use:")
println("   --> run_amod(out_name, par_file, par2change) to execute the model")
println("   --> plot_amod(experiment, variables) to plot variables of an experiment")
println("#### -------- ####")




