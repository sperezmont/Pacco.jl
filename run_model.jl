# =============================
#     Program: run_model.jl
#     Aim: prepare Julia for simulations of AMOD
# =============================
println("#### Preparing Julia to run AMOD ... ####")
using Pkg
Pkg.activate("amod_env")                    # -- activate amod virtual environment
global amod_path = pwd()                    # -- determine amod path. This line assumes we are working from the main AMOD directory

using NCDatasets, DataStructures            # -- import dependencies

include("./libs/load_args.jl")              # -- include libraries
include("./libs/nc.jl")

include("./par/earth_const.jl")             # -- include earth constants

include("./src/lib_amod.jl")                # -- include amod functions
include("./src/amod_defs.jl")
include("./src/amod_update.jl")
include("./src/amod_radiative.jl")
include("./src/amod_orbital.jl")
include("./src/amod_dynamics.jl")
include("./src/amod_thermodynamics.jl")

println("Done!")
println("Now use AMOD.run_amod(output, input) to execute the model as many times you want")
println("#### -------- ####")

