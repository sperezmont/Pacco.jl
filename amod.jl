# =============================
#     Program: amod.jl
#     Aim: prepare Julia for simulations of AMOD
# =============================
println("Getting Julia ready to run AMOD ...")
using Pkg
Pkg.activate("amod_env")                    # -- activate amod virtual environment
global amod_path = pwd()                    # -- determine amod path. This line assumes we are working from the main AMOD directory

# -- import dependencies
# ---- amod
using NCDatasets        # to make outputs and manage inputs (?)
using DataStructures    # to create OrderedDict's 
using Insolation        # to load orbital parameters

# ---- plotting
using CairoMakie            # plotting interface
using DSP                   # spectral analysis 
using Wavelets              # wavelet analysis
using ContinuousWavelets    # ""
using Statistics            # to make some minor calculations
using Interpolations        # ""

# -- import model libraries and functions
include("./libs/misc.jl")              # -- include libraries
include("./libs/nc.jl")
include("./libs/plot_lib.jl")

include("./par/earth_const.jl")             # -- include earth constants

include("./src/lib_amod.jl")                # -- include amod functions
include("./src/amod_defs.jl")
include("./src/amod_update.jl")
include("./src/amod_radiative.jl")
include("./src/amod_orbital.jl")
include("./src/amod_diag.jl")
include("./src/amod_dynamics.jl")
include("./src/amod_thermodynamics.jl")
include("./src/amod_derivative.jl")
println("Done!")

# AMOD header
printstyled("================================================================== \n", color=:light_blue)
printstyled("| AMOD v0.2                                                      | \n", color=:bold)
printstyled("|----------------------------------------------------------------| \n")
printstyled("|    To run model:                                               | \n")
printstyled("|      --> run_amod(out_name, par_file, par2change)              | \n")
printstyled("|      --> run_ensemble(par2per, out_name, par_file)             | \n")
printstyled("|                                                                | \n")
printstyled("|    To plot results:                                            | \n")
printstyled("|      --> plot_amod(experiment, variables)                      | \n")
printstyled("|      --> plot_wavelet(experiment, variable)                    | \n")
printstyled("================================================================== \n", color=:light_blue)




