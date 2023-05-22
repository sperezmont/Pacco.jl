# =============================
#     Program: Pacco.jl
#     Aim: prepare Julia for simulations of PACCO (default_params.jl)
# =============================
println("Getting Julia ready to run PACCO ...")
using Pkg
Pkg.activate(".")                         # -- activate pacco virtual environment
global pacco_path = pwd()                    # -- determine pacco path. This line assumes we are working from the main PACCO directory

# -- import dependencies
# ---- pacco
using NCDatasets        # to make outputs and manage inputs (?)
using Insolation        # to load orbital parameters
using OrdinaryDiffEq        # solver for ODEs
using JLD2                  # to save Julia objects
using NaNMath

# ---- plotting
using CairoMakie            # plotting interface
using DSP                   # spectral analysis 
using Wavelets              # wavelet analysis
using ContinuousWavelets    # ""
using Statistics            # to make some minor calculations
using StatsBase             # ""
using Interpolations        # ""
using LatinHypercubeSampling    # ""
using Images: findlocalmaxima   # we only need findlocalmaxima() (IMPORTANT!!)

# -- import model libraries and functions
include("./par/default_params.jl")           # -- includes source parameter file 

include("./libs/misc.jl")              # -- includes libraries
include("./libs/nc.jl")
include("./libs/plot_lib.jl")
include("./libs/analysis.jl")

include("./src/lib_pacco.jl")                # -- includes pacco functions
include("./src/pacco_defs.jl")
include("./src/forcing/orbital.jl")          # -- forcing functions
include("./src/physics/climate.jl")          # -- climate functions         
include("./src/physics/geometry.jl")         # -- cryosphere functions
include("./src/physics/dynamics.jl")         # 
include("./src/physics/thermodynamics.jl")   #

include("PaccoTests.jl")     # -- include tests.jl
println("Done!")

# PACCO header
printstyled("================================================================== \n", color=:light_blue)
printstyled("| PACCO v0.5                                                     | \n", color=:bold)
printstyled("|----------------------------------------------------------------| \n")
printstyled("|    To run model:                                               | \n")
printstyled("|      --> run_pacco(experiment; p)                              | \n")
printstyled("|      --> run_ensemble(experiment, params2per)                  | \n")
printstyled("|      --> run_pacco_lhs(experiment, params2per, nsim)           | \n")
printstyled("|                                                                | \n")
printstyled("|    To plot results:                                            | \n")
printstyled("|      --> plot_pacco(experiment; vars2plot)                     | \n")
printstyled("|      --> plot_pacco_states(experiment)                         | \n")
printstyled("|      --> plot_pacco_comp_states(experiment)                    | \n")
printstyled("|      --> plot_wavelet(experiment; vars2plot)                   | \n")
printstyled("|                                                                | \n")
printstyled("|    Shortcuts:                                                  | \n")
printstyled("|      --> runplot_pacco(experiment; p)                          | \n")
printstyled("|      --> runplot_pacco_ensemble(experiment, params2per)        | \n")
printstyled("|      --> runplot_pacco_lhs(experiment, params2per, nsim)       | \n")
printstyled("================================================================== \n", color=:light_blue)

