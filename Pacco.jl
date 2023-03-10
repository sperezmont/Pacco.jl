# =============================
#     Program: pacco.jl
#     Aim: prepare Julia for simulations of PACCO
# =============================
println("Getting Julia ready to run PACCO ...")
using Pkg
Pkg.activate("env")                         # -- activate pacco virtual environment
global pacco_path = pwd()                    # -- determine pacco path. This line assumes we are working from the main PACCO directory

# -- import dependencies
# ---- pacco
using NCDatasets        # to make outputs and manage inputs (?)
using DataStructures    # to create OrderedDict's 
using Insolation        # to load orbital parameters
using ProgressBars

# ---- plotting
using CairoMakie            # plotting interface
using DSP                   # spectral analysis 
using Wavelets              # wavelet analysis
using ContinuousWavelets    # ""
using Statistics            # to make some minor calculations
using StatsBase             # ""
using Interpolations        # ""
using LatinHypercubeSampling    # ""

# -- import model libraries and functions
include("./libs/misc.jl")              # -- include libraries
include("./libs/nc.jl")
include("./libs/plot_lib.jl")
include("./libs/analysis.jl")

include("./par/earth_const.jl")             # -- include earth constants

include("./src/lib_pacco.jl")                # -- include pacco functions
include("./src/pacco_defs.jl")
include("./src/pacco_update.jl")
include("./src/forcing/radiative.jl")
include("./src/forcing/orbital.jl")
include("./src/physics/diagnostic.jl")
include("./src/physics/dynamics.jl")
include("./src/physics/thermodynamics.jl")
println("Done!")

# PACCO header
printstyled("================================================================== \n", color=:light_blue)
printstyled("| PACCO v0.4                                                     | \n", color=:bold)
printstyled("|----------------------------------------------------------------| \n")
printstyled("|    To run model:                                               | \n")
printstyled("|      --> run_pacco(;experiment, par_file, par2change)          | \n")
printstyled("|      --> run_pacco_lhs(par2per, nsim; experiment, par_file)    | \n")
printstyled("|      --> run_ensemble(par2per; experiment, par_file)           | \n")
printstyled("|                                                                | \n")
printstyled("|    To plot results:                                            | \n")
printstyled("|      --> plot_pacco(;experiment/experiments, vars2plot)        | \n")
printstyled("|      --> plot_wavelet(;experiment, var2plot)                   | \n")
printstyled("================================================================== \n", color=:light_blue)

