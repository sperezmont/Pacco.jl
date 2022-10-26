# =============================
#     Program: amod.jl
#     Aim: AMOD executable
#     How to run: julia amod.jl [OUT] [PARF] [NPARS]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [PARF]  --> parameters filename (without "par/", default=amod_default.jl)
#             [NPARS] --> new parameters (changes)
# =============================

# First, activate amod environment and save pwd
using Pkg
Pkg.activate("amod_env")
amod_path = pwd()

# Second, load dependencies (modules)
include("./libs/load_args.jl")
include("./libs/nc.jl")

include("./src/amod_defs.jl")
include("./src/amod_initialization.jl")
#include("./src/amod_orbital")
#include("./src/amod_dynamics")
#include("./src/amod_thermodynamics")

# Now we load command line arguments
output_path = load_out(amod_path, ARGS)
load_parf(amod_path, output_path, ARGS)

ctl, inicond, radpar, orbpar, geopar, icepar = assign_parameters() # assing parameters

# Initialization
now = amod_init()

# Open output nc file
open_nc()




