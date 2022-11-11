# =============================
#     Program: plot_amod.jl
#     Aim: Plot AMOD output
#     How to run: julia amod.jl [OUT] [VARS]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [VARS]  --> variables to plot (up to 4 variables)
# =============================

## First, load external libraries, activate amod environment and save pwd 
using Pkg
if Pkg.project().name != "amod_env" # check if amod_env is activated # RETHINK -- spm
    Pkg.activate("amod_env")
end
amod_path = pwd()

include("./libs/nc.jl")
include("./libs/plot_lib.jl")

## Check arguments
if ARGS == []
    experiment = "test_default"
    vars = ["ins", "T_sl", "H", "Hsed"]
else
    experiment = ARGS[1]
    vars = ARGS[2:end]
end

## Load output data 
data = load_nc(amod_path*"/output/"*experiment*"/amod.nc", vars)

## Calculate spectra


## Plot
#plot_varspectra(data, vars, amod_path*"/output/"*experiment*"/")