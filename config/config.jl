# =============================
#   Program: config.jl
#   Aim: This script configures the dependencies for using PACCO
#   Author: Sergio Pérez Montero
#   Date: 25.10.2022
# =============================
using Pkg

env_name = "env"

# Check if environment exists
if isdir(env_name)
    Pkg.activate(env_name)
    Pkg.instantiate()
else
    # Environment generation
    Pkg.generate(env_name)
    Pkg.activate(env_name)

    # Adding dependencies ... 
    display("** Adding dependencies ... **")
    packages = ["OrdinaryDiffEq", "NaNMath", "JLD2",
                "NCDatasets", "Insolation",
                "CairoMakie", "DSP", "Wavelets", "ContinuousWavelets", "Statistics",
                "StatsBase", "Interpolations", "LatinHypercubeSampling",
                "Images"]
    for i in packages
        Pkg.add(i)
    end
end

# Check status and precompile
Pkg.precompile()
display("**** PACCO ready ****")

# Check if output/ directory is present
if isdir("output") == false
    mkdir("output")
end

