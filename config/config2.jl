# =============================
#   Program: config.jl
#   Aim: This script configures the dependencies for using AMOD
#   Author: Sergio Pérez Montero
#   Date: 25.10.2022
# =============================

# Check if environment exists
if isdir("amod_env") # check if amod_env is activated 
    rm("amod_env", recursive=true)
end

# Environment generation
using Pkg
Pkg.generate("amod_env")
Pkg.activate("amod_env")

# Adding dependencies ... 
display("** Adding dependencies ... **")
packages = ["NCDatasets", "DataStructures", "Insolation", "ProgressBars"
            "CairoMakie", "DSP", "Wavelets", "ContinuousWavelets", "Statistics",
            "StatsBase", "Interpolations", "LatinHypercubeSampling"]
for i in packages
    Pkg.add(i)
end

# Check status
Pkg.precompile()
Pkg.instantiate()
display("**** AMOD ready ****")

