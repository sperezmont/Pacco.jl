# =============================
#   Program: config.jl
#   Aim: This script configures the dependencies for using PACCO
#   Author: Sergio Pérez Montero
#   Date: 25.10.2022
# =============================

# Check if environment exists
if isdir("pacco_env") # check if pacco_env is activated 
    rm("pacco_env", recursive=true)
end

# Environment generation
using Pkg
Pkg.generate("pacco_env")
Pkg.activate("pacco_env")

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
display("**** PACCO ready ****")

