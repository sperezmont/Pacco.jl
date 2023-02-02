
"""
Author: Sergio Pérez Montero\n
Date: 25.10.2022\n

Aim: This script configures the dependencies for using AMOD, by Jorge Alvarez-Solas (julia version)\n

"""

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
packages = ["NCDatasets", "DataStructures", "Insolation",
            "CairoMakie", "DSP", "Wavelets", "ContinuousWavelets", "Statistics", "StatsBase", "Interpolations", "LatinHypercubeSampling"]
for i in packages
    Pkg.add(i)
end

# Precompile -- last very long... move to REPL
# using PackageCompiler
# create_sysimage(packages, sysimage_path="sys_amod.so", precompile_execution_file=["amod.jl", "plot_amod.jl"])

# Check status
Pkg.precompile()
Pkg.instantiate()
display("**** AMOD ready ****")

