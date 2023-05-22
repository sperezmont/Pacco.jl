# =============================
#   Program: config.jl
#   Aim: This script configures the dependencies for using PACCO
#   Author: Sergio Pérez Montero
#   Date: 25.10.2022
# =============================
using Pkg

# Check if Manifest.toml and Project.toml exist
if isfile("Manifest.toml")
    display("** Checking Manifest.toml ... **")
    Pkg.activate(".")
    Pkg.resolve()
    Pkg.instantiate()
else
    # Environment generation
    Pkg.activate(".")

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

# Check if data/ directory is present
if isdir("data") == false
    println("`data` directory is not present")
end

