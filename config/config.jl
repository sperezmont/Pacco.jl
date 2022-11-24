
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
Pkg.add("NCDatasets")
Pkg.add("DataStructures")
Pkg.add("Insolation")
Pkg.add("CairoMakie")
Pkg.add("DSP")
Pkg.add("FFTW")
Pkg.add("Wavelets")

# Check status
Pkg.precompile()
Pkg.instantiate()
display("**** AMOD ready ****")

