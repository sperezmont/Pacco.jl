
"""
Author: Sergio Pérez Montero\n
Date: 25.10.2022\n

Aim: This script configures the dependencies for using AMOD, by Jorge Alvarez-Solas (julia version)\n

"""

# Environment generation
using Pkg
Pkg.generate("amod_env")
Pkg.activate("amod_env")

# Adding dependencies ... 
display("** Adding dependencies ... **")
Pkg.add("NCDatasets")
Pkg.add("DataStructures")

# Check status
Pkg.status()
display("**** AMOD ready ****")

