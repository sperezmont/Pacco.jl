# =============================
#     Program: run_model.jl
#     Aim: prepare Julia for simulations of AMOD
# =============================
println("#### Preparing Julia to run AMOD ... ####")
using Pkg
Pkg.activate("amod_env")        # -- activate amod virtual environment
amod_path = pwd()               # -- get amod directory
include("./src/AMOD.jl")        # -- include AMOD module

println("Done!")
println("Now use AMOD.run_amod(output, input) to execute the model as many times you want")
println("#### -------- ####")

