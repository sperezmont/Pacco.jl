# =============================
#     Program: runner.jl
#     Aim: runs AMOD from Julia REPL
#     How to run: julia runner.jl [FILE]
#         where
#             [FILE] --> text file where each line represents an experiment
# =============================
using Pkg
Pkg.activate("amod_env")
exp_name, parf_name = ARGS[1], ARGS[2]
run(`julia amod.jl`)

