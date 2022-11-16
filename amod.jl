# =============================
#     Program: amod.jl
#     Aim: AMOD executable
#     How to run: julia amod.jl [OUT] [PARF]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [PARF]  --> parameters filename (without "par/", default=amod_default.jl)
# =============================

## First, load external libraries, activate amod environment and save pwd 
using Pkg
Pkg.activate("amod_env")
amod_path = pwd()

## Second, load dependencies
include("./libs/load_args.jl")
include("./libs/nc.jl")
include("./src/amod_update.jl")
include("./src/lib_amod.jl")
include("./src/amod_radiative.jl")
include("./src/amod_orbital.jl")

## Now we load command line arguments
output_path, outfldr = load_out(amod_path, ARGS)
parf = load_parf(amod_path, output_path, ARGS)

# Assign parameters
# -- assign parameters, CTL, INCOND, PAR, amod_INCOND
include("./src/amod_defs.jl")

# -- load Earth constants 
include("./par/earth_const.jl")

## Open out.out
# if out.out exists remove it
isfile(output_path * "/out.out") && rm(output_path * "/out.out")
f = open(output_path * "/out.out", "w")
write(f, "**** Starting AMOD " * outfldr * " ... ****\n")

## Initialize
NOW = copy(amod_INCOND)
OUT = update_amod_out(OUT, NOW) # update output
write(f, "time = " * string(NOW["time"]) * " --> " * "ins = " * string(NOW["ins"]) * " --> " * "T_sl = " * string(NOW["T_sl"]) * " --> " * "H = " * string(NOW["H"]) * "\n")

## Let's run!
OUT = amod_loop(NOW, OUT, PAR, CTL, f)

## Create output nc file
# if outfile exists remove it
isfile(output_path * "/amod.nc") && rm(output_path * "/amod.nc")
genout_nc(output_path, "amod.nc", OUT, out_precc, out_attr);

write(f, "**** AMOD " * outfldr * " done! ****" * "\n")
close(f)

## Done!
