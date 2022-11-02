# =============================
#     Program: amod.jl
#     Aim: AMOD executable
#     How to run: julia amod.jl [OUT] [PARF] [NPARS]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [PARF]  --> parameters filename (without "par/", default=amod_default.jl)
#             [NPARS] --> new parameters (changes)
# =============================

## First, load external libraries, activate amod environment and save pwd 
using Pkg
if Pkg.project().name != "amod_env" # check if amod_env is activated # RETHINK -- spm
    Pkg.activate("amod_env")
end
amod_path = pwd()

## Second, load dependencies
include("./libs/load_args.jl")
include("./libs/nc.jl")
include("./src/amod_update.jl")
include("./src/lib_amod.jl")

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
global f = open(output_path * "/out.out", "w")
write(f, "**** Starting AMOD " * outfldr * " ... ****\n")

## Initialize
# NOW contains the values of time step = time
global NOW = amod_INCOND 
OUT = update_amod_out(OUT, NOW)
write(f, "time = " * string(NOW["time"]) * " --> " * "H = " * string(NOW["H"]) * "\n")

## Let's run!
time_length = ceil((CTL["time_end"] - CTL["time_init"]) / CTL["dt"])
for n in 1:time_length 

    # update simulation time
    NOW["time"] = CTL["time_init"] + n * CTL["dt"]

    # update contour variables

    # run AMOD
    run_amod(CTL, PAR, NOW)

    # only update output variable at desired frequency
    if mod(NOW["time"], CTL["dt_out"]) == 0 
        global OUT = update_amod_out(OUT, NOW)
        write(f, "time = " * string(NOW["time"]) * " --> " * "H = " * string(NOW["H"]) * "\n")
    end
end

## Create output nc file
# if outfile exists remove it
isfile(output_path * "/amod.nc") && rm(output_path * "/amod.nc")    
genout_nc(output_path, "amod.nc", OUT, out_precc, out_attr);

write(f, "**** AMOD " * outfldr * " done! ****"*"\n")
close(f)

## Done!
