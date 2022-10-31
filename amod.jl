# =============================
#     Program: amod.jl
#     Aim: AMOD executable
#     How to run: julia amod.jl [OUT] [PARF] [NPARS]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [PARF]  --> parameters filename (without "par/", default=amod_default.jl)
#             [NPARS] --> new parameters (changes)
# =============================

# First, load external libraries, activate amod environment and save pwd 
using Pkg
if Pkg.project().name != "amod_env" # check if amod_env is activated
    Pkg.activate("amod_env")
end
amod_path = pwd()

# Second, load dependencies
include("./libs/load_args.jl")
include("./libs/nc.jl")
include("./src/amod_update.jl")

# Now we load command line arguments
output_path, outfldr = load_out(amod_path, ARGS)
parf = load_parf(amod_path, output_path, ARGS)

# Assign parameters
include("./src/amod_defs.jl") # assign parameters, CTL, INCOND, PAR, amod_INCOND

# open out.out
f = open(output_path*"/out.out", "w")
write(f, "**** Starting AMOD "*outfldr*" ... ****\n")

# Initialize
NOW = amod_INCOND # NOW contains the values of time step = time
OUT = update_amod_out(OUT, NOW)
write(f, "time = "*string(NOW["time"])*" --> "*"H = "*string(NOW["H"])*"\n") 

time_length = ceil((CTL["time_end"] - CTL["time_init"])/CTL["dt"])
for n in 1:time_length
    NOW["time"] = CTL["time_init"] + n*CTL["dt"] # update simulation time

    # Update contour variables

    # Run AMOD
    include("./src/lib_amod.jl")
    run_amod(CTL, PAR, NOW)
    
    if mod(NOW["time"], CTL["dt_out"]) == 0 # only updates output variable at desired frequency
        global OUT = update_amod_out(OUT, NOW)
        write(f, "time = "*string(NOW["time"])*" --> "*"H = "*string(NOW["H"])*"\n") 
    end
end

# Create output nc file
isfile(output_path*"/amod.nc") && rm(output_path*"/amod.nc")    # if outfile exists remove it
genout_nc(output_path, "amod.nc", OUT, out_precc, out_attr);

write(f, "**** AMOD "*outfldr*" done! ****")
close(f)


