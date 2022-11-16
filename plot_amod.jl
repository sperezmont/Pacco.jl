# =============================
#     Program: plot_amod.jl
#     Aim: Plot AMOD output
#     How to run: julia amod.jl [OUT] [VARS]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [VARS]  --> variables to plot (up to 4 variables)
# =============================

## First, load external libraries, activate amod environment and save pwd 
using Pkg
Pkg.activate("amod_env")
amod_path = pwd()

include("./libs/nc.jl")
include("./libs/plot_lib.jl")
include("./libs/spectrum_lib.jl")

## Check arguments
if ARGS == []
    experiment = "test_default"
    vars = ["ins", "T_sl", "H", "Hsed"]
else
    experiment = ARGS[1]
    vars = ARGS[2:end]
end

# -- include namelist file of experiment
include("./output/"*experiment*"/namelist.jl")

## Load output data 
data, time = load_nc(amod_path*"/output/"*experiment*"/amod.nc", vars)

## Calculate spectra
G_data, freqs_data = [], []
for v in 1:length(vars)
    new_data = data[v][2:end] - data[v][1:end-1] # detrending
    
    G_v, freqs_v = calc_spectrum(new_data, 1 / dt_out)
    N = length(G_v)
    if mod(N, 2) == 0
        N0 = Int(N / 2) + 2
    else
        N0 = Int((N-1) / 2) + 2
    end
    G_v, freqs_v = G_v[N0:end], freqs_v[N0:end]

    push!(G_data, G_v ./ (sum(G_v)))  # normalization through the entire integral (only positive values)
    push!(freqs_data, freqs_v)
end

## Plot
plot_spectrum(time, data, freqs_data, G_data, vars, [:orange, :red, :blue, :green], amod_path*"/output/"*experiment*"/"*"amod_spectra.png");
