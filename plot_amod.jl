# =============================
#     Program: plot_amod.jl
#     Aim: Plot AMOD output
#     How to run: julia plot_amod.jl [OUT] [VARS]
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
    vars = ["ins_norm", "T_sl", "T_surf", "TMB", "H", "B", "Z", "U", "Hsed"]
elseif length(ARGS) == 1
    experiment = ARGS[1]
    vars = ["ins", "TMB", "H", "Hsed"]
else
    experiment = ARGS[1]
    vars = ARGS[2:end]
end

# -- include namelist file of experiment
include("./output/" * experiment * "/namelist.jl")

## Load output data 
data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", vars)

## Calculate spectra
#window = 50
G_data, freqs_data, G_array, freqs_array = [], [], [], []
#G_array = Array{Float64}(undef, Int(length(vars)), Int(length(data[1]) / window), Int(window / 2 - 1))
for v in 1:length(vars)
    # -- detrending
    new_data = data[v][2:end] - data[v][1:end-1]
    new_data = new_data .- sum(new_data) / length(new_data) # eliminate mean value

    # -- spectrum
    G_v, freqs_v = calc_spectrum(new_data, 1 / dt_out; mode="fft")  # calculates fft or blackman tuckey method

    push!(G_data, G_v ./ (sum(G_v)))  # normalization through the entire integral (only positive values)
    push!(freqs_data, freqs_v)  # save frequencies

    # -- time vs freq map WORK IN PROGRESS -- spm 2022.11.18
    # calc_TimeFreqArr(new_data, window)

end

## Plot
colors = [:orange, :red, :firebrick, :maroon, :blue, :navy, :green, :black, :purple]
plot_spectrum(time, data, freqs_data, G_data, vars, colors, amod_path * "/output/" * experiment * "/" * "amod_spectra.png");
