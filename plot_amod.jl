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
    vars = ["ins", "T_sl", "TMB", "H", "B", "Hsed"]
elseif length(ARGS) == 1
    experiment = ARGS[1]
    vars = ["ins", "TMB", "H", "Hsed"]
else
    experiment = ARGS[1]
    vars = ARGS[2:end]
end

# -- include namelist file of experiment
include("./output/"*experiment*"/namelist.jl")

## Load output data 
data, time = load_nc(amod_path*"/output/"*experiment*"/amod.nc", vars)

## Calculate spectra
G_data, freqs_data, G_array, freqs_array = [], [], [], []
for v in 1:length(vars)
    # -- detrending
    new_data = data[v][2:end] - data[v][1:end-1]
    new_data = new_data .- sum(new_data)/length(new_data)  
    
    # -- spectrum
    F_v, freqs_v = calc_spectrum(new_data, 1 / dt_out)
    G_v = (abs.(F_v).^2) ./2

    push!(G_data, G_v ./ (sum(G_v)))  # normalization through the entire integral (only positive values)
    push!(freqs_data, freqs_v)

    # -- time vs freq map
    # ---- divide the time series in chunks and calculate fft
    G_v_array, freq_v_array = [], []
    idx, window = 1, 50

    for chunk in 1:(0.1*length(new_data))
        new_data_chunk = new_data[idx:(idx+window)]
        F_v_chunk, freqs_v_chunk = calc_spectrum(new_data_chunk, 1 / dt_out)
        G_v_chunk = (abs.(F_v_chunk).^2) ./2

        push!(G_v_array, G_v_chunk ./ (sum(G_v_chunk)))  # normalization through the entire integral (only positive values)
        push!(freq_v_array, freqs_v_chunk)
        idx = idx + 1
    end
    push!(G_array, G_v_array)
    push!(freqs_array, freq_v_array)

end
display(G_array[1][1])

## Plot
colors = [:orange, :red, :blue, :green, :black, :purple]
plot_spectrum(time, data, freqs_data, G_data, vars, colors, amod_path*"/output/"*experiment*"/"*"amod_spectra.png");
