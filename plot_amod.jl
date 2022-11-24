# =============================
#     Program: plot_amod.jl
#     Aim: Plot AMOD output
#     How to run: julia plot_amod.jl [OUT] [VARS]
#         where
#             [OUT]   --> output name (without "output/", default=test_default)
#             [VARS]  --> variables to plot
# =============================

## First, load external libraries, activate amod environment and save pwd 
using Pkg
Pkg.activate("amod_env")

using DSP
using FFTW
using Wavelets
using Statistics

amod_path = pwd()

include("./libs/nc.jl")
include("./libs/plot_lib.jl")

## Check arguments
default_vars = ["ins_norm", "T_sl", "SMB", "H", "Hsed"]
if ARGS == []
    experiment = "test_default"
    vars = default_vars
elseif length(ARGS) == 1
    experiment = ARGS[1]
    vars = default_vars
else
    experiment = ARGS[1]
    vars = ARGS[2:end]
end

## Load output data 
data, time = load_nc(amod_path * "/output/" * experiment * "/amod.nc", vars)

## Calculate spectra
#window = 50
G_data, freqs_data = [], []
for v in 1:length(vars)
    new_data = copy(data[v])

    # -- spectrum
    if false
        # ---- fft
        F, N = fft(new_data), length(new_data)
        Nmax = Int(ceil(N/2))
        kvec = 1:Nmax
        G = abs.(F)*2/N
        G = sqrt.(G[1:Nmax])
        freq = kvec ./ (N*(time[2] - time[1])) 
        G, freq = G[freq .>= 1/150.5e3], freq[freq .>= 1/150e3] # we eliminate values above and below Milankovitch cycles
        G, freq = G[freq .<= 1/22e3], freq[freq .<= 1/22e3]
        G = G ./ sum(G)
        if var(new_data) < 1e-1
            G = zeros(length(G))
        end
    else
        # ---- blackman tuckey
        N = length(new_data)
        Nmax = Int(ceil(N/2))
        P = periodogram(new_data, onesided=false, fs=1/(time[2] - time[1]), window=blackman(N))
        G, freq = P.power, P.freq
        G, freq = G[1:Nmax] .* 2, freq[1:Nmax]
        G = sqrt.(G) / (N/2)
        G, freq = G[freq .>= 1/150.5e3], freq[freq .>= 1/150e3] # we eliminate values above and below Milankovitch cycles
        G, freq = G[freq .<= 1/22e3], freq[freq .<= 1/22e3]
        G = G ./ sum(G)
        if var(new_data) < 1e-1
            G = zeros(length(G))
        end
    end
    push!(G_data, G)
    push!(freqs_data, freq)  # save frequencies

    # -- time vs freq map WORK IN PROGRESS -- spm 2022.11.18
    # calc_TimeFreqArr(new_data, window)

end

## Plot
clrmp_dict = Dict("time"=>:grays1,
    "H" => cgrad([:royalblue, :royalblue4]),#cgrad([:snow4, :royalblue, :royalblue4, :navy]),
    "Hsed" => cgrad([:firebrick, :maroon]),#cgrad([:maroon, :black]),
    "T" => :lajolla,
    "A" => :lajolla,
    "T_sl" => cgrad([:purple, :purple]),
    "TMB" => cgrad(:seismic, [0.0, 0.45, 0.55, 1.0], rev = true),
    "SMB" => cgrad(:seismic, [0.0, 0.45, 0.55, 1.0], rev = true),
    "Z" => :dense,
    "B" => :dense,
    "M" => :amp,
    "Acc" => cgrad(:ice, rev = true),
    "U_d" => :corkO,
    "U_b" => :corkO,
    "U" => :corkO,
    "T_surf" => :lajolla,
    "tau_b" => :thermal,
    "tau_d" => :thermal,
    "Q_dif" => :lajolla,
    "Q_difup" => :lajolla,
    "Q_difdown" => :lajolla,
    "Q_drag" => :lajolla,
    "alpha" => :lajolla,
    "Q_adv" => :lajolla,
    "fstream" => :grays1,
    "fstream_ref" => :grays1,
    "Hdot" => :ice,
    "Hseddot" => :copper,
    "Bdot" => :dense,
    "Tdot" => :lajolla,
    "fstreamdot" => :grays1,
    "ins" => cgrad([:black, :black]),#cgrad(:seismic, [0.0, 0.49, 0.51, 1.0]),
    "ins_norm" => cgrad([:black, :black]),#cgrad(:seismic, [0.0, 0.49, 0.51, 1.0]),
    "co2" => :thermal,
    "P" => :dense,
    "exc" => :thermal,
    "long_peri" => :thermal,
    "obl" => :thermal
)

colors = [clrmp_dict[i] for i in vars]
plot_spectrum(time, data, freqs_data, G_data, vars, colors, amod_path * "/output/" * experiment * "/" * "amod_results.png");
