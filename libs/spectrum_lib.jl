# =============================
#     Program: spectrum_lib.jl
#     Aim: functions to calculate the spectrum of a time series
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================
using FFTW
#using FourierAnalysis

# blackman TUCKEY -- mirar como posible metodo espectral, jas 2022.17.11
@doc """
    calc_spectrum uses FFTW.jl to calculate the periodogram of d 
        d :: Vector --> time series
        t :: Vector --> vector of times
        fs :: Real  --> sampling rate (Hz)

    Note: Algorithm based on https://www.matecdev.com/posts/julia-fft.html
"""
function calc_spectrum(d::Vector, fs::Real; mode="fft")
    if mode == "fft"
        F = fftshift(fft(d))
        freqs = fftshift(fftfreq(length(d), fs))
        F, freqs = F[freqs.>0], freqs[freqs.>0] # now eliminate freqs < 0
        G = (abs.(F) .^ 2) ./ 2 # calculate the power density
    # elseif mode == "fourier"
    #     S = spectra(d, Int(fs * 31536000), Int(fs * 31536000); tapering=rectangular)
    #     display(S.flabels)
    #     F, freqs = S.y, S.flabels
    #     F, freqs = F[freqs.>0], freqs[freqs.>0] # now eliminate freqs < 0
    elseif mode == "tuckey"
        error("tuckey not implemented yet")
    else
        display("mode not recognized, applying fft ...")
        F, freqs = F[freqs.>0], freqs[freqs.>0] # now eliminate freqs < 0
        F = fftshift(fft(d))
        freqs = fftshift(fftfreq(length(d), fs))
        G = (abs.(F) .^ 2) ./ 2 # calculate the power density
    end

    return G, freqs   # power, frequencies
end

@doc """
    calc_TimeFreqArr: calculates the Time-Frequency array of a time series given the window in which Fourier transform is applied
"""
function calc_TimeFreqArr(d::Vector, w::Real)
    nchunks = length(d) / w    # number of chunks in which we divide the time series
    idx = 1 # auxiliar index
    for chunk in 1:nchunks
        d_chunk = d[idx:(2*w+idx-1)]   # each chunk is centered in chunk*window

        idx = 3
    end
end


# C = [480 + 100 * (
#                                              1 * 1 * cos(2.0*pi*t / 41e3) +
#                                              1 * 0.1 * cos(2.0*pi*t / 23e3) + 
#                                              1 * 0.1 * cos(2.0*pi*t / 100e3)) for t in 1:500000] 

# fC, GC = calc_spectrum(C, 1);
# fftshift(1:10)'
# fftfreq(10, 10)