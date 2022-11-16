# =============================
#     Program: spectrum_lib.jl
#     Aim: functions to calculate the spectrum of a time series
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================
using FFTW 

@doc """
    calc_spectrum uses FFTW.jl to calculate the periodogram of d 
        d :: Vector --> time series
        t :: Vector --> vector of times
        fs :: Real  --> sampling rate (Hz)

    Note: Algorithm based on https://www.matecdev.com/posts/julia-fft.html
"""
function calc_spectrum(d::Vector, fs::Real) 
    F = fftshift(fft(d))
    freqs = fftshift(fftfreq(length(d), fs))  

    return (abs.(F).^2) ./2, freqs   # power, frequencies
end

# C = [480 + 100 * (
#                                              1 * 1 * cos(2.0*pi*t / 41e3) +
#                                              1 * 0.1 * cos(2.0*pi*t / 23e3) + 
#                                              1 * 0.1 * cos(2.0*pi*t / 100e3)) for t in 1:500000] 

# fC, GC = calc_spectrum(C, 1);
# fftshift(1:10)'
# fftfreq(10, 10)