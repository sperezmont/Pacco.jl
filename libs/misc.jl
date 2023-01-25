# =============================
#     Program: misc.jl
#     Aim: multiple functions for different tasks 
#     Author: Sergio Pérez-Montero, 2022.11.25
# =============================

@doc """
    change_namelist
        Takes a dictionary of new values for parameters and change them in the namelist.jl file
        Adapted from https://stackoverflow.com/questions/58013970/how-to-edit-a-line-of-a-file-in-julia
"""
function change_namelist(outpath, file, params)
    (tmppath, tmpio) = mktemp(outpath)  # create temporary files
    open(outpath * "/" * file) do io
        for line in eachline(io, keep=true) # keep so the new line isn't chomped
            for param in params.keys
                if occursin(param, line)
                    splitted_line = split(line) # -- split in elements
                    splitted_line[3] = string(params[param]) # -- modify third element (value)
                    line = join(splitted_line, " ") * "\n" # -- rewrite line 
                end
            end
            write(tmpio, line)
        end
    end
    close(tmpio)
    mv(tmppath, outpath * "/" * file, force=true)   # rewrite param file
end

@doc """
    load_out
        Loads the name of output directory and creates it (if necessary)
"""
function load_out(main_path, outfldr)
    # Now check and create outfldr, if necessary
    isdir(main_path * "/output/" * outfldr) || mkdir(main_path * "/output/" * outfldr)

    output_path = main_path * "/output/" * outfldr * "/"
    return output_path
end

@doc """
    load_parf
        Copy important par/ files to output directory
"""
function load_parf(main_path, out_path, file_param)
    # Load parameters and Earth constants
    cp(main_path * "/par/" * file_param, out_path * "namelist.jl", force=true)
    cp(main_path * "/par/earth_const.jl", out_path * "earth_const.jl", force=true)
end

@doc """
    calc_comb
        takes a list of vectors and returns the number of possible combinations
"""
function calc_comb(d::Vector)
    # returns number of combinations between elements of d
    n = 1
    for i in 1:length(d)
        n = n * length(d[i])
    end
    return n
end

@doc """
    calc_permutations  
        Takes a dictionary of variables => [list of values] and
        returns a list of dictionaries with all the possible permutations
"""
function calc_permutations(d::OrderedDict)
    # calculate lengths
    lengths, d_keys, d_values = [], collect(keys(d)), collect(values(d))
    for (key, value) in d
        push!(lengths, length(value))
    end
    maxlen, minlen, nvars = maximum(lengths), minimum(lengths), length(d)

    # generate a list with the permutations
    ncomb = calc_comb(d_values)
    list_of_perm = Array{Any}(undef, (ncomb, nvars)) # ncomb, nvariables
    for i in 1:length(d_values)
        ncombi = calc_comb(d_values[i:end])
        list_of_perm[:, i] = repeat(d_values[i], inner=Int(ncombi / length(d_values[i])), outer=Int(ncomb / ncombi))
    end
    list_of_perm = [list_of_perm[i, :] for i in 1:size(list_of_perm)[1]] # transform to list of lists

    # create the list of dictionaries
    perm = []
    for i in 1:ncomb
        push!(perm, OrderedDict(d_keys .=> list_of_perm[i]))
    end
    return perm
end

# @doc """
#     calc_ndn:
#         calculates the nearest 2^n number to nmbr by Jan Swierczek-Jereczek
# """
# function calc_ndn(nmbr::Int)
#     dyadic_list = 2.0 .^ (1:50)
#     n = argmin((Float64(nmbr) .- dyadic_list) .^ 2.0)
#     return 2^n
# end

# @doc """
#     vector2dyadic:
#         Takes a vector and transforms it to dyadic size
# """
# function vector2dyadic(d)
#     # define ld
#     ld = length(d)

#     # create grid and interpolant
#     grid = 1:ld
#     d_interp = linear_interpolation(grid, d)

#     # calculate nearest dyadic size
#     nds = calc_ndn(ld)
#     grid_dyadic = range(1, stop=ld, length=nds)

#     # return new time series
#     return d_interp.(grid_dyadic)
# end

@doc """


"""
function calc_spectrum(d, fs)
    # -- spectrum blackman tuckey
    N = length(d)
    Nmax = Int(ceil(N / 2))
    P = periodogram(d, fs=fs, window=blackman(N))
    G, freq = P.power, P.freq
    return G, freq
end

@doc """
    calc_wavelet:
        Generate an array with the values of the wavelet applied to d
"""
function calc_wavelet(d, fs; sigma=π)
    wvt = ContinuousWavelets.wavelet(Morlet(sigma), s=8, boundary=ZPBoundary(), averagingType=NoAve(), β=1)
    S = ContinuousWavelets.cwt(d, wvt)                                           
    freq= getMeanFreq(ContinuousWavelets.computeWavelets(length(d), wvt)[1], fs)  
    S = abs.(S) .^ 2 
    S = S ./ sum(S, dims=2)

    # # -- cone of influence
    # s = 8
    # taus = sqrt(2) * s
    # coi = 1 ./ exp(t ./ taus)
    return S, freq
end

@doc """
    gen_lhs:
        Generate a Latin Hypercube Sampling for AMOD ensembles using LatinHypercubeSampling.jl
        par2per --> Dictionary with parameters to permute, key => (min, max)
        nsim    --> number of simulations (number of sample points)
        ngns -> number of generations

        Based on example from:
            https://mrurq.github.io/LatinHypercubeSampling.jl/stable/man/lhcoptim/
"""
function gen_lhs(par2per::Dict, nsim::Int; ngens=10)
    # Plan the LHS
    nds = length(par2per)
    plan, _ = LHCoptim(nsim, nds, ngens)

    # Generate a vector with [min, max] values
    minmax = [par2per[i] for i in keys(par2per)]

    # Scale plan
    scaled_plan = scaleLHC(plan, minmax)
    scaled_plan_list = [scaled_plan[i, :] for i in 1:size(scaled_plan)[1]]       # -- convert to list of lists
    scaled_plan_dict = [Dict(keys(par2per) .=> scaled_plan_list[i]) for i in 1:length(scaled_plan_list)] # -- convert to dictionary

    return scaled_plan, scaled_plan_dict
end