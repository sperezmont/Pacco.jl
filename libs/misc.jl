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
            for param in keys(params)
                if occursin(param*" ", line)    # this space is needed!
                    splitted_line = split(line) # -- split in elements
                    if typeof(params[param]) == String
                        splitted_line[3] = "\""*string(params[param])*"\""
                    else
                        splitted_line[3] = string(params[param]) # -- modify third element (value)
                    end
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
    unpack_lhs:
        takes a matrix row and generate a dictionary with keys(d)=orig_keys and typeof(d[key]) = types
"""
function unpack_lhs(m, row, types, orig_keys)
    mpos = m[row, :]
    d = []
    for j in eachindex(mpos)
        push!(d, types[j](mpos[j]))
    end
    d = Dict(orig_keys .=> d)
    return d
end

@doc """
    gen_lhs:
        Generate a Latin Hypercube Sampling for PACCO ensembles using LatinHypercubeSampling.jl
        par2per     --> Dictionary with parameters to permute, key => (min, max)
        nsim        --> number of simulations (number of sample points)
        ngns        --> number of generations
        pars_type   --> mode of operation: 1 = all continuous, 2 = mix, 3 = read pars_type_list
        pars_type_list --> List of parameters type: Continuous() or Categorical(ncatvals, catWeight)

        Based on example from:
            https://mrurq.github.io/LatinHypercubeSampling.jl/stable/man/lhcoptim/
"""
function gen_lhs(par2per::Dict, nsim::Int; ngens=10, pars_type=1, pars_type_list=[], ncatvals=2, catWeigth=0.0025)
    # Generate a vector with [min, max] values
    minmax = [par2per[i] for i in keys(par2per)]

    # Get types of parameters
    types = [typeof.(par2per[key][1]) for (key, value) in par2per]

    if pars_type == 1  # assume all parameters are Continuous()        
        # Plan the LHS
        nds = length(par2per)
        plan, _ = LHCoptim(nsim, nds, ngens)
        types = [Float64 for (key, value) in par2per]
    else    # assume continuous and boolean categorical parameters
        if pars_type == 2
            dims::Vector{LHCDimension} = []
            for (ky, vl) in par2per
                if (vl[1] in [true, false]) && (vl[2] in [true, false])
                    ncatvals = 2
                    push!(dims, Categorical(ncatvals, catWeigth))
                else
                    push!(dims, Continuous())
                end
            end 
        else
            dims = pars_type_list
        end
        # Plan the LHS
        nds = length(par2per)
        initialSample = randomLHC(nsim, dims)
        plan = LHCoptim!(initialSample, ngens; dims=dims)[1]     
    end

    # Scale plan
    scaled_plan = scaleLHC(plan, minmax)
    scaled_plan_dict = [unpack_lhs(scaled_plan, i, types, keys(par2per)) for i in 1:size(scaled_plan)[1]]
    
    return scaled_plan, scaled_plan_dict
end
