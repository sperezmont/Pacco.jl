# =============================
#     Program: load_args.jl
#     Aim: functions to load arguments from different functions
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
function calc_comb(data::Vector)
    # returns number of combinations between elements of data
    n = 1
    for i in 1:length(data)
        n = n * length(data[i])
    end
    return n
end

@doc """
    calc_permutations  
        Takes a dictionary of variables => [list of values] and
        returns a list of dictionaries with all the possible permutations
"""
function calc_permutations(data::OrderedDict)
    # calculate lengths
    lengths, data_keys, data_values = [], collect(keys(data)), collect(values(data))
    for (key, value) in data
        push!(lengths, length(value))
    end
    maxlen, minlen, nvars = maximum(lengths), minimum(lengths), length(data)

    # generate a list with the permutations
    ncomb = calc_comb(data_values)
    list_of_perm = Array{Any}(undef, (ncomb, nvars)) # ncomb, nvariables
    for i in 1:length(data_values)
        ncombi = calc_comb(data_values[i:end])
        list_of_perm[:, i] = repeat(data_values[i], inner=Int(ncombi / length(data_values[i])), outer=Int(ncomb / ncombi))
    end
    list_of_perm = [list_of_perm[i, :] for i in 1:size(list_of_perm)[1]] # transform to list of lists

    # create the list of dictionaries
    perm = []
    for i in 1:ncomb
        push!(perm, OrderedDict(data_keys .=> list_of_perm[i]))
    end
    return perm
end


# A = [[1, 2, 3], [20, 30], [1000, 2000, 3000, 4000], [3, 15]]
# ncomb = 1
# for i in 1:length(A)
#     ncomb = ncomb * length(A[i])
# end
# B = Array{Any}(undef, (ncomb, length(A))) # ncomb, nvariables

# for i in 1:length(A)
#     ncombi = 1
#     for j in i:length(A)
#         ncombi = ncombi * length(A[j])
#     end
#     B[:, i] = repeat(A[i], inner=Int(ncombi / length(A[i])), outer=Int(ncomb / ncombi))
# end