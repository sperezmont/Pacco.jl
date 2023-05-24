# =============================
#     Program: misc.jl
#     Aim: multiple functions for different tasks 
#     Author: Sergio Pérez-Montero, 2022.11.25
# =============================

"""
    write_run_info(outpath, expname, params)
writes the running information in `run-info.txt`
"""
function write_run_info(outpath::String, expname::String, params)
    isfile(outpath * "/run-info.txt") && rm(outpath * "/run-info.txt")
    f = open(outpath * "/run-info.txt", "w")
    write(f, "$(expname)\n")
    write(f, "# Parameters used: \n")

    nms = fieldnames(typeof(params))
    for i in eachindex(nms)
        valuei = getfield(params, nms[i])
        write(f, "  $(nms[i]) --> $(valuei)\n")
    end
    close(f)
end

"""
    get_index_from_time(time, times2find; tol)
finds the index required in time series

## Attributes
* time        --> time vector
* times2find  --> times to find in time vector (ONLY two elements)
"""
function get_index_from_time(time::Vector, times2find::Vector; tol = 1e-5)
    idx1, idx2 = 1, 2

    for i in eachindex(time)
        error = abs.(time[i] .- times2find) ./ abs.(times2find .+ tol)  # to avoid 1/0

        if error[1] < tol
            idx1 = i
        end
        if error[2] < tol
            idx2 = i
        end
    end
    return (idx1, idx2)
end

"""
    get_runs(experiment)
gets the experiment or experiments given a string or a vector of strings
"""
function get_runs(experiment)
    if typeof(experiment) == String
        if split(experiment, "/")[1] == "output"
            out_path = pwd() * "/" * experiment * "/"
            new_experiment = join(split(experiment, "/")[2:end], "/")
        else
            out_path = pwd() * "/output/" * experiment * "/"
            new_experiment = experiment
        end

        if isfile(out_path * "/pacco.nc")   # one run
            labels = [new_experiment]
            desired_runs = out_path .* ["/pacco.nc"]
        else    # ensemble of runs
            elements = readdir(out_path .* "/runs/")
            desired_runs = out_path .* "/runs/" .* elements .* "/pacco.nc"
            labels = [split.(desired_runs[i], "/")[end-1] for i in eachindex(desired_runs)]
        end
    elseif typeof(experiment) == Vector{String}
        out_path = pwd() .* "/output/" .* experiment .* "/"
        labels = experiment
        desired_runs = out_path .* "/pacco.nc"
        new_experiment = experiment
    end
    return desired_runs, labels, new_experiment
end

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
function gen_lhs(par2per::Dict, nsim::Int; ngens=5, pars_type=1, pars_type_list=[], ncatvals=2, catWeigth=0.0025)
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
function calc_permutations(d::Dict)
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
        push!(perm, Dict(d_keys .=> list_of_perm[i]))
    end
    return perm
end

##########################################################################################################
# OLD CODE
##########################################################################################################

"""
    get_parameters_from_good_runs(experiment, filename)
returns the parameter names and values of an ensemble
"""
function get_parameters_from_runs(experiment, filename; all_runs=true)
    path_to_results = get_path_to_results_or_runs(experiment, "results")

    if all_runs
        runs = readlines(path_to_results * "/$(filename)") .* "/pacco.nc"
    else
        runs = get_good_runs_from_file(path_to_results * "/$(filename)") .* "/pacco.nc"
    end

    length_to_remove = length("/pacco.nc")
    vector_of_run_info = []
    for e in eachindex(runs)
        if runs[e][1:3] == "NS_"
            new_name = runs[e][4:end]
            new_name = new_name[1:end-length_to_remove]
            push!(vector_of_run_info, get_path_to_results_or_runs(experiment, "runs") * "/$(new_name)/run-info.txt")
        else
            push!(vector_of_run_info, runs[e][1:end-length_to_remove] * "/run-info.txt")
        end
    end
    par_names, parameters = [], []
    for i in eachindex(vector_of_run_info)
        open(vector_of_run_info[i]) do file
            lines = readlines(file)
            for line in lines
                if line[1:11] == "OrderedDict"
                    new_line = split(line[13:end-1], " ")
                    pars = []

                    for e in new_line[1:3:end]
                        if i == 1
                            push!(par_names, e[2:end-1])
                        end
                    end
                    
                    for e in new_line[3:3:end]
                        if e[end] == ','
                            e2save = e[1:end-1]
                        else
                            e2save = e[1:end]
                        end
                        push!(pars, parse(Float64, e2save))
                    end
                    push!(parameters, pars)
                end
            end
        end
    end
    return par_names, parameters
end

"""
    detect_deglaciation(x, dt)
detects complete deglaciations in a time series

## Arguments
* `x` time series (> 0) to analyze (e.g ice thickness or ice volume)
* `t` time dimension

## Return
`number_of_deglaciations` and `timing`
"""
function detect_deglaciation(x::Vector, t::Vector)
    dx = x[2:end] .- x[1:end-1]
    dt = t[2:end] .- t[1:end-1]

    dxdt = (dx ./ dt)

    is_deglaciation = zeros(length(dxdt))
    for i in eachindex(dxdt)
        if x[i+1] == 0.0    # no ice
            if sign(dxdt[i]) < 0   # glacial termination found
                is_deglaciation[i] = true
            end
        else
            is_deglaciation[i] = false
        end
    end
    number_of_deglaciations = sum(is_deglaciation)
    timing = vcat(false, is_deglaciation) # in this way, we extend to the same length as original data
    return number_of_deglaciations, timing
end

function read_dir_good_runs(outpath)
    elements = readdir(outpath)
    nchars, k = length("good_runs"), 1
    only_txt_files = []
    for e in elements
        if e[1:nchars] == "good_runs"
            push!(only_txt_files, e)
        end
    end
    return only_txt_files
end

function write_filtered_ensemble(outpath, mask, rule, file)
    if file != []
        only_txt_files = read_dir_good_runs(outpath)
        if only_txt_files == []
            file_name = "good_runs_r$(rule).txt"
        else
            file_name = "$(file[1:end-4])r$(rule).txt"
        end
    else
        file_name = "good_runs_r$(rule).txt"
    end
    
    f = open(outpath * "/$(file_name)", "w")
    path_to_runs = join(split(outpath, "/")[1:end-2], "/") * "/runs/"
    runs = readdir(path_to_runs)
    for i in eachindex(mask)
        if mask[i] == true
            write(f, "$(path_to_runs)$(runs[i])\n")
        else
            write(f, "NS_$(runs[i])\n") # NOT SELECTED RUN
        end
    end
    close(f)
end

function get_exp_labels(runs)
    labels = []
    for run in runs
        elements = split(run, "/")
        if run[1:3] == "NS_"
            push!(labels, elements[end-1][4:end])
        else
            push!(labels, elements[end-1])
        end
    end
    return labels     
end

function get_good_runs_from_file(path_to_file)
    runs, good_runs = readlines(path_to_file), []
    for i in eachindex(runs)
        if runs[i][1:3] != "NS_"
            push!(good_runs, runs[i])
        end
    end
    return good_runs
end

function get_path_to_results_or_runs(experiment, dir)
    exp_path = pwd() * "/output/" * experiment
    elements = readdir(exp_path)
    if dir in elements
        return exp_path * "/$(dir)/"
    else
        return exp_path
    end 
end

function is_ensemble(exp)
    elems = readdir(pwd() * "/output/" * exp * "/")
    if "pacco.nc" in elems
        is_ens = false
        new_new_elems = elems
    else
        is_ens = true

        new_elems = readdir(pwd() * "/output/" * exp * "/runs/")
        new_new_elems = []
        for e in new_elems # drop "/results/" directory
            push!(new_new_elems, e)
        end
    end
    return new_new_elems, is_ens
end

@doc """
    readdir_store_sims:
"""
function readdir_store_sims(path_to_experiment)
    elements = []
    for element in readdir(path_to_experiment)
        if element[end-3:end] != ".png"
            push!(elements, element)
        end
    end
    return elements
end


@doc """
    nans_detector:
        This function will scream if any nan or missing is found in the ensemble runs
"""
function nans_detector(; experiment::String="test_default_ens", variable="H_n")
    # Define some local variables
    locdir = pwd() * "/output/" * experiment * "/"

    # Read directory
    elements = readdir_store_sims(locdir)

    # Look for nans and scream their names, times and variables!
    k = 0
    for e in elements
        d = NCDataset(locdir * "/" * e * "/pacco.nc")[variable]
        isnan_list, ismissing_list = isnan.(d), ismissing.(d)
        if sum(isnan_list) > 0
            printstyled("NaN found in $e", color=:red)
            k += 1
        end
        if sum(ismissing_list) > 0
            printstyled("missing found in $e", color=:red)
            k += 1
        end
    end
    if k == 0
        printstyled("no NaN or missing found", color=:green)
    end
end

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

