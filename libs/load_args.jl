# =============================
#     Program: load_args.jl
#     Aim: functions to load arguments from run_amod
#     Author: Sergio Pérez-Montero, 2022.11.25
# =============================
function change_namelist(main_path, file, vals)
    cp(main_path*"/par/namelist.jl", out_path*"namelist.txt")
    
end

function load_out(main_path, outfldr)
    # Now check and create outfldr, if necessary
    isdir("./output/"*outfldr) || mkdir("output/"*outfldr)

    output_path = main_path*"/output/"*outfldr*"/"
    return output_path
end

function load_parf(main_path, out_path, file_param)
    # Load parameters and Earth constants
    cp(main_path*"/par/"*file_param, out_path*"namelist.jl", force=true)
    cp(main_path*"/par/earth_const.jl", out_path*"earth_const.jl", force=true)
end

function calc_permutations(data::OrderedDict)
    # calculate lengths
    lengths, data_keys, data_values = [], collect(keys(data)), collect(values(data))
    for (key, value) in data
        push!(lengths, length(value))
    end
    maxlen, minlen, nvars = maximum(lengths), minimum(lengths), length(data)

    # generate a list with the permutations
    perm = copy(data_values[1])
    for i in 2:nvars
        perm = Iterators.product(perm, data_values[i])
    end
    perm = vec(collect(perm))
    perm_list = []
    # AHORA TENGO QUE DESEMPAQUETAR LAS TUPLES

    # group by permutation and variable name
    # new_perm = []
    # for i in 1:nvars:length(perm)-nvars
    #     new_permi = []
    #     for j in 1:nvars
    #         push!(new_permi, perm[i+j])
    #     end
    #     push!(new_perm, new_permi)
    # end

    # perm_keys = []  # -- generate keys
    # for i in 1:perm

    # end

    return perm
end