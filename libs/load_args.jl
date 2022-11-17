# =============================
#     Program: load_args.jl
#     Aim: functions to load arguments from the command line
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

function load_out(main_path, args)
    # Now check and create [OUT], if necessary
    (args == []) ? (outfldr = "test_default") : (outfldr = args[1])
    isdir("./output/"*outfldr) || mkdir("output/"*outfldr)

    output_path = main_path*"/output/"*outfldr*"/"
    return output_path, outfldr
end

function load_parf(main_path, out_path, args)
    # Load parameters and Earth constants
    if (args == []) || (length(args) == 1)
        file_param = "amod_default.jl"
    else 
        file_param = args[2]
    end

    cp(main_path*"/par/"*file_param, output_path*"namelist.jl", force=true)
    cp(main_path*"/par/earth_const.jl", output_path*"earth_const.jl", force=true)

    include(out_path*"namelist.jl")
    include(out_path*"earth_const.jl")
end