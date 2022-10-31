# =============================
#     Program: load_args.jl
#     Aim: functions to load arguments from the command line
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
    (args == []) ? (file_param = "amod_default.jl") : (file_param = args[2])

    cp(main_path*"/par/"*file_param, output_path*file_param, force=true)
    cp(main_path*"/par/earth_const.jl", output_path*"earth_const.jl", force=true)

    include(out_path*file_param)
    include(out_path*"earth_const.jl")
end