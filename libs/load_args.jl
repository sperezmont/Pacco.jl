# =============================
#     Program: load_args.jl
#     Aim: functions to load arguments from run_amod
#     Author: Sergio Pérez-Montero, 2022.11.25
# =============================
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