# =============================
#     Program: amod_diag.jl
#     Aim: This program contains functions to calculate diagnostic variables
# =============================
@doc """
    calc_E: calculates ice extension
"""
function calc_E(now_di, par_di)
    for hm in par_di["hemisphere"]
        now_di["E"*hm] = par_di["E_ref"*hm] * (now_di["T"*hm] - par_di["T_ref"*hm]) / par_di["A_te"*hm]
    end
    return now_di
end

@doc """
    calc_V: calculates ice volume from predicted mean H and diagnosed E
"""
function calc_V(now_di, par_di)
    for hm in par_di["hemisphere"]
        now_di["V"*hm] = now_di["E"*hm] * now_di["H"*hm] / par_di["L"*hm]
    end
    return now_di
end





