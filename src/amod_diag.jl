# =============================
#     Program: amod_diag.jl
#     Aim: This program contains functions to calculate diagnostic variables
# =============================
@doc """
    calc_E: calculates ice extension
"""
function calc_E(now_di, par_di)
    for hm in par_di["hemisphere"]
        if par_di["active_climate"]
            Temp = now_di["T_"*hm]
        else    
            Temp = now_di["T_sl_"*hm]
        end
        now_di["E_"*hm] = par_di["E_ref_"*hm] * (Temp - par_di["T_ref_"*hm]) / par_di["A_te_"*hm]
    end
    return now_di
end

@doc """
    calc_V: calculates ice volume from predicted mean H and diagnosed E
"""
function calc_V(now_di, par_di)
    for hm in par_di["hemisphere"]
        now_di["V_"*hm] = now_di["E_"*hm] * now_di["H_"*hm] / 1e3   # km3
        now_di["V_"*hm] = now_di["V_"*hm] * rhoi / rhow * 1e3 / A_oc    # m sle
    end
    return now_di
end





