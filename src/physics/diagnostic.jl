# =============================
#     Program: diagnostic.jl
#     Aim: This program contains functions to calculate diagnostic variables
# =============================
"""
    calc_E(now, par)
calculates ice-sheet extension

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_E(now, par)
    for hm in par["hemisphere"]
        if par["active_climate"]
            Temp = now["T_"*hm]
        else    
            Temp = now["T_sl_"*hm]
        end
        now["E_"*hm] = par["E_ref_"*hm] * (Temp - par["T_ref0_"*hm]) / par["A_te_"*hm]
    end
    return now
end

"""
    calc_V(now, par) 
calculates ice-sheet volume from prognostic mean H and diagnosed E

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_V(now, par)
    for hm in par["hemisphere"]
        now["V_"*hm] = now["E_"*hm] * now["H_"*hm] / 1e3   # km3
        now["V_"*hm] = now["V_"*hm] * rhoi / rhow * 1e3 / A_oc    # m sle
    end
    return now
end





