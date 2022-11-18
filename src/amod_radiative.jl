# =============================
#     Program: amod_radiative.jl
#     Aim: This program contains functions to calculate radiative parameters
# =============================

@doc """
"""
function calc_rad_co2(CO2)
    CO2_0 = 280.0
    RCO2_fac = 5.35
    return RCO2_fac * log(CO2 / CO2_0)
end

@doc """
"""
function calc_Tsl(now_r, par_r)
    # First, calculate radiative contribution of current CO2 level
    if par_r["active_radco2"]
        now_r["co2"] = calc_rad_co2(now_r["co2"])     # it does nothing for the moment -- spm 2022.11.08
    end

    # Second, calculate insolation and normalize it
    now_r["ins"] = calc_insol_day(now_r, par_r)
    ins_norm = (now_r["ins"] - par_r["ins_min"]) / (par_r["ins_max"] - par_r["ins_min"])     # between 0 and 1, norm = 1
    ins_norm = 2.0 * ins_norm - 1.0                                                          # between 1 and -1, norm = 2

    # Third, calculate sea-level temperature 
    now_r["T_sl"] = T_ref + par_r["A_t"] * ins_norm
    return now_r
end




