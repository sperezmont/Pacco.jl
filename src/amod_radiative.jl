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
    now_r = calc_insol_day(now_r, par_r)
    insnorm = (now_r["ins"] - par_r["ins_min"]) / (par_r["ins_max"] - par_r["ins_min"])     # between 0 and 1, norm = 1
    now_r["ins_norm"] = 2.0 * insnorm - 1.0                                                          # between 1 and -1, norm = 2

    # Third, compute anthropogenic forcing if time >= +2000 yrs
    if (now_r["time"] >= par_r["time_ant"]) 
        T_ant = par_r["A_ant"] / exp((now_r["time"] - par_r["time_ant"]) / par_r["tau_ant"])
    else
        T_ant = 0.0
    end

    # Finally, calculate sea-level temperature 
    now_r["T_sl"] = T_ref + par_r["A_t"] * now_r["ins_norm"] + T_ant
    return now_r
end




