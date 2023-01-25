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
    calc_albedo_ref: calculates reference value for albedo
"""
function calc_albedo_ref(now_r, par_r)
    for hm in par_r["hemisphere"]
        if now_r["H_"*hm] == 0.0
            now_r["ice_time_"*hm] = 0.0 # -- no ice
            now_r["albedo_ref_"*hm] = par_r["albedo_land"]
        else
            n = 1 # exponent number in albedo-age parameterisation
            now_r["albedo_ref_"*hm] = max(par_r["albedo_newice"] - par_r["albedo_slope"] * now_r["ice_time_"*hm]^n, par_r["albedo_land"])
        end
    end
    return now_r
end

@doc """
    calc_rf: calculates radiative forcing for each hemisphere
"""
function calc_rf(now_r, par_r)
    for hm in par_r["hemisphere"]
        # -- insolation
        now_r = calc_ins(now_r, par_r)
        now_r["ins_anom_"*hm] = now_r["ins_"*hm] - par_r["ins_ref_"*hm]
        # -- radiative forcing from co2
        radco2 = calc_rad_co2(now_r["co2_"*hm])
        # -- total        
        now_r["rf_"*hm] = now_r["ins_anom_"*hm] + radco2

        #if (now_r["time"] >= par_r["time_anth"]) # -- anthropogenic forcing
        #    now_r["rf_"*hm] += par_r["Ac_anth"] / exp((now_r["time"] - par_r["time_anth"]) / par_r["tau_anth"])
        #end
    end
    return now_r
end

@doc """
    calc_Tsl: calculates sea level temperature
"""
function calc_Tsl(now_r, par_r)
    # First, calculate insolation
    now_r = calc_ins(now_r, par_r)

    for hm in par_r["hemisphere"]
        # -- normalize insolation
        insnorm = (now_r["ins_"*hm] - par_r["ins_min"]) / (par_r["ins_max"] - par_r["ins_min"])     # between 0 and 1, norm = 1
        now_r["ins_norm_"*hm] = 2.0 * insnorm - 1.0                                                          # between 1 and -1, norm = 2

        # Second, compute anthropogenic forcing if time >= +2000 yrs
        if (now_r["time"] >= par_r["time_anth"])
            T_anth = par_r["At_anth"] / exp((now_r["time"] - par_r["time_anth"]) / par_r["tau_anth"])
        else
            T_anth = 0.0
        end

        # Finally, calculate sea-level temperature 
        now_r["T_sl_"*hm] = par_r["T_ref_"*hm] + par_r["A_t"] * now_r["ins_norm_"*hm] + T_anth
    end
    return now_r
end




