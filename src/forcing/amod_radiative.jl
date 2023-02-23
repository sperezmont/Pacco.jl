# =============================
#     Program: amod_radiative.jl
#     Aim: This program contains functions to calculate radiative parameters
# =============================
@doc """
    calc_rad_co2:
        calculates the radiative forcing due to co2 (in W/m²)
"""
function calc_rad_co2(CO2::Real)
    CO2_0 = 280.0
    RCO2_fac = 5.35
    return RCO2_fac * log(CO2 / CO2_0)
end

@doc """
    calc_albedo_ref:
        calculates reference value for albedo
"""
function calc_albedo_ref(now_r::OrderedDict, par_r::OrderedDict)
    for hm in par_r["hemisphere"]
        if now_r["H_"*hm] == 0.0
            now_r["ice_time_"*hm] = 0.0 # -- no ice
            now_r["albedo_ref_"*hm] = par_r["albedo_land"]
        else
            n = 1 # exponent number in albedo-age parameterisation
            now_r["albedo_ref_"*hm] = max(par_r["albedo_newice"] - par_r["albedo_slope"] * now_r["ice_time_"*hm]^n, par_r["albedo_oldice"])
        end
    end
    return now_r
end

@doc """
    calc_T_rf:
        calculates radiative forcing for each hemisphere
"""
function calc_T_rf(now_r::OrderedDict, par_r::OrderedDict)
    for hm in par_r["hemisphere"]
        # -- insolation
        now_r = calc_ins(now_r, par_r)
        now_r["ins_anom_"*hm] = now_r["ins_"*hm] - par_r["ins_ref_"*hm]
        # -- radiative forcing from co2
        radco2 = calc_rad_co2(now_r["co2_"*hm])
        # -- total        
        now_r["T_rf_"*hm] = par_r["csi"] * now_r["ins_anom_"*hm] + par_r["cs"] * radco2

    end
    return now_r
end

@doc """
    calc_T_sl:
        calculates sea level temperature
"""
function calc_T_sl(now_r::OrderedDict, par_r::OrderedDict)
    # First, calculate insolation
    now_r = calc_ins(now_r, par_r)

    for hm in par_r["hemisphere"]
        # -- normalize insolation
        insnorm = (now_r["ins_"*hm] - par_r["ins_min"]) / (par_r["ins_max"] - par_r["ins_min"])     # between 0 and 1, norm = 1
        now_r["ins_norm_"*hm] = 2.0 * insnorm - 1.0                                                          # between 1 and -1, norm = 2

        # Second, compute anthropogenic forcing if time >= time_anth
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

#############################
# Time derivatives
#############################
@doc """
    calc_Tdot:
        calculates regional temperature derivative
"""
function calc_Tdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]

        # Anthropogenic forcing?
        if now_dt["time"] < par_dt["time_anth"] # unperturbed climate
            temp_ref = par_dt["T_ref_"*hm]
        else    # perturbed climate
            temp_ref = par_dt["T_ref_"*hm] + par_dt["cco2"] * now_dt["co2_"*hm] / par_dt["co2_ref"]
        end

        if par_dt["height_temp"] == "useH"
            HTF = now_dt["H_"*hm]
        elseif par_dt["height_temp"] == "useZ"
            HTF = now_dt["Z_"*hm]
        else
            printstyled("dev par must be removed!", color=:red)
        end

        now_dt["Tdot_"*hm] = (temp_ref - now_dt["T_"*hm]
                              +
                              now_dt["T_rf_"*hm]
                              -
                              par_dt["csz"] * HTF) / par_dt["tau_rf_"*hm]
    end
    return now_dt
end

@doc """
    calc_albedodot:
        calculates albedo derivative
"""
function calc_albedodot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now_dt["albedodot_"*hm] = (now_dt["albedo_ref_"*hm] - now_dt["albedo_"*hm]) / par_dt["tau_albedo"]

        if now_dt["H_"*hm] == 0.0
            now_dt["albedo_"*hm] = par_dt["albedo_land"]
            now_dt["albedodot_"*hm] = 0.0
        elseif now_dt["ice_time_"*hm] < 10.0 # First ice
            now_dt["albedo_"*hm] = par_dt["albedo_newice"]
            now_dt["albedodot_"*hm] = 0.0
        end
    end
    return now_dt
end

@doc """
    calc_co2dot:
        calculates co2 derivative
"""
function calc_co2dot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]

        # -- anthropogenic forcing?
        if now_dt["time"] < par_dt["time_anth"] # unperturbed climate
            co2ref = par_dt["co2_ref"]
        else    # perturbed climate
            actual_diftime = now_dt["time"] - par_dt["time_anth"]
            co2ref = par_dt["co2_ref"] + par_dt["co2_anth"] * (
                0.75 / exp(actual_diftime / 365.0)  # ocean invasion
                + 0.135 / exp(actual_diftime / 5500.0)  # sea floor CaCO3 neutralization 
                + 0.035 / exp(actual_diftime / 8200.0)  # terrestrial floor CaCO3 neutralization
                + 0.08 / exp(actual_diftime / 200000.0) # silicate weathering
            )   # Archer 1997
        end

        now_dt["co2dot_"*hm] = (co2ref - now_dt["co2_"*hm] + par_dt["ktco2"] * (now_dt["T_"*hm] - par_dt["T_ref_"*hm])) / par_dt["tau_co2"]
    end
    return now_dt
end


