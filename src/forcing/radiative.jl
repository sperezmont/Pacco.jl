# =============================
#     Program: radiative.jl
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
        now_r["T_sl_"*hm] = par_r["T_ref0_"*hm] + par_r["A_t"] * now_r["ins_norm_"*hm] + T_anth
    end
    return now_r
end

"""
    calc_temp_and_tempref(now, par, hm)
selects temperature to use and computes reference temperature desired

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters
* `hm` Hemisphere in which we calculate ("_n" or "_s")

## Return
`temp` and `now["T_ref_"*hm]` to use
"""
function calc_temp_and_tempref(now::OrderedDict, par::OrderedDict, hm)
    if par["active_climate"]
        if par["height_temp"] == "useH"
            temp = now["T_surf_"*hm]
        elseif par["height_temp"] == "useZ"
            temp = now["T_"*hm]
        else
            printstyled("dev par must be removed!", color=:red)
        end

    else
        temp = now["T_surf_"*hm]  # if only dynamics, we take into account the lapse rate here
    end

    # Anthropogenic forcing?
    if now["time"] < par["time_anth"] # unperturbed climate
        if par["height_temp"] == "useH"
            now["T_ref_"*hm] = par["T_ref0_"*hm] - grad * now["Z_"*hm]
        elseif par["height_temp"] == "useZ"
            now["T_ref_"*hm] = par["T_ref0_"*hm]
        end
    else    # perturbed climate
        if par["height_temp"] == "useH"
            now["T_ref_"*hm] = par["T_ref0_"*hm] - grad * now["Z_"*hm]
            now["T_ref_"*hm] += par["cco2"] * now["co2_"*hm] / par["co2_ref"]
        elseif par["height_temp"] == "useZ"
            now["T_ref_"*hm] = par["T_ref0_"*hm] + par["cco2"] * now["co2_"*hm] / par["co2_ref"]
        end
    end

    return temp, now["T_ref_"*hm]
end

#############################
# Time derivatives
#############################
@doc """
    calc_Tdot:
        calculates regional temperature derivative
"""
function calc_Tdot(now, par)
    for hm in par["hemisphere"]

        # Anthropogenic forcing?
        if now["time"] < par["time_anth"] # unperturbed climate
            now["T_ref_"*hm] = par["T_ref0_"*hm]
        else    # perturbed climate
            now["T_ref_"*hm] = par["T_ref0_"*hm] + par["cco2"] * now["co2_"*hm] / par["co2_ref"]
        end

        if par["height_temp"] == "useH"
            HTF = now["H_"*hm]
        elseif par["height_temp"] == "useZ"
            HTF = now["Z_"*hm]
        else
            printstyled("dev par must be removed!", color=:red)
        end

        now["Tdot_"*hm] = (now["T_ref_"*hm] - now["T_"*hm]
                              +
                              now["T_rf_"*hm]
                              -
                              par["csz"] * HTF) / par["tau_rf_"*hm]
    end
    return now
end

@doc """
    calc_albedodot:
        calculates albedo derivative
"""
function calc_albedodot(now, par)
    for hm in par["hemisphere"]
        now["albedodot_"*hm] = (now["albedo_ref_"*hm] - now["albedo_"*hm]) / par["tau_albedo"]

        if now["H_"*hm] == 0.0
            now["albedo_"*hm] = par["albedo_land"]
            now["albedodot_"*hm] = 0.0
        elseif now["ice_time_"*hm] < 10.0 # First ice
            now["albedo_"*hm] = par["albedo_newice"]
            now["albedodot_"*hm] = 0.0
        end
    end
    return now
end

@doc """
    calc_co2dot:
        calculates co2 derivative
"""
function calc_co2dot(now, par)
    for hm in par["hemisphere"]

        # -- anthropogenic forcing?
        if now["time"] < par["time_anth"] # unperturbed climate
            co2ref = par["co2_ref"]
        else    # perturbed climate
            actual_diftime = now["time"] - par["time_anth"]
            co2ref = par["co2_ref"] + par["co2_anth"] * (
                0.75 / exp(actual_diftime / 365.0)  # ocean invasion
                + 0.135 / exp(actual_diftime / 5500.0)  # sea floor CaCO3 neutralization 
                + 0.035 / exp(actual_diftime / 8200.0)  # terrestrial floor CaCO3 neutralization
                + 0.08 / exp(actual_diftime / 200000.0) # silicate weathering
            )   # Archer 1997
        end

        temp, now["T_ref_"*hm] = calc_temp_and_tempref(now, par, hm)
        now["co2dot_"*hm] = (co2ref - now["co2_"*hm] + par["ktco2"] * (temp - now["T_ref_"*hm])) / par["tau_co2"]
    end
    return now
end


