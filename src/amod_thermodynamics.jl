# =============================
#     Program: amod_thermodynamics.jl
#     Aim: This program contains functions to calculate thermodynamics
# =============================

@doc """
    calc_T_surf: calculates surface temperature
"""
function calc_T_surf(now_t, par_t)
    if par_t["tsurf_case"] == "linear"
        return now_t["T_sl"] - grad * now_t["S"]
    else
        error("ERROR, T_surf option not recognized")
    end
end

@doc """
    calc_snowfall: calculates snowfall rate
"""
function calc_snowfall(now_t, par_t)
    # First, calculate the saturation vapor pressure e_s
    if par_t["cc_case"] == "cc"
        T_0 = t_0 + degK
        e_s = 0.6113 * exp(Lv / Rv * (1 / T_0 - 1 / now_t["T_surf"]))          # Clausius-Clapeyron differential equation direct approximation
    elseif par_t["cc_case"] == "AERKi"
        t = now_t["T_surf"] - degK    # conversion to ºC
        e_s = 6.1121 * exp((22.587 * t) / (t + 273.86))     # Alduchov and Eskridge (1996)  
    else
        error("ERROR, Clausius-Clapeyron option not recognized")
    end

    # now_t, compute specific humidity
    eps = Rd / Rv
    q_s = (eps * e_s) / (now_t["P"] - e_s * (1 - eps))
    q = q_s * par_t["RH"]

    # Then, calculate precipitation
    pr = (1 + k_pr * now_t["S"] / par_t["L"]) * q / par_t["tau_w"]  # REMBO by Robinson et al. (2010)

    # Calculate the fraction of snow
    if now_t["T_surf"] <= (par_t["T_snow"]) # if below t_snow, full snowfall
        snf = pr
    elseif now_t["T_surf"] >= (par_t["T_rain"]) # if above t_rain, full rain
        snf = 0
    else # smooth transition
        f_snow = (now_t["T_surf"] - par_t["T_rain"]) / (par_t["T_snow"] - par_t["T_rain"])  # assume linear transition
        snf = f_snow * pr
    end
    return snf
end

@doc """
    calc_surfmelt: calculates surface melting rate
"""
function calc_surfmelt(now_t, par_t)
    if par_t["sm_case"] == "PDD"    # positive degree day method, as in Robinson et al. 2010
        if now_t["T_surf"] >= (par_t["melt_offset"])
            return par_t["lambda"] * (now_t["T_surf"] + par_t["melt_offset"])
        else
            return 0.0
        end
    elseif par_t["sm_case"] == "ITM"
        error("ERROR, surface melt option not implemented yet")
    else
        error("ERROR, surface melt option not recognized")
    end
end

@doc """
    calc_SMB: calculates surface mass balance
"""
function calc_SMB(now_t, par_t)
    # First, calculates Accumulation
    now_t["Acc"] = calc_snowfall(now_t, par_t)  # accumulation is assumed to be full snowfall

    # Second, calculate Melting
    now_t["M"] = calc_surfmelt(now_t, par_t)

    # Third, return SMB
    return now_t["Acc"] + now_t["M"]
end

@doc """
    calc_Qdif: calculates diffusive heat
"""
function calc_Qdif(now_t, par_t, ann_kt)
    # Uupdate air diffusion
    now_t["Q_difup"] = -2 * (now_t["T"] - now_t["T_surf"]) / (now_t["H"]^2) * ann_kt / (par_t["c"] * rho)

    # Update diffusion from the Mantle
    now_t["Q_difdown"] = -2 * (now_t["T"] - par_t["T_mantle"]) / (par_t["H_mantle"]^2) * ann_kt / (par_t["c"] * rho)

    # Update diffusion from geothermal flux

    # Compute the total
    return now_t["Q_difup"] + now_t["Q_difdown"]

end