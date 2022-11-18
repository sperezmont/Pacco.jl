# =============================
#     Program: amod_thermodynamics.jl
#     Aim: This program contains functions to calculate thermodynamics
# =============================

@doc """
    calc_T_surf: calculates surface temperature
"""
function calc_T_surf(now_t, par_t)
    if par_t["tsurf_case"] == "linear"
        return now_t["T_sl"] - grad * now_t["Z"]
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
        e_s = 611.3 * exp(Lv / Rv * (1 / degK - 1 / now_t["T_surf"]))        # Clausius-Clapeyron differential equation direct approximation http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    elseif par_t["cc_case"] == "AERKi"
        t = now_t["T_surf"] - degK    # conversion to ºC
        e_s = 6.1121 * exp((22.587 * t) / (t + 273.86))     # Alduchov and Eskridge (1996)  
    else
        error("ERROR, Clausius-Clapeyron option not recognized")
    end

    # Compute specific humidity
    eps = Rd / Rv
    q_s = (eps * e_s) / (now_t["P"] - e_s * (1 - eps))
    q = q_s * par_t["RH"]

    # Then, calculate precipitation
    pr = (1 + par_t["k_pr"] * now_t["Z"] / par_t["L"]) * q / par_t["tau_w"]  # [kg water yr⁻¹?] REMBO by Robinson et al. (2010)

    # Now, change to meters of ice equivalent per year
    surface = pi * par_t["L"]^2 # assume radial ice sheet ?? -- spm
    pr = pr * rhoi / (rhow)# * surface)

    # Calculate the fraction of snow
    if now_t["T_surf"] <= (par_t["T_snow"]) # if below t_snow, full snowfall
        snf = pr
    elseif now_t["T_surf"] >= (par_t["T_rain"]) # if above t_rain, full rain
        snf = 0.0
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
            return par_t["lambda"] * (now_t["T_surf"] - par_t["melt_offset"])
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
    snf_now = calc_snowfall(now_t, par_t)
    now_t["Acc"] = snf_now + 0.0   # accumulation is assumed to be full snowfall

    # Second, calculate Melting
    now_t["M"] = calc_surfmelt(now_t, par_t)

    # Third, return SMB
    now_t["SMB"] = now_t["Acc"] - now_t["M"]
    return now_t
end

@doc """
    calc_Qdif: calculates diffusive heat
"""
function calc_Qdif(now_t, par_t, ann_kt)
    # Uupdate air diffusion
    now_t["Q_difup"] = -2 * (now_t["T"] - now_t["T_surf"]) / (now_t["H"]^2) * ann_kt / (par_t["c"] * rhoi)

    # Update diffusion from the Mantle
    now_t["Q_difdown"] = -2 * (now_t["T"] - par_t["T_mantle"]) / (par_t["H_mantle"]^2) * ann_kt / (par_t["c"] * rhoi)

    # Update diffusion from geothermal flux

    # Compute the total
    now_t["Q_dif"] = now_t["Q_difup"] + now_t["Q_difdown"]
    return now_t
end

@doc """
    calc_Qdrag: calculates drag heat
"""
function calc_Qdrag(now_t, par_t)
    now_t["Q_drag"] = now_t["fstream"] * now_t["tau_b"] * now_t["U_b"] / (par_t["c"] * rhoi)
    return now_t
end