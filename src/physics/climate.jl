# =============================
#     Program: climate.jl
#     Aim: This program contains functions to calculate climate variables
# =============================

########################
# Prognostic variables 
########################
"""
    calcdot_regional_temperature(u, p)
calculates air thermal relaxation through 
    dT/dt = ((Tref + RI + RCO2 - cZ * Z) - T) / tauT
"""
function calcdot_regional_temperature(u::Vector, p::Params)
    if p.insol_case == "ISI"
        RI = p.cISI * (u[10] - p.insol_ref)
    elseif p.insol_case == "caloric"
        RI = p.cCAL * (u[10] - p.insol_ref)
    else
        RI = p.cI * (u[10] - p.insol_ref)
    end

    RCO2 = p.cC * calc_carbon_dioxide_rad(u[2])
    return (u[12] + RI + RCO2 - p.cZ * u[13] - u[1]) / p.tauT
end

"""
    calcdot_carbon_dioxide(u, p, t)
calculates carbon dioxide (C) derivative through 
    dC/dt = ((Cref + kTC * (T - Tref)) - C) / tauC
"""
function calcdot_carbon_dioxide(u::Vector, p::Params, t::Real)
    if p.carbon_case in ["dynamic", "trended"]
        if t < p.time_anth # unperturbed climate
            reference_carbon = p.Cref + p.kTC * (u[1] - u[12]) 
        else    # perturbed climate
            actual_diftime = t - p.time_anth
            Cref = p.Cref + p.C_anth * (
                0.75 / exp(actual_diftime / 365.0)  # ocean invasion
                + 0.135 / exp(actual_diftime / 5500.0)  # sea floor CaCO3 neutralization 
                + 0.035 / exp(actual_diftime / 8200.0)  # terrestrial floor CaCO3 neutralization
                + 0.08 / exp(actual_diftime / 200000.0) # silicate weathering
            )   # Archer 1997
            reference_carbon = Cref + p.kTC * (u[1] - u[12])
        end

        if p.carbon_case == "dynamic"
            return (reference_carbon - u[2]) / p.tauC
        elseif p.carbon_case == "trended"
            return p.kC + (reference_carbon - u[2]) / p.tauC
        end

    elseif p.carbon_case == "constant"
        return 0.0
    end
end

"""
    calcdot_albedo(u, p)
calculates albedo derivative through 
    dα/dt = (αref - α) / τα
"""
function calcdot_albedo(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        return 0.0
    elseif u[3] < 10.0 # first ice
        return 0.0
    else
        return (u[16] - u[4]) / p.tau_albedo
    end
end

"""
    calcdot_iceage(u)
calculates ice age
    dA/dt = 1.0
"""
function calcdot_iceage()
    return 1.0
end

########################
# Diagnostic variables 
########################
"""
    calc_sealevel_temperature!(u, p, t)
calculates sea level temperature
    Tsl = Tref₀ + Aₜ ⋅ Iₙₒᵣₘ
"""
function calc_sealevel_temperature!(u::Vector, p::Params, t::Real)
    # First, normalize insolation
    normalized_insol = 2.0 * (u[10] - p.insol_min) / (p.insol_max - p.insol_min) - 1.0   # between 1 and -1, norm = 2

    # Second, compute anthropogenic forcing if time >= time_anth
    if (t >= p.time_anth)
        u[11] = p.Tref0 + p.At * normalized_insol + p.AT_anth / exp((t - p.time_anth) / p.tau_anth)
    else
        u[11] = p.Tref0 + p.At * normalized_insol
    end
    return nothing
end

"""
    calc_reference_temperature!(u, p, t)
calculates climatic reference temperature
"""
function calc_reference_temperature!(u::Vector, p::Params, t::Real)
    # Anthropogenic forcing?
    if t < p.time_anth # unperturbed climate
        u[12] = p.Tref0
    else    # perturbed climate
        u[12] = p.Tref0 + p.kCT * u[2] / p.Cref
    end
    return nothing
end

"""
    calc_reference_albedo!(u, p)
calculates reference value for albedo
    αref = αnewice - kα ⋅ A   
"""
function calc_reference_albedo!(u::Vector, p::Params)
    if p.active_aging == true
        if u[5] == 0.0  # no ice
            u[16] = p.albedo_land  # albedoref is albedo_land
        else
            n = 1 # exponent number in albedo-age parameterisation
            #albedoslope = (1+ (10-1)*u[6]) * p.albedo_slope
            u[16] = max(p.albedo_newice - p.k_albedo * u[3]^n, p.albedo_oldice) 
        end
    elseif p.active_aging == false
        if u[5] == 0.0  # no ice
            u[16] = p.albedo_land  # albedoref is albedo_land
        else
            u[16] = p.albedo_newice
        end
    end
    return nothing
end

"""
    calc_surface_temperature!(u, p)
calculates air temperature at surface level
"""
function calc_surface_temperature!(u::Vector, p::Params)
    surface = u[5] + u[6] + u[7]    # always, even if there is no ice

    if p.active_climate
        u[17] = u[1] - p.Γ * surface    # reference is regional T
    else
        u[17] = u[11] - p.Γ * surface    # reference is Tsl
    end
    return nothing
end

"""
    calc_s!(u, p)
calculates snowfall accumulation rate
"""
function calc_snowfall!(u::Vector, p::Params)
    if p.snowfall_case == "ins"
        normalized_inso = 2.0 * (u[10] - p.insol_min) / (p.insol_max - p.insol_min) - 1.0
        pr = p.pr_ref + p.Apr * normalized_inso

        # Calculate the fraction of snow
        if u[17] <= p.Tsnow # if below t_snow, full snowfall
            snf = pr
        elseif u[17] >= p.Train # if above t_rain, full rain
            snf = 0.0
        else # smooth transition
            fsnow = (u[17] - p.Train) / (p.Tsnow - p.Train)  # assume linear transition
            snf = fsnow * pr
        end
        u[18] = max(snf, 0.0)
    elseif p.snowfall_case == "linear"
        if p.active_climate
            u[18] = max(p.sref + p.ks * (u[1] - u[12]), 0.0)   # sref + ks * (T - Tref)
        else
            u[18] = max(p.sref + p.ks * (u[17] - u[12]), 0.0)   # sref + ks * (Tsurf - Tref)
        end
    end
    return nothing
end

"""
    calc_ablation!(u, p)
calculates ablation rate
"""
function calc_ablation!(u::Vector, p::Params)
    if p.ablation_case == "PDD"    # positive degree day method, as in Robinson et al. 2010
        if u[17] >= (p.Tthreshold)
            if p.active_climate
                u[19] = p.lambda * (u[1] - p.Tthreshold)   # lambda(T - Tthreshold)
            else
                u[19] = p.lambda * (u[17] - p.Tthreshold)   # lambda(Tsurf - Tthreshold)
            end
        else
            u[19] = 0.0
        end
    elseif p.ablation_case == "PDD-LIN" # PDD linear with sea level temperature
        if u[11] >= (p.Tthreshold)
            u[19] = p.lambda * (u[11] - p.Tthreshold)   # lambda(Tsl - Tthreshold)
        else
            u[19] = 0.0
        end
    elseif p.ablation_case == "ITM"
        if p.active_snow_on_ice # let snow cover hide the old ice (if smb is ≥ 0)
            if (u[18] - u[19]) <= 0   # if m = s - a <= 0, use albedo(t)
                u[19] = p.km + p.kI * max((1 - u[4]) * (u[10] - p.insol_ref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)    # albedo
            else
                u[19] = p.km + p.kI * max((1 - p.albedo_newice) * (u[10] - p.insol_ref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)   # albedoₙ
            end
        else    # use old ice albedo (time evolving albedo)
            u[19] = p.km + p.kI * max((1 - u[4]) * (u[10] - p.insol_ref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)    # albedo
        end

    else
        error("ERROR, surface melting option not recognized")
    end
    return nothing
end

########################
# Other variables 
########################
"""
    calc_carbon_dioxide_rad(pCO2)
calculates the radiative forcing due to carbon dioxide (in W/m²), Myhre et al. (1998)
"""
function calc_carbon_dioxide_rad(pCO2::Real)
    return 5.35 * NaNMath.log(pCO2 / 280.0)
end