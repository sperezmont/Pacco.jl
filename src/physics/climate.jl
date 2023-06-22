# =============================
#     Program: climate.jl
#     Aim: This program contains functions to calculate climate variables
# =============================

########################
# Prognostic variables 
########################
"""
    calc_Tdot(u, p)
calculates air thermal relaxation through dTdt = ((Tref + RI + RCO2 - cz * Z) - T) / tauT
"""
function calc_Tdot(u::Vector, p::Params)
    if p.I_case == "ISI"
        RI = p.cisi * (u[10] - p.Iref)
    elseif p.I_case == "caloric"
        RI = p.ccal * (u[10] - p.Iref)
    else
        RI = p.ci * (u[10] - p.Iref)
    end

    RCO2 = p.cc * calc_rad_pCO2(u[2], p)
    return (u[12] + RI + RCO2 - p.cz * u[13] - u[1]) / p.tauT
end

"""
    calc_pCO2dot(u, p, t)
calculates pCO2 derivative through dpCO2dt = ((pCO2ref + ktpCO2 * (T - Tref)) - pCO2) / taupCO2
"""
function calc_pCO2dot(u::Vector, p::Params, t::Real)
    # -- anthropogenic forcing?
    if t < p.time_anth # unperturbed climate
        return (p.pCO2ref + p.ktpCO2 * (u[1] - u[12]) - u[2]) / p.taupCO2
    else    # perturbed climate
        actual_diftime = t - p.time_anth
        pCO2ref = p.pCO2ref + p.pCO2anth * (
            0.75 / exp(actual_diftime / 365.0)  # ocean invasion
            + 0.135 / exp(actual_diftime / 5500.0)  # sea floor CaCO3 neutralization 
            + 0.035 / exp(actual_diftime / 8200.0)  # terrestrial floor CaCO3 neutralization
            + 0.08 / exp(actual_diftime / 200000.0) # silicate weathering
        )   # Archer 1997
        new_pCO2ref = pCO2ref + p.ktpCO2 * (u[1] - u[12])
        return (new_pCO2ref - u[2]) / p.taupCO2
    end

end

"""
    calc_alphadot(u, p)
calculates alpha derivative through dalphadt = (alpharef - alpha) / taualpha
"""
function calc_alphadot(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        return 0.0
    elseif u[3] < 10.0 # first ice
        return 0.0
    else
        return (u[16] - u[4]) / p.taualpha
    end
end

"""
    calc_iceagedot(u)
calculates ice age
"""
function calc_iceagedot()
    return 1.0
end

########################
# Diagnostic variables 
########################
"""
    calc_Tsl!(u, p, t)
calculates sea level temperature
"""
function calc_Tsl!(u::Vector, p::Params, t::Real)
    # First, normalize insolation
    Inorm = 2.0 * (u[10] - p.Imin) / (p.Imax - p.Imin) - 1.0   # between 1 and -1, norm = 2

    # Second, compute anthropogenic forcing if time >= time_anth
    if (t >= p.time_anth)
        u[11] = p.Tref0 + p.At * Inorm + p.Atanth / exp((t - p.time_anth) / p.tauanth)
    else
        u[11] = p.Tref0 + p.At * Inorm
    end
    return nothing
end

"""
    calc_Tref!(u, p, t)
calculates climatic reference temperature
"""
function calc_Tref!(u::Vector, p::Params, t::Real)
    # Anthropogenic forcing?
    if t < p.time_anth # unperturbed climate
        u[12] = p.Tref0
    else    # perturbed climate
        u[12] = p.Tref0 + p.cpCO2 * u[2] / p.pCO2ref
    end
    return nothing
end

"""
    calc_alpha_ref!(u, p)
calculates reference value for albedo
"""
function calc_alpha_ref!(u::Vector, p::Params)
    if p.active_aging == true
        if u[5] == 0.0  # no ice
            u[16] = p.alphaland  # alpharef is alphaland
        else
            n = 1 # exponent number in albedo-age parameterisation
            #alphaslope = (1+ (10-1)*u[6]) * p.alpha_slope
            u[16] = max(p.alphanewice - p.kalpha * u[3]^n, p.alphaoldice) # alphaₙ - kalpha * t
        end
    elseif p.active_aging == false
        if u[5] == 0.0  # no ice
            u[16] = p.alphaland  # alpharef is alphaland
        else
            u[16] = p.alphanewice
        end
    end
    return nothing
end

"""
    calc_Tsurf!(u, p)
calculates air temperature at ice sheet surface level
"""
function calc_Tsurf!(u::Vector, p::Params)
    if p.active_climate
        u[17] = u[1] - p.Γ * u[13]    # reference is regional T
    else
        u[17] = u[11] - p.Γ * u[13]    # reference is Tsl
    end
    return nothing
end

########################
# Other variables 
########################
"""
    calc_rad_pCO2(CO2)
calculates the radiative forcing due to pCO2 (in W/m²), Myhre et al. (1998)
"""
function calc_rad_pCO2(CO2::Real, p)
    return 5.35 * NaNMath.log(CO2 / p.pCO2ref)
end