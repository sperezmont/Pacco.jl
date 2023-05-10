# =============================
#     Program: climate.jl
#     Aim: This program contains functions to calculate climate variables
# =============================

########################
# Prognostic variables 
########################
"""
    calc_Tdot(u, p)
calculates air thermal relaxation through dTdt = ((Tref + R - cz * Z) - T) / τᵢ
"""
function calc_Tdot(u::Vector, p::Params)
    return ((u[13] + u[11] - p.cz * u[14]) - u[1]) / p.tau_T
end

"""
    calc_co2dot(u, p, t)
calculates co2 derivative through dco2dt = ((co2ref + ktco2 * (T - Tref)) - co2) / τco2
"""
function calc_co2dot(u::Vector, p::Params, t::Real)
    # -- anthropogenic forcing?
    if t < p.time_anth # unperturbed climate
        return (p.co2_ref + p.ktco2 * (u[1] - u[13]) - u[2]) / p.tau_co2
    else    # perturbed climate
        actual_diftime = t - p.time_anth
        co2ref = p.co2_ref + p.co2_anth * (
            0.75 / exp(actual_diftime / 365.0)  # ocean invasion
            + 0.135 / exp(actual_diftime / 5500.0)  # sea floor CaCO3 neutralization 
            + 0.035 / exp(actual_diftime / 8200.0)  # terrestrial floor CaCO3 neutralization
            + 0.08 / exp(actual_diftime / 200000.0) # silicate weathering
        )   # Archer 1997
        new_co2ref = co2ref + p.ktco2 * (u[1] - u[13])
        return (new_co2ref - u[2]) / p.tau_co2
    end

end

"""
    calc_alphadot(u, p)
calculates α derivative through dαdt = (αref - α) / τα
"""
function calc_alphadot(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        return 0.0
    elseif u[3] < 10.0 # first ice
        return 0.0
    else
        return (u[17] - u[4]) / p.tau_alpha
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
    calc_R!(u, p)
calculates radiative forcing
"""
function calc_R!(u::Vector, p::Params)
    u[11] = p.ci * (u[10] - p.I_ref) + p.cc * calc_rad_co2(u[2])
    return nothing
end

"""
    calc_Tsl!(u, p, t)
calculates sea level temperature
"""
function calc_Tsl!(u::Vector, p::Params, t::Real)
    # First, normalize insolation
    I_norm = 2.0 * (u[10] - p.I_min) / (p.I_max - p.I_min) - 1.0   # between 1 and -1, norm = 2

    # Second, compute anthropogenic forcing if time >= time_anth
    if (t >= p.time_anth)
        u[12] = p.Tref0 + p.At * I_norm + p.At_anth / exp((t - p.time_anth) / p.tau_anth)
    else
        u[12] = p.Tref0 + p.At * I_norm
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
        u[13] = p.Tref0
    else    # perturbed climate
        u[13] = p.Tref0 + p.cco2 * u[2] / p.co2_ref
    end
    return nothing
end


"""
    calc_alpha_ref!(u, p)
calculates reference value for albedo
"""
function calc_alpha_ref!(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        u[17] = p.alpha_land  # alpha_ref is alpha_land
    else
        n = 1 # exponent number in albedo-age parameterisation
        u[17] = max(p.alpha_newice - p.alpha_slope * u[3]^n, p.alpha_oldice) # αₙ - αₛ * t
    end
    return nothing
end

"""
    calc_Tsurf!(u, p)
calculates air temperature at ice sheet surface level
"""
function calc_Tsurf!(u::Vector, p::Params)
    if p.active_climate
        u[18] = u[1] - p.lapse_rate * u[14]    # reference is regional T
    else
        u[18] = u[12] - p.lapse_rate * u[14]    # reference is Tsl
    end
    return nothing
end

########################
# Other variables 
########################
"""
    calc_rad_co2(CO2)
calculates the radiative forcing due to co2 (in W/m²)
"""
function calc_rad_co2(CO2::Real)
    return 5.35 * log(CO2 / 280.0)
end