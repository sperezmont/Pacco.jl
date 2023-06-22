# =============================
#     Program: thermodynamics.jl
#     Aim: This program contains functions to calculate ice-sheet thermodynamics
# =============================

########################
# Prognostic variables 
########################
"""
    calcdot_ice_temperature(u)
calculates ice temperature derivative through
    dTice/dt = Qdif + Qdrag
"""
function calcdot_ice_temperature(u::Vector)
    return u[26] + u[27]
end

########################
# Diagnostic variables 
########################
"""
    calc_diffusional_heat!(u, p)
calculates diffusion heat between the ice sheet and the atmosphere
"""
function calc_diffusional_heat!(u::Vector, p::Params)
    kt_ann = p.kth * p.sec_year
    #qgeo_ann = p.Qgeo * p.sec_year * 1e-3

    if u[5] < 10.0  # -- check if there is no ice
        u[26] = 0.0
    else
        # -- update ice-air diffusion -- CHECK units!! spm 2022.12.07
        Q_difup = -2 * ((u[8] - u[17]) / (u[5]^2)) * (kt_ann / (p.cice * p.rhoi))
        # -- update ice-mantle diffusion
        Q_difdown = -2 * ((u[8] - p.Tmantle) / (p.Hmantle^2)) * (kt_ann / (p.cice * p.rhom))
        # -- update diffusion from geothermal flux

        # -- total
        u[26] = Q_difup + Q_difdown
    end
    return nothing
end

"""
    calc_dragging_heat!(u, p)
calculates drag heating
"""
function calc_dragging_heat!(u::Vector, p::Params)
    if u[5] < 10.0  # -- check if there is no ice
        u[27] = 0.0
    else
        u[27] = u[9] * u[21] * u[23] / (p.cice * p.rhoi) #/ p.L # -- spm 2022.11.24
    end
    return nothing
end