# =============================
#     Program: geometry.jl
#     Aim: This program contains functions to calculate ice geometry variables
# =============================

########################
# Diagnostic variables 
########################
"""
    calc_icesheet_elevation!(u, p)
calculates ice surface elevation
"""
function calc_icesheet_elevation!(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        u[13] = 0.0
    elseif p.active_iso
        u[13] = max(u[5] + u[7], # H + Beq
            0.0 + (1 - (p.rhoi / p.rhow) * u[5]))     # 0 + (1- rhoi/rhow) * H (Pattyn 2017, Robinson 2020)
    else
        u[13] = u[5] + p.Beq  # H + Beq
    end
    return nothing
end

"""
    calc_icesheet_area!(u, p)
calculates ice-sheet area (Surf)
"""
function calc_icesheet_area!(u::Vector, p::Params)
    if p.active_climate
        u[14] = (p.L / 1e3)^2 * (u[1] - u[12]) / p.Ath # L² * (T - Tref) / Ath
    else
        u[14] = (p.L / 1e3)^2 * (u[11] - u[12]) / p.Ath # L² * (Tsl - Tref) / Ath
    end
    return nothing
end

"""
    calc_icesheet_volume!(u, p) 
calculates ice-sheet volume (Vol) from prognostic mean H and diagnosed Surf
"""
function calc_icesheet_volume!(u::Vector, p::Params)
    u[15] = u[14] * u[5] * p.rhoi / p.rhow / p.Surfoc # S * H * rhoi / rhow / Aoc
    return nothing
end