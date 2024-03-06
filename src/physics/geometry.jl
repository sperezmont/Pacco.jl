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
    if u[5] == p.ice_exists_thr  # no ice
        u[13] = 0.0
    elseif p.active_iso
        u[13] = u[5] + u[6] + u[7]
    else
        u[13] = u[5] + u[6] + p.Beq  # H + Beq
    end
    return nothing
end

"""
    calc_icesheet_area!(u, p)
calculates ice-sheet area (Surf)
"""
function calc_icesheet_area!(u::Vector, p::Params)
    if p.active_climate
        u[14] = 8 * pi * (p.L / 1e3)^2 * ((u[1] - u[12]) / p.Ath - u[23] / p.hrz_vel_scale)
    else
        u[14] = 8 * pi * (p.L / 1e3)^2 * ((u[11] - u[12]) / p.Ath - u[23] / p.hrz_vel_scale)
    end
    return nothing
end

"""
    calc_icesheet_volume!(u, p) 
calculates ice-sheet volume (Vol) from prognostic mean H and diagnosed Surf
"""
function calc_icesheet_volume!(u::Vector, p::Params)
    u[15] = u[14] * u[5] * p.rhoice / p.rhowater / p.Surfoc # S * H * rhoice / rhowater / Surfoc
    return nothing
end