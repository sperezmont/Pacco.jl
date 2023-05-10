# =============================
#     Program: geometry.jl
#     Aim: This program contains functions to calculate ice geometry variables
# =============================

########################
# Diagnostic variables 
########################
"""
    calc_Z!(u, p)
calculates ice surface elevation
"""
function calc_Z!(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        u[14] = 0.0
    elseif p.active_iso
        u[14] = max(u[5] + u[7], # H + B
            0.0 + (1 - (p.rhoi / p.rhow) * u[5]))     # 0 + (1- rhoi/rhow) * H (Pattyn 2017, Robinson 2020)
    else
        u[14] = u[5] + p.Beq  # H + Beq
    end
    return nothing
end

"""
    calc_E!(u, p)
calculates ice-sheet extension
"""
function calc_E!(u::Vector, p::Params)
    if p.active_climate
        u[15] = p.Eref * (u[1] - u[13]) / p.Ate # Eref * (T - Tref) / Ate
    else
        u[15] = p.Eref * (u[12] - u[13]) / p.Ate # Eref * (Tsl - Tref) / Ate
    end
    return nothing
end

"""
    calc_V!(u, p) 
calculates ice-sheet volume from prognostic mean H and diagnosed E
"""
function calc_V!(u::Vector, p::Params)
    u[16] = u[15] * u[5] * p.rhoi / p.rhow / p.Aoc # E * H * rhoi / rhow / Aoc
    return nothing
end