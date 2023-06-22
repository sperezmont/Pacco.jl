# =============================
#     Program: geometry.jl
#     Aim: This program contains functions to calculate ice geometry variables
# =============================

########################
# Diagnostic variables 
########################
"""
    calc_z!(u, p)
calculates ice surface elevation
"""
function calc_z!(u::Vector, p::Params)
    if u[5] == 0.0  # no ice
        u[13] = 0.0
    elseif p.active_iso
        u[13] = max(u[5] + u[7], # H + zb
            0.0 + (1 - (p.rhoi / p.rhow) * u[5]))     # 0 + (1- rhoi/rhow) * H (Pattyn 2017, Robinson 2020)
    else
        u[13] = u[5] + p.zbeq  # H + Beq
    end
    return nothing
end

"""
    calc_A!(u, p)
calculates ice-sheet area
"""
function calc_A!(u::Vector, p::Params)
    if p.active_climate
        u[14] = (p.L / 1e3)^2 * (u[1] - u[12]) / p.Ata # L² * (T - Tref) / Ata
    else
        u[14] = (p.L / 1e3)^2 * (u[11] - u[12]) / p.Ata # L² * (Tsl - Tref) / Ata
    end
    return nothing
end

"""
    calc_V!(u, p) 
calculates ice-sheet volume from prognostic mean H and diagnosed A
"""
function calc_V!(u::Vector, p::Params)
    u[15] = u[14] * u[5] * p.rhoi / p.rhow / p.Aoc # A * H * rhoi / rhow / Aoc
    return nothing
end