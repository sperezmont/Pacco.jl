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
function calc_icesheet_elevation!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if u[H_idx] == p.ice_exists_thr  # no ice
        u[z_idx] = 0.0
    elseif p.active_iso
        u[z_idx] = u[H_idx] + u[Hsed_idx] + u[B_idx]
    else
        u[z_idx] = u[H_idx] + u[Hsed_idx] + p.Beq  # H + Beq
    end
    return nothing
end

"""
    calc_icesheet_area!(u, p)
calculates ice-sheet area (Surf)
"""
function calc_icesheet_area!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.active_climate
        u[Surf_idx] = pi * (u[L_idx] / 1e3)^2 * ((u[T_idx] - u[Tref_idx]) / p.Ath - u[vb_idx] / 50.0)
    else
        u[Surf_idx] = pi * (u[L_idx] / 1e3)^2 * ((u[Tsl_idx] - u[Tref_idx]) / p.Ath - u[vb_idx] / p.hrz_vel_scale)
    end
    return nothing
end

"""
    calc_icesheet_volume!(u, p) 
calculates ice-sheet volume (Vol) from prognostic mean H and diagnosed Surf
"""
function calc_icesheet_volume!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    u[Vol_idx] = u[Surf_idx] * u[H_idx] * p.rhoice / p.rhowater / p.Surfoc # S * H * rhoice / rhowater / Surfoc
    return nothing
end

"""
    calc_aspect_ratio!(u, p)
calculates ice-sheet aspect ratio/horizontal scale (L)
"""
function calc_aspect_ratio!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if u[H_idx] <= p.ice_is_big_thr  # no ice
        u[L_idx] = 0.0
    else
        if p.ice_hrz_scale_case == "fixed"
            u[L_idx] = p.L
        elseif p.ice_hrz_scale_case == "dynamic"    # following Verbitsky et al. (2018)
            #u[L_idx] = min(max(p.hrz_scale_coeff * u[H_idx]^p.hrz_scale_exp, p.hrz_scale_lb), p.hrz_scale_ub)
            u[L_idx] = min(max(p.hrz_scale_coeff * u[z_idx]^p.hrz_scale_exp, p.hrz_scale_lb), p.hrz_scale_ub)
        end
    end
    return nothing
end