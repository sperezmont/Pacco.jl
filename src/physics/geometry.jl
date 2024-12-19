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
        #u[Surf_idx] = 8 * pi * (p.L / 1e3)^2 * ((u[T_idx] - u[Tref_idx]) / p.Ath - u[vb_idx] / p.hrz_vel_scale) # old
        #u[Surf_idx] = -pi * (u[L_idx] / 1e3)^2 * (-(u[T_idx] - u[Tref_idx]) / p.Ath + u[vb_idx] / p.hrz_vel_scale) # after reviews
        #u[Surf_idx] = -pi * (u[L_idx] / 1e3)^2 * (1 .+ u[vb_idx] / max(u[v_idx], 1e-3)) 

        # opt1, opt2 = true, false
        # if true
        #     nu = p.hrz_vel_scale
        # else
        #     nu = u[v_idx]
        # end


        # kdef = max(- (u[T_idx] - u[Tref_idx]) / p.Ath, 0.0) # only allow positive values
        # (nu != 0.0) ? (kbas = u[vb_idx] / nu) : (kbas = 0.0)  # if there are basal velocities, the potentially glaciated surface is bigger

        # sdef = pi * kdef * (u[L_idx])^2 / 1e6 # km2 # via an anomaly in the surface

        # if opt1 # too much signal from H
        #     u[Surf_idx] = sdef + kbas * (pi * p.hrz_scale_ub^2 / 1e6 - sdef)
        # elseif opt2 # Vol is huge, postMPT is much bigger than preMPT
        #     n = 2
        #     u[Surf_idx] = (1 + kbas)^n * sdef
        # else
        #     u[Surf_idx] = 0.0
        # end

        u[Surf_idx] = pi * max(- (u[T_idx] - u[Tref_idx]) / p.Ath, 0.0) * (u[L_idx])^2 / 1e6 + u[vb_idx] / p.hrz_vel_scale * (pi * p.hrz_scale_ub^2 / 1e6 - pi * max(- (u[T_idx] - u[Tref_idx]) / p.Ath, 0.0) * (u[L_idx])^2 / 1e6)

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
    u[Vol_idx] = u[Surf_idx] * u[H_idx] * -1 * p.rhoice / p.rhowater / p.Surfoc # S * H * rhoice / rhowater / Surfoc
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