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
function calcdot_ice_temperature(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if u[H_idx] <= p.ice_is_big_thr  # -- check if there is not enough ice
        return 0.0
    else
        return p.sec_year * (u[Pe_idx]^(1/2)) / (p.rhoice * p.cice * u[H_idx]) * (u[hcond_idx] + u[hvadv_idx] + u[hdrag_idx] + u[hgeo_idx]) 
    end
end

########################
# Diagnostic variables 
########################
"""
    calc_conductive_heat!(u, p)
calculates diffusion heat between the ice sheet and the atmosphere
"""
function calc_conductive_heat!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.active_climate
        temp2use = u[T_idx]  # T
    else
        temp2use = u[Tsurf_idx]    # Tsurf
    end

    if u[H_idx] <= p.ice_is_big_thr  # -- check if there is not enough ice
        u[hcond_idx] = 0.0
    elseif temp2use >= u[Tice_idx]  # stop extracting heat
        u[hcond_idx] = 0.0
    else
        u[hcond_idx] = min(p.kthr * (temp2use - u[Tice_idx]) / (u[H_idx]*(1 - (u[Pe_idx]^(-1/2)))), 0.0)
    end
    return nothing
end

"""
    calc_vertical_heat_advection!(u, p)
    calculates vertical heat advection
"""
function calc_vertical_heat_advection!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.active_climate
        temp2use = u[T_idx]  # T
    else
        temp2use = u[Tsurf_idx]    # Tsurf
    end

    if u[H_idx] < p.ice_is_big_thr  # -- check if there is not enough ice
        u[hvadv_idx] = 0.0
    elseif temp2use >= u[Tice_idx]  # stop advecting heat
        u[hvadv_idx] = 0.0
    else   
        u[hvadv_idx] = min(p.rhoice * p.cice * u[s_idx] * (temp2use - u[Tice_idx]) / (u[Pe_idx]^(1/2) - 1) / p.sec_year, 0.0)
    end
    return nothing
end

"""
    calc_dragging_heat!(u, p)
calculates drag heating
"""
function calc_dragging_heat!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if u[H_idx] < p.ice_is_big_thr  # -- check if there is not enough ice
        u[hdrag_idx] = 0.0
    else
        u[hdrag_idx] = u[taub_idx] * u[vb_idx] / p.sec_year
    end
    return nothing
end

"""
    calc_geothermal_heat!(u, p)
    calculates the geothermal heat contributed to the ice sheet by the geothermal flux
"""
function calc_geothermal_heat!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if u[H_idx] < p.ice_is_big_thr  # -- check if there is not enough ice
        u[hgeo_idx] = 0.0
    else
        u[hgeo_idx] = p.hgeo     
    end
    return nothing
end

"""
    calc_peclet_number!(u, p)
    calculates the Peclet number Pe
"""
function calc_peclet_number!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    kt_ann = p.kthr * p.sec_year # J yr⁻¹ m⁻¹ K⁻¹
    offset = 1e-1

    if u[H_idx] == p.ice_exists_thr  # -- no ice
        u[Pe_idx] = 0.0
    elseif u[H_idx] < p.ice_is_big_thr  # -- check if there is not enough ice
        u[Pe_idx] = 1.0 + offset
    else
        u[Pe_idx] = max(u[H_idx] / kt_ann * p.rhoice * p.cice * u[s_idx], 1/0.9^2) # 1 - Pe^-1/2 ≥ 0.1
    end
    return nothing
end

"""
    calc_temperate_layer!(u, p)
    calculates the temperate layer thickness H_b
"""
function calc_temperate_layer!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.temperate_layer_case == "fixed"
        u[Hb_idx] = p.basal_scale
    elseif p.temperate_layer_case == "verb2018"
        u[Hb_idx] = max((1 / sqrt(p.Pe)) * u[H_idx], 0.0)
    end
    return nothing
end