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
function calcdot_ice_temperature(u::Vector, p::Params)
    if u[5] <= p.ice_is_big_thr  # -- check if there is not enough ice
        return 0.0
    else
        return p.sec_year / (p.rhoice * p.cice * p.basal_scale) * (u[26] + u[27] + u[28] + u[29])  # syr/(ρcH_b) * (hcond + hvadv + hdrag + hgeo)
    end
end

########################
# Diagnostic variables 
########################
"""
    calc_conductive_heat!(u, p)
calculates diffusion heat between the ice sheet and the atmosphere
"""
function calc_conductive_heat!(u::Vector, p::Params)
    kt_ann = p.kthr * p.sec_year # J yr⁻¹ m⁻¹ K⁻¹

    if p.active_climate
        temp2use = u[1]  # T
    else
        temp2use = u[17]    # Tsurf
    end

    if u[5] <= p.ice_is_big_thr  # -- check if there is not enough ice
        u[26] = 0.0
    elseif temp2use - u[8] >= 0 # the environment is warmer than the bedrock
        u[26] = 0.0
    else
        u[26] = p.kthr * (temp2use - u[8]) / u[5]
        # if p.diffusion_case == "2pts"   # assume only two points in the column (simplest case)
        #     #u[26] = ((temp2use - u[8]) / (u[5]^2)) * (kt_ann / (p.cice * p.rhoice)) 
        #      # W/m²
        # elseif p.diffusion_case == "3pts"   # assume three points in the column (based on Moreno-Parada et al., 2024, TC)
        #     third_point_temp = u[8]
        #     #u[26] = ((temp2use -2 * third_point_temp + u[8]) / (u[5]^2)) * (kt_ann / (p.cice * p.rhoice)) 
        # else
        #     error("Diffusion case case not recognized")
        # end
    end
    return nothing
end

"""
    calc_vertical_heat_advection!(u, p)
    calculates vertical heat advection
"""
function calc_vertical_heat_advection!(u::Vector, p::Params)
    if p.active_climate
        temp2use = u[1]  # T
    else
        temp2use = u[17]    # Tsurf
    end

    if u[5] <= p.ice_is_big_thr  # -- check if there is not enough ice
        u[27] = 0.0
    elseif u[8] - temp2use <= 0 # the environment is warmer than the bedrock
        u[27] = 0.0
    else    
        u[27] = 0.0#p.rhoice * p.cice * u[18] / p.sec_year * p.basal_scale / u[5] * (temp2use - u[8]) 
    end
    return nothing
end

"""
    calc_dragging_heat!(u, p)
calculates drag heating
"""
function calc_dragging_heat!(u::Vector, p::Params)
    if u[5] <= p.ice_is_big_thr  # -- check if there is not enough ice
        u[28] = 0.0
    else
        #u[28] = - u[21] * u[23] / (p.kthr * p.sec_year) * p.vrt_vel_scale  # - vb * τb / kthr * ω
        u[28] = u[21] * u[23] / p.sec_year 
    end
    return nothing
end

"""
    calc_geothermal_heat!(u, p)
    calculates the geothermal heat contributed to the ice sheet by the geothermal flux
"""
function calc_geothermal_heat!(u::Vector, p::Params)
    #u[28] = - p.qgeo / p.kthr * p.vrt_vel_scale    # - qgeo / kthr * ω
    u[29] = p.hgeo 
    return nothing
end
