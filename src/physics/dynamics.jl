# =============================
#     Program: dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================

########################
# Prognostic variables 
########################
"""
    calcdot_icethickness(u, p)
calculates ice thickness derivative through 
    dH/dt = (ṡ - ȧ) - v * H / Locn
"""
function calcdot_icethickness(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.active_ice
        if p.ice_discharge_case == "fixed"
            return (u[s_idx] - u[a_idx]) - u[v_idx] * u[H_idx] / p.L # Hdot = (s - a) - v * H / L
        elseif p.ice_discharge_case == "dynamic"
            return (u[s_idx] - u[a_idx]) - u[v_idx] * u[H_idx] / u[L_idx] # Hdot = (s - a) - v * H / L(H)
        end
    else
        return u[s_idx] - u[a_idx]    # Hdot = s - a
    end
end

"""
    calcdot_sediment_thickness(u, p)
calculates sediment layer thickness derivative through 
    dHsed/dt = -fv * v + fa * a
"""
function calcdot_sediment_thickness(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    return -p.fv * u[v_idx] + p.fa * u[a_idx] # Hseddot = -fv * v + fa * a
end

"""
    calcdot_bedrock_elevation(u, p)
calculates bed elevation derivative through
    dB/dt = ((Beq - H * ρi/ρm) - B) / τB
"""
function calcdot_bedrock_elevation(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.active_iso
        return ((p.Beq - u[H_idx] * p.rhoice / p.rhobed) - u[B_idx]) / p.taubedrock
    else
        return 0.0
    end
end

"""
    calcdot_streaming_fraction(u, p)
calculates stream fraction derivative through 
    dfstr/dt = (fstr_ref - fstr) / taukin
"""
function calcdot_streaming_fraction(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.streaming_case == "fixed"
        return (u[fstr_ref_idx] - u[fstr_idx]) / p.taukin
    elseif p.streaming_case == "dynamic"
        return (u[fstr_ref_idx] - u[fstr_idx]) / (u[L_idx]/u[v_idx])     # Based on Johannesson et al. (1989, Journal of Glaciology)
    else
        error("Streaming fraction case not recognized")
    end
end

########################
# Diagnostic variables 
########################
"""
    calc_driving_stress!(u, p)
calculates driving stress
"""
function calc_driving_stress!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if u[H_idx] <= p.ice_is_big_thr  # no ice
        u[taud_idx] = 0.0
    else
        if (p.deformationalflow_case == "profile")
            u[taud_idx] = p.rhoice * p.g * u[H_idx] * u[z_idx] / u[L_idx]  # rhoice * g * H * Z / L
        else
            error("Ice-sheet dynamics case not recognized")
        end
    end
    return nothing
end

"""
    calc_basal_stress!(u, p)
calculates basal stress
"""
function calc_basal_stress!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.plugflow_case == "profile"
        u[taub_idx] = u[taud_idx]  # rhoice * g * H * Z / L
    elseif p.plugflow_case == "constant"
        u[taub_idx] = p.rhoice * p.g * u[H_idx] * p.sintheta # rhoice * g * H * sin(θ)
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_deformational_velocity!(u, p)
calculates deformational velocity
"""
function calc_deformational_velocity!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if (p.deformational_case == "glen58")
        u[vd_idx] = (2.0 * p.Aflow * u[H_idx] * (u[taud_idx]^p.glen_n)) / (p.glen_n + 2)
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_beta!(u, p)
calculates basal velocity coefficient
"""
function calc_beta!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    u[Hsed_idx] = max(p.Hsed_min, u[Hsed_idx])
    u[beta_idx] = max(p.beta_min, min(1.0, (u[Hsed_idx] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min)))

    return nothing
end

"""
    calc_basal_velocity!(u, p)
calculates basal velocity
"""
function calc_basal_velocity!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.basal_case == "weertmanq"
        # u[vb_idx] = u[fstr_idx] * max(0.0, min(1.0, (u[Hsed_idx] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min))) * p.Cs * u[taub_idx]^2.0            # Pollard and DeConto (2012): vb = Cs' ⋅ τb²  
        
        u[vb_idx] = u[fstr_idx] * u[beta_idx] * p.Cs * u[taub_idx]^2.0            # Pollard and DeConto (2012): vb = Cs' ⋅ τb²  
    
    elseif p.basal_case == "mb"    
        if (u[s_idx] - u[a_idx]) < 0   # if m = s - a <= 0, use (t)
            u[vb_idx] = u[fstr_idx] * max(0.0, min(1.0, (u[Hsed_idx] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min))) * p.Cs * u[taub_idx]^2.0
        else
            u[vb_idx] = p.fstrmin * max(0.0, min(1.0, (u[Hsed_idx] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min))) * p.Cs * u[taub_idx]^2.0
        end

    else
        error("Ice-sheet basal velocity case not recognized")
    end
    return nothing
end

"""
    calc_reference_streaming_fraction!(u, p)
calculates reference value for streaming fraction
"""
function calc_reference_streaming_fraction!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    # Streaming inland propagation 
    if p.ref_streaming_case == "theo"   # to check if streaming facilitates deglaciations
        if u[H_idx] < 10 # 
            propagation_coef = 0.0
        else    
            propagation_coef = 1.0 #.* u[H_idx] ./ 1000
        end
    elseif p.ref_streaming_case == "thermo" # apply PACCO thermodynamics to compute fstr
        if (u[Tice_idx] < p.Tstr) # the base is frozen
            propagation_coef = 0.0
        else    # the base becomes temperate and the ice streams grow
            propagation_coef = (u[Tice_idx] - p.Tstr) / (p.Tmp - p.Tstr)
        end
    else
        error("Reference streaming fraction case not recognized")
    end
    u[fstr_ref_idx] = min((p.fstrmax - p.fstrmin) * propagation_coef + p.fstrmin, p.fstrmax)  # fstrref goes from fstrmin to fstrmax 
    return nothing
end