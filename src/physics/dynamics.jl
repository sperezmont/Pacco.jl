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
    dH/dt = (ṡ - ȧ) - v * H / L
"""
function calcdot_icethickness(u::Vector, p::Params)
    if p.active_ice
        return (u[18] - u[19]) - u[25] * u[5] / p.L # Hdot = (s - a) - v * H / L
    else
        return u[18] - u[19]    # Hdot = s - a
    end
end

"""
    calcdot_sediment_thickness(u, p)
calculates sediment layer thickness derivative through 
    dHsed/dt = -fv * v + fa * a
"""
function calcdot_sediment_thickness(u::Vector, p::Params)
    return -p.fv * u[25] + p.fa * u[19] # Hseddot = -fv * v + fa * a
end

"""
    calcdot_bedrock_elevation(u, p)
calculates bed elevation derivative through
    dB/dt = ((Beq - H * ρi/ρm) - B) / τB
"""
function calcdot_bedrock_elevation(u::Vector, p::Params)
    if p.active_iso
        return ((p.Beq - u[5] * p.rhoice / p.rhobed) - u[7]) / p.taubedrock
    else
        return 0.0
    end
end

"""
    calcdot_streaming_fraction(u, p)
calculates stream fraction derivative through 
    dfstr/dt = (fstr_ref - fstr) / taukin
"""
function calcdot_streaming_fraction(u::Vector, p::Params)
    if p.streaming_case == "fixed"
        return (u[24] - u[9]) / p.taukin
    elseif p.streaming_case == "dynamic"
        return (u[24] - u[9]) / (p.L/u[25])     # Based on Johannesson et al. (1989, Journal of Glaciology)
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
function calc_driving_stress!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[20] = p.rhoice * p.g * u[5] * u[13] / p.L  # rhoice * g * H * Z / L
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_basal_stress!(u, p)
calculates basal stress
"""
function calc_basal_stress!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[21] = u[20]  # rhoice * g * H * Z / L
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_deformational_velocity!(u, p)
calculates deformational velocity
"""
function calc_deformational_velocity!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[22] = (2.0 * p.Aflow * u[5] * (u[20]^p.glen_n)) / (p.glen_n + 2)
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_basal_velocity!(u, p)
calculates basal velocity
"""
function calc_basal_velocity!(u::Vector, p::Params)
    if p.basal_case == "weertmanq"
        u[23] = u[9] * max(0.0, min(1.0, (u[6] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min))) * p.Cs * u[21]^2.0#max(0.0, min(1.0, u[6])) * p.Cs * u[21]^2.0            # Pollard and DeConto (2012): vb = Cs' ⋅ τb²  
    
    elseif p.basal_case == "mb"    
        if (u[18] - u[19]) < 0   # if m = s - a <= 0, use (t)
            u[23] = u[9] * max(0.0, min(1.0, (u[6] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min))) * p.Cs * u[21]^2.0
        else
            u[23] = p.fstrmin * max(0.0, min(1.0, (u[6] - p.Hsed_min) / (p.Hsed_max - p.Hsed_min))) * p.Cs * u[21]^2.0
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
function calc_reference_streaming_fraction!(u::Vector, p::Params)
    # Streaming inland propagation 
    if p.ref_streaming_case == "theo"   # to check if streaming facilitates deglaciations
        if u[5] < 10 # 
            propagation_coef = 0.0
        else    
            propagation_coef = 1.0 #.* u[5] ./ 1000
        end
    elseif p.ref_streaming_case == "thermo" # apply PACCO thermodynamics to compute fstr
        if (u[8] < p.Tstr) # the base is frozen
            propagation_coef = 0.0
        else    # the base becomes temperate and the ice streams grow
            propagation_coef = (u[8] - p.Tstr) / (p.Tmp - p.Tstr)
        end
    else
        error("Reference streaming fraction case not recognized")
    end
    u[24] = min((p.fstrmax - p.fstrmin) * propagation_coef + p.fstrmin, p.fstrmax)  # fstrref goes from fstrmin to fstrmax 
    return nothing
end