# =============================
#     Program: dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================

########################
# Prognostic variables 
########################
"""
    calc_Hdot(u, p)
calculates ice thickness derivative through dHdt = (s - a) - v * H / L
"""
function calc_Hdot(u::Vector, p::Params)
    if p.active_ice
        return (u[18] - u[19]) - u[25] * u[5] / p.L # Hdot = (s - a) - v * H / L
    else
        return u[18] - u[19]    # Hdot = s - a
    end
end

"""
    calc_Hseddot(u, p)
calculates sediment layer thickness derivative through dHseddt = -f1 * v + f2 * a
"""
function calc_Hseddot(u::Vector, p::Params)
    return -p.f1 * u[25] + p.f2 * u[19] # Hseddot = -f1 * v + f2 * a
end

"""
    calc_zbdot(u, p)
calculates bed elevation derivative through dzbdt = ((zbeq - H * rhoi/rhom) - zb) / taubed
"""
function calc_zbdot(u::Vector, p::Params)
    if p.active_iso
        return ((p.zbeq - u[5] * p.rhoi / p.rhom) - u[7]) / p.taubed
    else
        return 0.0
    end
end

"""
    calc_fstrdot(u, p)
calculates stream fraction derivative through dfstrdt = (fstr_ref - fstr) / taukin
"""
function calc_fstrdot(u::Vector, p::Params)
    return (u[24] - u[9]) / p.taukin
end

########################
# Diagnostic variables 
########################
"""
    calc_taud!(u, p)
calculates driving stress
"""
function calc_taud!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[20] = p.rhoi * p.g * u[5] * u[13] / p.L  # rhoi * g * H * z / L
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_taub!(u, p)
calculates basal stress
"""
function calc_taub!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[21] = u[20]  # rhoi * g * H * z / L
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_vd!(u, p)
calculates deformational velocity
"""
function calc_vd!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[22] = (2.0 * p.Aflow * u[5] * (u[20]^p.glen_n)) / (p.glen_n + 2)
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_vb!(u, p)
calculates basal velocity
"""
function calc_vb!(u::Vector, p::Params)
    if p.basal_case == "weertmanq"
        u[23] = max(0.0, min(1.0, u[6])) * p.Cs * u[21]^2.0            # Pollard and DeConto (2012) take m = 2 so tau_b^(m-1) * tau_b = tau_b^2  
    else
        error("Ice-sheet basal velocity case not recognized")
    end
    return nothing
end

"""
    calc_fstr_ref!(u, p)
calculates reference value for streaming fraction
"""
function calc_fstr_ref!(u::Vector, p::Params)
    # Streaming inland propagation
    propagation_coef = max(0.0, min(1.0, (u[8] + p.Tsb) / (p.Tsb)))
    u[24] = (p.fstrmax - p.fstrmin) * propagation_coef + p.fstrmin  # the reference value of the streaming fraction is linear with p (thus temperature) -- jas
    return nothing
end