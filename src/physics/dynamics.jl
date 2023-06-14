# =============================
#     Program: dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================

########################
# Prognostic variables 
########################
"""
    calc_Hdot(u, p)
calculates ice thickness derivative through dHdt = (A - M) - U * H / L
"""
function calc_Hdot(u::Vector, p::Params)
    if p.active_ice
        return (u[19] - u[20]) - u[26] * u[5] / p.L # Hdot = (A - M) - U * H / L
    else
        return u[19] - u[20]    # Hdot = A - M
    end
end

"""
    calc_Hseddot(u, p)
calculates sediment layer thickness derivative through dHseddt = -f1 * U + f2 * M
"""
function calc_Hseddot(u::Vector, p::Params)
    return -p.f1 * u[26] + p.f2 * u[20] # Hseddot = -f1 * U + f2 * M
end

"""
    calc_Bdot(u, p)
calculates bed elevation derivative through dBdt = ((Beq - H * rhoi/rhom) - B) / τB
"""
function calc_Bdot(u::Vector, p::Params)
    if p.active_iso
        return ((p.Beq - u[5] * p.rhoi / p.rhom) - u[7]) / p.tau_bed
    else
        return 0.0
    end
end

"""
    calc_fstreamdot(u, p)
calculates stream fraction derivative through dfstreamdt = (fstream_ref - fstream) / tau_kin
"""
function calc_fstreamdot(u::Vector, p::Params)
    return (u[25] - u[9]) / p.tau_kin
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
        u[21] = p.rhoi * p.g * u[5] * u[14] / p.L  # rhoi * g * H * Z / L
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
        u[22] = u[21]  # rhoi * g * H * Z / L
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_Ud!(u, p)
calculates deformational velocity
"""
function calc_Ud!(u::Vector, p::Params)
    if p.dyn_case == "SIA"
        u[23] = (2.0 * p.A_flow * u[5] * (u[21]^p.glen_n)) / (p.glen_n + 2)
    else
        error("Ice-sheet dynamics case not recognized")
    end
    return nothing
end

"""
    calc_Ub!(u, p)
calculates basal velocity
"""
function calc_Ub!(u::Vector, p::Params)
    if p.basal_case == "weertmanq"
        beta = max(0.0, min(1.0, u[6]))       # sediments Csprime weight
        Csprime = beta * p.Cs                  # sliding par calculation based on amount of sediments
        u[24] = Csprime * u[22]^2.0            # Pollard and DeConto (2012) take m = 2 so tau_b^(m-1) * tau_b = tau_b^2  
    else
        error("Ice-sheet basal velocity case not recognized")
    end
    return nothing
end

"""
    calc_fstream_ref!(u, p)
calculates reference value for streaming fraction
"""
function calc_fstream_ref!(u::Vector, p::Params)
    # Streaming inland propagation
    propagation_coef = max(0.0, min(1.0, (u[8] + p.Tsb) / (p.Tsb)))
    u[25] = (p.fstream_max - p.fstream_min) * propagation_coef + p.fstream_min  # the reference value of the streaming fraction is linear with alpha (thus temperature) -- jas
    return nothing
end