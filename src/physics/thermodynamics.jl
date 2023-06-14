# =============================
#     Program: thermodynamics.jl
#     Aim: This program contains functions to calculate ice-sheet thermodynamics
# =============================

########################
# Prognostic variables 
########################
"""
    calc_Ticedot(u)
calculates ice temperature derivative through dTicedt = Qdif + Qdrag
"""
function calc_Ticedot(u::Vector)
    return u[27] + u[28]
end

########################
# Diagnostic variables 
########################
"""
    calc_A!(u, p)
calculates accumulation rate
"""
function calc_A!(u::Vector, p::Params)
    if p.A_case == "ins"
        Inorm = 2.0 * (u[10] - p.I_min) / (p.I_max - p.I_min) - 1.0
        pr = p.pr_ref + p.A_pr * Inorm

        # Calculate the fraction of snow
        if u[18] <= p.Tsnow # if below t_snow, full snowfall
            snf = pr
        elseif u[18] >= p.Train # if above t_rain, full rain
            snf = 0.0
        else # smooth transition
            fsnow = (u[18] - p.Train) / (p.Tsnow - p.Train)  # assume linear transition
            snf = fsnow * pr
        end
        u[19] = max(snf, 0.0)
    elseif p.A_case == "linear"
        if p.active_climate
            u[19] = max(p.Aref + p.ka * (u[1] - u[13]), 0.0)   # Aref + ka * (T - Tref)
        else
            u[19] = max(p.Aref + p.ka * (u[18] - u[13]), 0.0)   # Aref + ka * (Tsurf - Tref)
        end
    end
    return nothing
end

"""
    calc_M!(u, p)
calculates melting rate
"""
function calc_M!(u::Vector, p::Params)
    if p.M_case == "PDD"    # positive degree day method, as in Robinson et al. 2010
        if u[18] >= (p.Tthreshold)
            if p.active_climate
                u[20] = p.lambda * (u[1] - p.Tthreshold)   # λ(T - Tthreshold)
            else
                u[20] = p.lambda * (u[18] - p.Tthreshold)   # λ(Tsurf - Tthreshold)
            end
        else
            u[20] = 0.0
        end
    elseif p.M_case == "ITM"
        if (u[19] - u[20]) <= 0   # if MB = A - M <= 0, use α(t)
            u[20] = p.km + p.ki * max((1 - u[4]) * (u[10] - p.I_ref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)    # α
        else
            u[20] = p.km + p.ki * max((1 - p.alpha_newice) * (u[10] - p.I_ref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)   # αₙ
        end

    else
        error("ERROR, surface melting option not recognized")
    end
    return nothing
end

"""
    calc_Qdif!(u, p)
calculates diffusion heat between the ice sheet and the atmosphere
"""
function calc_Qdif!(u::Vector, p::Params)
    kt_ann = p.kt * p.sec_year
    #qgeo_ann = p.Qgeo * p.sec_year * 1e-3

    if u[5] < 10.0  # -- check if there is no ice
        u[27] = 0.0
    else
        # -- update ice-air diffusion -- CHECK units!! spm 2022.12.07
        Q_difup = -2 * ((u[8] - u[18]) / (u[5]^2)) * (kt_ann / (p.c * p.rhoi))
        # -- update ice-mantle diffusion
        Q_difdown = -2 * ((u[8] - p.Tmantle) / (p.Hmantle^2)) * (kt_ann / (p.c * p.rhoi))
        # -- update diffusion from geothermal flux

        # -- total
        u[27] = Q_difup + Q_difdown
    end
    return nothing
end

"""
    calc_Qdrag!(u, p)
calculates drag heating
"""
function calc_Qdrag!(u::Vector, p::Params)
    if u[5] < 10.0  # -- check if there is no ice
        u[28] = 0.0
    else
        u[28] = u[9] * u[22] * u[24] / (p.c * p.rhoi) #/ p.L # -- spm 2022.11.24
    end
    return nothing
end