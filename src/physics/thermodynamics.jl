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
    return u[26] + u[27]
end

########################
# Diagnostic variables 
########################
"""
    calc_s!(u, p)
calculates snowfall accumulation rate
"""
function calc_s!(u::Vector, p::Params)
    if p.s_case == "ins"
        Inorm = 2.0 * (u[10] - p.Imin) / (p.Imax - p.Imin) - 1.0
        pr = p.pr_ref + p.Apr * Inorm

        # Calculate the fraction of snow
        if u[17] <= p.Tsnow # if below t_snow, full snowfall
            snf = pr
        elseif u[17] >= p.Train # if above t_rain, full rain
            snf = 0.0
        else # smooth transition
            fsnow = (u[17] - p.Train) / (p.Tsnow - p.Train)  # assume linear transition
            snf = fsnow * pr
        end
        u[18] = max(snf, 0.0)
    elseif p.s_case == "linear"
        if p.active_climate
            u[18] = max(p.Aref + p.ks * (u[1] - u[12]), 0.0)   # Aref + ks * (T - Tref)
        else
            u[18] = max(p.Aref + p.ks * (u[17] - u[12]), 0.0)   # Aref + ks * (Tsurf - Tref)
        end
    end
    return nothing
end

"""
    calc_a!(u, p)
calculates ablation rate
"""
function calc_a!(u::Vector, p::Params)
    if p.a_case == "PDD"    # positive degree day method, as in Robinson et al. 2010
        if u[17] >= (p.Tthreshold)
            if p.active_climate
                u[19] = p.lambda * (u[1] - p.Tthreshold)   # lambda(T - Tthreshold)
            else
                u[19] = p.lambda * (u[17] - p.Tthreshold)   # lambda(Tsurf - Tthreshold)
            end
        else
            u[19] = 0.0
        end
    elseif p.a_case == "ITM"
        if (u[18] - u[19]) <= 0   # if m = s - a <= 0, use alpha(t)
            u[19] = p.km + p.ki * max((1 - u[4]) * (u[10] - p.Iref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)    # alpha
        else
            u[19] = p.km + p.ki * max((1 - p.alphanewice) * (u[10] - p.Iref), 0.0) + p.lambda * max(u[1] - p.Tthreshold, 0.0)   # alphaₙ
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
        u[26] = 0.0
    else
        # -- update ice-air diffusion -- CHECK units!! spm 2022.12.07
        Q_difup = -2 * ((u[8] - u[17]) / (u[5]^2)) * (kt_ann / (p.c * p.rhoi))
        # -- update ice-mantle diffusion
        Q_difdown = -2 * ((u[8] - p.Tmantle) / (p.Hmantle^2)) * (kt_ann / (p.c * p.rhom))
        # -- update diffusion from geothermal flux

        # -- total
        u[26] = Q_difup + Q_difdown
    end
    return nothing
end

"""
    calc_Qdrag!(u, p)
calculates drag heating
"""
function calc_Qdrag!(u::Vector, p::Params)
    if u[5] < 10.0  # -- check if there is no ice
        u[27] = 0.0
    else
        u[27] = u[9] * u[21] * u[23] / (p.c * p.rhoi) #/ p.L # -- spm 2022.11.24
    end
    return nothing
end