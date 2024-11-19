# =============================
#     Program: climate.jl
#     Aim: This program contains functions to calculate climate variables
# =============================

########################
# Prognostic variables 
########################
"""
    calcdot_regional_temperature(u, p)
calculates air thermal relaxation through 
    dT/dt = ((Tref + RI + RCO2 - cZ * Z) - T) / tauT
"""
function calcdot_regional_temperature(u::Vector{T}, p::Params{T}, t::T) where {T<:AbstractFloat}
    if p.regtemp_case in ["dynamic", "trend", "comp"]
        if p.insol_case == "ISI"    # integrated summer insolation
            RI = p.cISI * (u[I_idx] - p.insol_ref)

        elseif p.insol_case == "caloric"    # seasons
            RI = p.cCAL * (u[I_idx] - p.insol_ref)

        elseif p.insol_case == "input"
            insol_file = split(p.insol_input, "/")[3]
            if (insol_file[1:3] == "ISI") || (insol_file[1:8] == "mean_ISI")
                RI = p.cISI * (u[I_idx] - p.insol_ref)
            elseif (insol_file[1:7] == "caloric") || (insol_file[1:12] == "mean_caloric")
                RI = p.cCAL * (u[I_idx] - p.insol_ref)
            else
                RI = p.cI * (u[I_idx] - p.insol_ref)
            end
        else    # SSI, summer solstice insolation
            RI = p.cI * (u[I_idx] - p.insol_ref)

        end

        RCO2 = p.cC * calc_carbon_dioxide_rad(u[C_idx])

        return (u[Tref_idx] + RI + RCO2 - p.cZ * u[z_idx] - u[T_idx]) / p.tauT

    elseif p.regtemp_case == "constant"
        return 0.0
    elseif p.regtemp_case == "slope"
        return p.kT 
    else
        error("ERROR, regional temperature case not recognized")
    end
end

"""
    calcdot_carbon_dioxide(u, p, t)
calculates carbon dioxide (C) derivative through 
    dC/dt = ((Cref + kTC * (T - Tref)) - C) / tauC
"""
function calcdot_carbon_dioxide(u::Vector{T}, p::Params{T}, t::T) where {T<:AbstractFloat}
    if p.carbon_case in ["dynamic", "trend", "comp"]
        reference_carbon = u[Cref_idx] + p.kTC * (u[T_idx] - u[Tref_idx])

        return (reference_carbon - u[C_idx]) / p.tauC

    elseif p.carbon_case == "constant"
        return 0.0
    elseif p.carbon_case == "slope"
        return p.kC
    else
        error("ERROR, carbon cycle case not recognized")
    end
end



"""
    calcdot_albedo(u, p)
calculates albedo derivative through 
    dα/dt = (αref - α) / τα
"""
function calcdot_albedo(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.albedo_case == "prognostic"
        if u[H_idx] == p.ice_exists_thr  # no ice        
            return 0.0
        elseif u[iceage_idx] < p.ice_is_old_thr # first ice
            return 0.0
        else
            return (u[albedo_ref_idx] - u[albedo_idx]) / p.tau_albedo
        end
    elseif p.albedo_case == "diagnostic"
        return 0.0
    elseif p.albedo_case == "diagnostic_relaxed"
        if u[H_idx] == p.ice_exists_thr  # no ice        
            return 0.0
        elseif u[iceage_idx] < p.ice_is_old_thr # first ice
            return 0.0
        else
            return (u[albedo_ref_idx] - u[albedo_idx]) / 1e3
        end
    end
end

"""
    calcdot_iceage(u)
calculates ice age
    dA/dt = 1.0
"""
function calcdot_iceage()
    return 1.0
end

########################
# Diagnostic variables 
########################
"""
    calc_sealevel_temperature!(u, p, t)
calculates sea level temperature
    Tsl = Tref₀ + Aₜ ⋅ Iₙₒᵣₘ
"""
function calc_sealevel_temperature!(u::Vector{T}, p::Params{T}, t::T) where {T<:AbstractFloat}
    # First, normalize insolation
    normalized_insol = 2.0 * (u[I_idx] - p.insol_min) / (p.insol_max - p.insol_min) - 1.0   # between 1 and -1, norm = 2

    # Second, compute anthropogenic forcing if time >= time_anth
    if (t >= p.time_anth)
        u[Tsl_idx] = p.Tref + p.At * normalized_insol + p.AT_anth / exp((t - p.time_anth) / p.tau_anth)
    else
        u[Tsl_idx] = p.Tref + p.At * normalized_insol
    end
    return nothing
end

"""
    calc_reference_temperature!(u, p, t)
calculates climatic reference temperature
"""
function calc_reference_temperature!(u::Vector{T}, p::Params{T}, t::T) where {T<:AbstractFloat}
    # Anthropogenic forcing?
    if t < p.time_anth # unperturbed climate
        if p.regtemp_case == "trend"
            u[Tref_idx] = p.Tref + p.deltaT + p.kT * (t - p.time_init)
        else
            u[Tref_idx] = p.Tref
        end
    else    # perturbed climate
        u[Tref_idx] = p.Tref + p.kCT * u[C_idx] / p.Cref
    end
    return nothing
end

"""
    calc_reference_carbon!(u, p, t)
calculates climatic reference temperature
"""
function calc_reference_carbon!(u::Vector{T}, p::Params{T}, t::T) where {T<:AbstractFloat}
    if t < p.time_anth # unperturbed climate

        if p.carbon_case == "trend"
            u[Cref_idx] = p.Cref + p.deltaC + p.kC * (t - p.time_init)
        else
            u[Cref_idx] = p.Cref
        end

    else    # perturbed climate
        actual_diftime = t - p.time_anth
        u[Cref_idx] = p.Cref + p.C_anth * (
            0.75 / exp(actual_diftime / 365.0)  # ocean invasion
            + 0.135 / exp(actual_diftime / 5500.0)  # sea floor CaCO3 neutralization 
            + 0.035 / exp(actual_diftime / 8200.0)  # terrestrial floor CaCO3 neutralization
            + 0.08 / exp(actual_diftime / 200000.0) # silicate weathering
        )   # Archer 1997
    end

    return nothing
end

"""
    calc_reference_albedo!(u, p)
calculates reference value for albedo
    αref = αnewice - kα ⋅ A   
"""
function calc_reference_albedo!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.albedo_case == "prognostic"
        if p.active_aging == true
            if u[H_idx] <= p.ice_exists_thr  # no ice
                u[albedo_ref_idx] = p.albedo_land  # albedoref is albedo_land
            else
                n = 1 # exponent number in albedo-age parameterisation
                #albedoslope = (1+ (10-1)*u[Hsed_idx]) * p.albedo_slope
                u[albedo_ref_idx] = max(p.albedo_newice - p.k_albedo * u[iceage_idx]^n, p.albedo_oldice) 
            end
        elseif p.active_aging == false
            if u[H_idx] == p.ice_exists_thr  # no ice
                u[albedo_ref_idx] = p.albedo_land  # albedoref is albedo_land
            else
                u[albedo_ref_idx] = p.albedo_newice
            end
        end
    elseif p.albedo_case in ["diagnostic", "diagnostic_relaxed"]
        if p.active_aging == true
            if u[H_idx] <= p.ice_exists_thr  # no ice
                u[albedo_ref_idx] = p.albedo_land  # albedoref is albedo_land
            else
                u[albedo_ref_idx] = max(p.albedo_newice - (p.albedo_newice - p.albedo_oldice) / p.tau_albedo * u[iceage_idx], p.albedo_oldice) 
            end
        elseif p.active_aging == false
            if u[H_idx] <= p.ice_exists_thr  # no ice
                u[albedo_ref_idx] = p.albedo_land  # albedoref is albedo_land
            else
                u[albedo_ref_idx] = p.albedo_newice
            end
        end
    end

    return nothing
end

"""
    calc_surface_temperature!(u, p)
calculates air temperature at surface level
"""
function calc_surface_temperature!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.active_climate
        u[Tsurf_idx] = u[T_idx] - p.Γ * (u[H_idx] + u[Hsed_idx] + u[B_idx])   # reference is regional T
    else
        u[Tsurf_idx] = u[Tsl_idx] - p.Γ * (u[H_idx] + u[Hsed_idx] + u[B_idx])    # reference is Tsl
    end
    return nothing
end

"""
    calc_s!(u, p)
calculates snowfall accumulation rate
"""
function calc_snowfall!(u::Vector{T}, p::Params{T}, t::T) where {T<:AbstractFloat}
    if p.snowfall_case == "ins"
        normalized_inso = 2.0 * (u[I_idx] - p.insol_min) / (p.insol_max - p.insol_min) - 1.0
        pr = p.pr_ref + p.Apr * normalized_inso

        # Calculate the fraction of snow
        if u[Tsurf_idx] <= p.Tsnow # if below t_snow, full snowfall
            snf = pr
        elseif u[Tsurf_idx] >= p.Train # if above t_rain, full rain
            snf = 0.0
        else # smooth transition
            fsnow = (u[Tsurf_idx] - p.Train) / (p.Tsnow - p.Train)  # assume linear transition
            snf = fsnow * pr
        end
        u[s_idx] = max(snf, 0.0)
    elseif p.snowfall_case == "linear"
        if p.active_climate
            u[s_idx] = max(p.sref + p.ks * (u[T_idx] - u[Tref_idx]), 0.0)   # sref + ks * (T - Tref)
        else
            u[s_idx] = max(p.sref + p.ks * (u[Tsurf_idx] - u[Tref_idx]), 0.0)   # sref + ks * (Tsurf - Tref)
        end
    elseif p.ablation_case == "PDD-LIN"
        u[s_idx] = max(p.sref + p.ks * (u[Tsl_idx] - u[Tref_idx]), 0.0)   # sref + ks * (Tsl - Tref)
    elseif p.snowfall_case == "variable"

        if t < -1250e3
            snow = 1 * p.sref
            ksnow = 1.8 * p.ks
        else
            snow = p.sref
            ksnow = p.ks
        end

        if p.active_climate
            u[s_idx] = max(snow + ksnow * (u[T_idx] - u[Tref_idx]), 0.0)   # sref + ks * (T - Tref)
        else
            u[s_idx] = max(snow + ksnow * (u[Tsurf_idx] - u[Tref_idx]), 0.0)   # sref + ks * (Tsurf - Tref)
        end
    end
    return nothing
end

"""
    calc_ablation!(u, p)
calculates ablation rate
"""
function calc_ablation!(u::Vector{T}, p::Params{T}) where {T<:AbstractFloat}
    if p.ablation_case == "PDD"    # positive degree day method, as in Robinson et al. 2010
        if u[Tsurf_idx] >= (p.Tthreshold)
            if p.active_climate
                u[a_idx] = p.lambda * (u[T_idx] - p.Tthreshold)   # lambda(T - Tthreshold)
            else
                u[a_idx] = p.lambda * (u[Tsurf_idx] - p.Tthreshold)   # lambda(Tsurf - Tthreshold)
            end
        else
            u[a_idx] = 0.0
        end
    elseif p.ablation_case == "PDD-LIN" # PDD linear with sea level temperature
        if u[Tsl_idx] >= (p.Tthreshold)
            u[a_idx] = p.lambda * (u[Tsl_idx] - p.Tthreshold)   # lambda(Tsl - Tthreshold)
        else
            u[a_idx] = 0.0
        end
    elseif p.ablation_case == "ITM"
        if p.active_snow_on_ice # let snow cover hide the old ice (if smb is ≥ 0)
            if (u[s_idx] - u[a_idx]) <= 0   # if m = s - a <= 0, use albedo(t)
                u[albedo_eff_idx] = u[albedo_idx]
                u[a_idx] = p.km + p.kI * max((1 - u[albedo_idx]) * (u[I_idx] - p.insol_ref), 0.0) + p.lambda * max(u[T_idx] - p.Tthreshold, 0.0)    # albedo
            else
                u[albedo_eff_idx] = p.albedo_newice
                u[a_idx] = p.km + p.kI * max((1 - p.albedo_newice) * (u[I_idx] - p.insol_ref), 0.0) + p.lambda * max(u[T_idx] - p.Tthreshold, 0.0)   # albedoₙ
            end
        else    # use old ice albedo (time evolving albedo)
            u[a_idx] = p.km + p.kI * max((1 - u[albedo_idx]) * (u[I_idx] - p.insol_ref), 0.0) + p.lambda * max(u[T_idx] - p.Tthreshold, 0.0)    # albedo
        end

    else
        error("ERROR, surface melting option not recognized")
    end
    return nothing
end

########################
# Other variables 
########################
"""
    calc_carbon_dioxide_rad(pCO2)
calculates the radiative forcing due to carbon dioxide (in W/m²), Myhre et al. (1998)
"""
function calc_carbon_dioxide_rad(pCO2::Real)
    return 5.35 * NaNMath.log(pCO2 / 280.0)
end