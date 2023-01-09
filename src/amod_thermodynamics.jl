# =============================
#     Program: amod_thermodynamics.jl
#     Aim: This program contains functions to calculate thermodynamics
# =============================
@doc """
    calc_P: calculates surface pressure
"""
function calc_P(now_t, par_t)
    for hm in par_t["hemisphere"]
        now_t["P_"*hm] = P_sl * exp((-now_t["Z_"*hm] * g) / (Rd * now_t["T_surf_"*hm]))   # http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-1/
    end
    return now_t
end

@doc """
    calc_T_surf: calculates surface temperature
"""
function calc_T_surf(now_t, par_t)
    for hm in par_t["hemisphere"]
        if par_t["active_climate"]
            refT = copy(now_t["T_"*hm])
        else
            refT = copy(now_t["T_sl_"*hm])
        end

        if par_t["tsurf_case"] == "linear"
            now_t["T_surf_"*hm] = refT - grad * now_t["Z_"*hm]
        else
            error("ERROR, T_surf option not recognized")
        end
    end
    return now_t
end

@doc """
    calc_Acc: calculates accumulation rate
"""
function calc_Acc(now_t, par_t)
    for hm in par_t["hemisphere"]
        if par_t["ac_case"] == "ins"
            pr = par_t["pr_ref"] + par_t["A_pr"] * now_t["ins_norm_"*hm]
            # Calculate the fraction of snow
            if now_t["T_surf_"*hm] <= (par_t["T_snow"]) # if below t_snow, full snowfall
                snf = pr
            elseif now_t["T_surf_"*hm] >= (par_t["T_rain"]) # if above t_rain, full rain
                snf = 0.0
            else # smooth transition
                f_snow = (now_t["T_surf_"*hm] - par_t["T_rain"]) / (par_t["T_snow"] - par_t["T_rain"])  # assume linear transition
                snf = f_snow * pr
            end
            now_t["Acc_"*hm] = snf
        elseif par_t["ac_case"] == "linear"
            now_t["Acc_"*hm] = par_t["Acc_ref_"*hm] + par_t["ka"] * (now_t["T_"*hm] - par_t["T_ref_"*hm])
            now_t["Acc_"*hm] = max(now_t["Acc_"*hm], 0.0)
        end
    end
    return now_t
end

@doc """
    calc_M: calculates surface melting rate
"""
function calc_M(now_t, par_t)
    for hm in par_t["hemisphere"]
        if par_t["sm_case"] == "PDD"    # positive degree day method, as in Robinson et al. 2010
            if now_t["T_surf_"*hm] >= (par_t["melt_offset"])
                now_t["M_"*hm] = par_t["lambda"] * (now_t["T_surf_"*hm] - par_t["melt_offset"])
            else
                now_t["M_"*hm] = 0.0
            end
        elseif par_t["sm_case"] == "ITM"
            if now_t["SMB_"*hm] <= 0
                # I have to test this without >0 condition -- spm 2022.12.19
                now_t["M_"*hm] = (par_t["km"]
                                 + par_t["ki"] * max((1 - now_t["albedo_"*hm]) * now_t["ins_anom_"*hm], 0.0)
                                 + par_t["lambda"] * max(now_t["T_"*hm] - par_t["T_ref_"*hm], 0.0))
            else
                # I have to test this without >0 condition -- spm 2022.12.19
                now_t["M_"*hm] = (par_t["km"]
                                 + par_t["ki"] * max((1 - par_t["albedo_newice"]) * now_t["ins_anom_"*hm], 0.0)
                                 + par_t["lambda"] * max(now_t["T_"*hm] - par_t["T_ref_"*hm], 0.0))
            end
        else
            error("ERROR, surface melt option not recognized")
        end
    end
    return now_t
end

@doc """
    calc_SMB: calculates surface mass balance
"""
function calc_SMB(now_t, par_t)
    # First, calculates Accumulation
    now_t = calc_Acc(now_t, par_t)

    # Second, calculate Melting
    now_t = calc_M(now_t, par_t)

    # Third, return SMB
    for hm in par_t["hemisphere"]
        now_t["SMB_"*hm] = now_t["Acc_"*hm] - now_t["M_"*hm]
    end
    return now_t
end

@doc """
    calc_TMB: calculates total mass balance
"""
function calc_TMB(now_t, par_t)
    for hm in par_t["hemisphere"]
        now_t["TMB_"*hm] = now_t["SMB_"*hm] + 0.0   # for now, TMB = SMB -- 2022.11.17 spm
    end
    return now_t
end

@doc """
    calc_Qdif: calculates diffusive heat
"""
function calc_Qdif(now_t, par_t)
    kt_ann = par_t["kt"] * sec_year
    #qgeo_ann = par["Q_geo"] * sec_year * 1e-3

    for hm in par_t["hemisphere"]
        if now_t["H_"*hm] < 10.0  # -- check if there is no ice
            now_t["Q_dif_"*hm] = 0.0
        else
            # -- update ice-air diffusion -- CHECK units!! spm 2022.12.07
            now_t["Q_difup_"*hm] = -2 * ((now_t["T_ice_"*hm] - now_t["T_surf_"*hm]) / (now_t["H_"*hm]^2)) *
                                   (kt_ann / (par_t["c"] * rhoi))
            # -- update ice-mantle diffusion
            now_t["Q_difdown_"*hm] = -2 * ((now_t["T_ice_"*hm] - par_t["T_mantle_"*hm]) / (par_t["H_mantle_"*hm]^2)) *
                                     (kt_ann / (par_t["c"] * rhoi))
            # -- update diffusion from geothermal flux

            # -- total
            now_t["Q_dif_"*hm] = now_t["Q_difup_"*hm] + now_t["Q_difdown_"*hm]
        end
    end
    return now_t
end

@doc """
    calc_Qdrag: calculates drag heat
"""
function calc_Qdrag(now_t, par_t)
    for hm in par_t["hemisphere"]
        if now_t["H_"*hm] < 10.0  # -- check if there is no ice
            now_t["Q_drag_"*hm] = 0.0
        else
            now_t["Q_drag_"*hm] = now_t["fstream_"*hm] * now_t["tau_b_"*hm] * now_t["U_b_"*hm] / (par_t["c"] * rhoi) #/ par_t["L"] # -- spm 2022.11.24
        end
    end
    return now_t
end