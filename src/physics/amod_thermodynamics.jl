# =============================
#     Program: amod_thermodynamics.jl
#     Aim: This program contains functions to calculate thermodynamics
# =============================
"""
    calc_P(now, par)
calculates surface pressure

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_P(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        now["P_"*hm] = P_sl * exp((-now["Z_"*hm] * g) / (Rd * now["T_surf_"*hm]))   # http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-1/
    end
    return now
end

"""
    calc_T_surf(now, par)
calculates air temperature at ice sheet surface level

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_T_surf(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        if par["active_climate"]
            refT = copy(now["T_"*hm])
        else
            refT = copy(now["T_sl_"*hm])
        end

        if par["tsurf_case"] == "linear"
            now["T_surf_"*hm] = refT - grad * now["Z_"*hm]
        else
            error("ERROR, T_surf option not recognized")
        end
    end
    return now
end

"""
    calc_Acc(now, par)
calculates accumulation rate

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_Acc(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        if par["ac_case"] == "ins"
            pr = par["pr_ref"] + par["A_pr"] * now["ins_norm_"*hm]
            # Calculate the fraction of snow
            if now["T_surf_"*hm] <= (par["T_snow"]) # if below t_snow, full snowfall
                snf = pr
            elseif now["T_surf_"*hm] >= (par["T_rain"]) # if above t_rain, full rain
                snf = 0.0
            else # smooth transition
                f_snow = (now["T_surf_"*hm] - par["T_rain"]) / (par["T_snow"] - par["T_rain"])  # assume linear transition
                snf = f_snow * pr
            end
            now["Acc_"*hm] = max(snf, 0.0)
        elseif par["ac_case"] == "linear"
            temp, now["T_ref_"*hm] = calc_temp_and_tempref(now, par, hm)
            now["Acc_"*hm] = par["Acc_ref_"*hm] + par["ka"] * (temp - now["T_ref_"*hm])
            now["Acc_"*hm] = max(now["Acc_"*hm], 0.0)
        end
    end
    return now
end

"""
    calc_M(now, par)
calculates surface melting rate

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_M(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        temp, now["T_ref_"*hm] = calc_temp_and_tempref(now, par, hm)
        if par["sm_case"] == "PDD"    # positive degree day method, as in Robinson et al. 2010

            if now["T_surf_"*hm] >= (par["T_threshold"])
                now["M_"*hm] = par["lambda"] * (temp - par["T_threshold"])
            else
                now["M_"*hm] = 0.0
            end

        elseif par["sm_case"] == "ITM"

            if now["SMB_"*hm] <= 0
                # I have to test this without >0 condition -- spm 2022.12.19
                albedo_to_use = now["albedo_"*hm]
            else
                # I have to test this without >0 condition -- spm 2022.12.19
                albedo_to_use = par["albedo_newice"]
            end

            now["M_"*hm] = (par["km"]
                              + par["ki"] * max((1 - albedo_to_use) * now["ins_anom_"*hm], 0.0)
                              + par["lambda"] * max(temp - par["T_threshold"], 0.0))

        else
            error("ERROR, surface melt option not recognized")
        end
    end
    return now
end

"""
    calc_SMB(now, par)
calculates surface mass balance

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_SMB(now::OrderedDict, par::OrderedDict)
    # First, calculates Accumulation
    now = calc_Acc(now, par)

    # Second, calculate Melting
    now = calc_M(now, par)

    # Third, return SMB
    for hm in par["hemisphere"]
        now["SMB_"*hm] = now["Acc_"*hm] - now["M_"*hm]
    end
    return now
end

"""
    calc_TMB(now, par)
calculates total mass balance

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_TMB(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        now["TMB_"*hm] = now["SMB_"*hm] + 0.0   # for now, TMB = SMB -- 2022.11.17 spm
    end
    return now
end

"""
    calc_Qdif(now, par)
calculates diffusive heat

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_Qdif(now::OrderedDict, par::OrderedDict)
    kt_ann = par["kt"] * sec_year
    #qgeo_ann = par["Q_geo"] * sec_year * 1e-3

    for hm in par["hemisphere"]
        if now["H_"*hm] < 10.0  # -- check if there is no ice
            now["Q_dif_"*hm] = 0.0
        else
            # -- update ice-air diffusion -- CHECK units!! spm 2022.12.07
            now["Q_difup_"*hm] = -2 * ((now["T_ice_"*hm] - now["T_surf_"*hm]) / (now["H_"*hm]^2)) *
                                   (kt_ann / (par["c"] * rhoi))
            # -- update ice-mantle diffusion
            now["Q_difdown_"*hm] = -2 * ((now["T_ice_"*hm] - par["T_mantle_"*hm]) / (par["H_mantle_"*hm]^2)) *
                                     (kt_ann / (par["c"] * rhoi))
            # -- update diffusion from geothermal flux

            # -- total
            now["Q_dif_"*hm] = now["Q_difup_"*hm] + now["Q_difdown_"*hm]
        end
    end
    return now
end

"""
    calc_Qdrag(now, par)
calculates drag heat

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_Qdrag(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        if now["H_"*hm] < 10.0  # -- check if there is no ice
            now["Q_drag_"*hm] = 0.0
        else
            now["Q_drag_"*hm] = now["fstream_"*hm] * now["tau_b_"*hm] * now["U_b_"*hm] / (par["c"] * rhoi) #/ par["L"] # -- spm 2022.11.24
        end
    end
    return now
end

#############################
# Time derivatives
#############################
"""
    calc_T_icedot(now, par)
calculates ice temperature derivative

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_T_icedot(now, par)
    for hm in par["hemisphere"]
        now["T_icedot_"*hm] = now["Q_dif_"*hm] + now["Q_drag_"*hm]
    end
    return now
end