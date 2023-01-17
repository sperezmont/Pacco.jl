# =============================
#     Program: amod_derivative.jl
#     Aim: This program contains functions for calculating derivatives
# =============================
@doc """
    calc_Hdot: calculates ice thickness derivative
"""
function calc_Hdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        # -- climatic term
        now_dt["Hdot_"*hm] = copy(now_dt["TMB_"*hm])
        # -- dynamic term
        if par_dt["active_ice"]
            now_dt["Hdot_"*hm] += -now_dt["U_"*hm] * now_dt["H_"*hm] / par_dt["L"]
        end
    end
    return now_dt
end

@doc """
    calc_Hseddot: calculates sediments thickness derivative
"""
function calc_Hseddot(now_dt, par_dt, ctl_dt)
    for hm in par_dt["hemisphere"]
        now_dt["Hseddot_"*hm] = -par_dt["f_1"] * now_dt["U_"*hm] + par_dt["f_2"] * now_dt["M_"*hm] / ctl_dt["dt"]
    end
    return now_dt
end

@doc """
    calc_Bdot: calculates bedrock elevation derivative
"""
function calc_Bdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        if par_dt["active_iso"]
            now_dt["Bdot_"*hm] = -(now_dt["B_"*hm] - par_dt["B_eq_"*hm] + now_dt["H_"*hm] * rhoi / rhom) / par_dt["tau_bed_"*hm] # needs further improvement -- spm 2022.11.17
        else
            now_dt["Bdot_"*hm] = 0.0
        end
    end
    return now_dt
end

@doc """
    calc_T_icedot: calculates ice temperature derivative
"""
function calc_T_icedot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now_dt["T_icedot_"*hm] = now_dt["Q_dif_"*hm] + now_dt["Q_drag_"*hm]
    end
    return now_dt
end

@doc """
    calc_Tdot: calculates regional temperature derivative
"""
function calc_Tdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        # I think it is more correct to use Z instead of H -- spm 2023.01.03
        now_dt["Tdot_"*hm] = (par_dt["T_ref_"*hm] - now_dt["T_"*hm]
                             + par_dt["cs"] * now_dt["rf_"*hm] 
                             - par_dt["csz"] * now_dt["Z_"*hm]) / par_dt["tau_rf_"*hm] 
    end
    return now_dt
end

@doc """
    calc_albedodot: calculates albedo derivative
"""
function calc_albedodot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now_dt["albedodot_"*hm] = (now_dt["albedo_ref_"*hm] - now_dt["albedo_"*hm]) / par_dt["tau_albedo"]

        if now_dt["H_"*hm] == 0.0
            now_dt["albedo_"*hm] = par_dt["albedo_land"]
            now_dt["albedodot_"*hm] = 0.0
        elseif now_dt["ice_time_"*hm] < 10.0 # First ice
            now_dt["albedo_"*hm] = par_dt["albedo_newice"]
            now_dt["albedodot_"*hm] = 0.0
        end
    end
    return now_dt
end

@doc """
    calc_co2dot: calculates co2 derivative
"""
function calc_co2dot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now_dt["co2dot_"*hm] = (par_dt["co2_ref"] - now_dt["co2_"*hm] + 10.0 * (now_dt["T_"*hm] - par_dt["T_ref_"*hm])) / 10.0
    end
    return now_dt
end








