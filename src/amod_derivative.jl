# =============================
#     Program: amod_derivative.jl
#     Aim: This program contains functions for calculating derivatives
# =============================
@doc """
    calc_Hdot:
        calculates ice thickness derivative
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
    calc_Hseddot:
        calculates sediments thickness derivative
"""
function calc_Hseddot(now_dt, par_dt, ctl_dt)
    for hm in par_dt["hemisphere"]
        now_dt["Hseddot_"*hm] = -par_dt["f_1"] * now_dt["U_"*hm] + par_dt["f_2"] * now_dt["M_"*hm] / ctl_dt["dt"]
    end
    return now_dt
end

@doc """
    calc_Bdot:
        calculates bedrock elevation derivative
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
    calc_T_icedot:
        calculates ice temperature derivative
"""
function calc_T_icedot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now_dt["T_icedot_"*hm] = now_dt["Q_dif_"*hm] + now_dt["Q_drag_"*hm]
    end
    return now_dt
end

@doc """
    calc_Tdot:
        calculates regional temperature derivative
"""
function calc_Tdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]

        # Anthropogenic forcing?
        if now_dt["time"] < par_dt["time_anth"] # unperturbed climate
            Temp_ref = par_dt["T_ref_"*hm]
        else    # perturbed climate
            Temp_ref = par_dt["T_ref_"*hm] + par_dt["cco2"] * now_dt["co2_"*hm] / par_dt["co2_ref"]
        end

        if par_dt["height_temp"] == "useH"
            HTF = now_dt["H_"*hm]
        elseif par_dt["height_temp"] == "useZ"
            HTF = now_dt["Z_"*hm]
        else
            printstyled("dev par must be removed!", color=:red)
        end

        now_dt["Tdot_"*hm] = (Temp_ref - now_dt["T_"*hm]
                              +
                              now_dt["T_rf_"*hm]
                              -
                              par_dt["csz"] * HTF) / par_dt["tau_rf_"*hm]
    end
    return now_dt
end

@doc """
    calc_albedodot:
        calculates albedo derivative
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
    calc_co2dot:
        calculates co2 derivative
"""
function calc_co2dot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]

        # -- anthropogenic forcing?
        if now_dt["time"] < par_dt["time_anth"] # unperturbed climate
            co2ref = par_dt["co2_ref"]
        else    # perturbed climate
            actual_diftime = now_dt["time"] - par_dt["time_anth"]
            co2ref = par_dt["co2_ref"] + par_dt["co2_anth"] * (
                0.75 / exp(actual_diftime / 365.0)
                + 0.135 / exp(actual_diftime / 5500.0)
                + 0.035 / exp(actual_diftime / 8200.0)
                + 0.08 / exp(actual_diftime / 200000.0)
            )   # Archer 1998
        end

        now_dt["co2dot_"*hm] = (co2ref - now_dt["co2_"*hm] + par_dt["ktco2"] * (now_dt["T_"*hm] - par_dt["T_ref_"*hm])) / par_dt["tau_co2"]
    end
    return now_dt
end








