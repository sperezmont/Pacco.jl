# =============================
#     Program: amod_derivative.jl
#     Aim: This program contains functions for calculating derivatives
# =============================
@doc """
    calc_Hdot: calculates ice thickness derivative
"""
function calc_Hdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now_dt["Hdot"*hm] = 0.0
        # -- climatic term
        if par_dt["active_climate"]
            now_dt["Hdot"*hm] = copy(now_dt["TMB"])
        end
        # -- dynamical term
        if par_dt["active_dynamics"]
            now_dt["Hdot"*hm] = now_dt["Hdot"*hm] 
                              - now_dt["U"*hm] * now_dt["H"*hm] / par_dt["L"]
        end
    end
    return now_dt
end

@doc """
    calc_Hseddot: calculates sediments thickness derivative
"""
function calc_Hdot(now_dt, par_dt, ctl_dt)
    for hm in par_dt["hemisphere"]
        now_dt["Hseddot"*hm] = -par_dt["f_1"] * now_dt["U"*hm] + par_dt["f_2"] * now_dt["M"*hm] / ctl_dt["dt"]
    end
    return now_dt
end

@doc """
    calc_Bdot: calculates bedrock elevation derivative
"""
function calc_Bdot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        if par_dt["active_iso"]
            now_dt["Bdot"*hm] = -(now_dt["B"*hm] - par_dt["B_eq"*hm] + now_dt["H"*hm] * rhoi / rhom) / par["tau_bed"*hm] # needs further improvement -- spm 2022.11.17
        else
            now["Bdot"*hm] = 0.0
        end
    end
    return now_dt
end

@doc """
    calc_T_icedot: calculates ice temperature derivative
"""
function calc_T_icedot(now_dt, par_dt)
    for hm in par_dt["hemisphere"]
        now["T_icedot"*hm] = now["Q_dif"*hm] + now["Q_drag"*hm]
    end
    return now_dt
end

@doc """
    calc_Tdot: calculates regional temperature derivative
"""
function calc_Tdot(now_dt, par_dt)
    for hm in par_d["hemisphere"]
        now_dt["Tdot"*hm] = (par_dt["T_ref"*hn] - now_dt["T"*hn] + par_dt["cs"]*now_dt["rf"*hm] - grad*now_dt["H"*hm]) / par_dt["tau_rf"*hm] # I have to test grad*now_dt["Z"*hm] -- spm 2022.12.20
    end
    return now_dt
end

@doc """
    calc_albedodot: calculates albedo derivative
"""
function calc_albedodot(now_dt, par_dt)
    for hm in par_d["hemisphere"]
        now_dt["albedodot"*hm] = (now_dt["albedo_ref"*hm] - now_dt["albedo"*hm]) / par_dt["tau_albedo"*hm]
        
        if now_dt["H"*hm] == 0.0
            now_dt["albedo"*hm] = par_dt["albedo_land"]
            now_dt["albedodot"*hm] = 0.0
        elseif now_dt["ice_time"*hm] < 10.0 # First ice
            now_dt["albedo"*hm] = par_dt["albedo_newice"*hm]
            now_dt["albedodot"*hm] = 0.0                
        end
    end
    return now_dt
end

@doc """
    calc_co2dot: calculates regional temperature derivative
"""
function calc_co2dot(now_dt, par_dt)
    for hm in par_d["hemisphere"]
        now_dt["co2dot"*hm] = (par_dt["co2_ref"] - now_dt["co2"*hm] + 10.0 * (now_dt["T"*hm] - par_dt["T_ref"*hm])) / 10.0
    end
    return now_dt
end








