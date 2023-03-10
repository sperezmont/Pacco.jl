# =============================
#     Program: dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================
"""
    calc_tau_d(now, par)
calculates total driving stress through: \n
`taud = rhoi * g * H * Z / L` 

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_tau_d(now, par)
    for hm in par["hemisphere"]
        now["tau_d_"*hm] = rhoi * g * now["H_"*hm] * now["Z_"*hm] / par["L"]
    end
    return now
end

"""
    calc_tau_b(now, par)
calculates total basal stress through: \n
`taub = rhoi * g * H * Z / L` \n
It is assumed that `tau_d = tau_b` as in the SIA \n
Shallow Ice Approximation (SIA): The basal shear stress balances out completely the driving stress

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_tau_b(now, par)
    for hm in par["hemisphere"]
        now["tau_b_"*hm] = rhoi * g * now["H_"*hm] * now["Z_"*hm] / par["L"]
    end
    return now
end

"""
    calc_U_d(now, par)
calculates the driving velocity through different options \n
`"glen"` case: \n
`U_d = 2.0 * A * H * tau_d^glen_n / (glen_n + 2)`

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_U_d(now, par)
    for hm in par["hemisphere"]
        if (par["ud_case"] == "sia")
            now["U_d_"*hm] = (2.0 * now["A"] * now["H_"*hm] * (now["tau_d_"*hm]^par["glen_n"])) / (par["glen_n"] + 2)
        else
            error("ERROR, velocity option not recognized")
        end
    end
    return now
end

"""
    calc_U_b(now, par)
calculates the basal velocity: \n
`"weertman"` case: Weertman sliding law\n

Note: sliding treatment assumes temperate base in streaming areas, so no dependence on temperature just on sediments
## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_U_b(now, par)
    for hm in par["hemisphere"]
        if par["ub_case"] == "weertmanq"      # weertman quadratic 
            beta = max(0.0, min(1.0, now["Hsed_"*hm]))       # sediments Csprime weight
            Csprime = beta * par["C_s"]                  # sliding par calculation based on amount of sediments
            now["U_b_"*hm] = Csprime * now["tau_b_"*hm]^2.0            # Pollard and DeConto (2012) take m = 2 so tau_b^(m-1) * tau_b = tau_b^2  
        else
            error("ERROR, basal velocity option not recognized")
        end
    end
    return now
end

"""
    calc_fstream(now, par, ctl)
calculates stream fraction

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters
* `ctl` Run control parameters

## Return
updated `now` dictionary
"""
function calc_fstream(now, par, ctl)
    tau_kin = par["L"] / par["v_kin"]   # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)

    for hm in par["hemisphere"]
        # streaming inland propagation
        now["alpha_"*hm] = max(0.0, min(1.0, (now["T_ice_"*hm] + par["T_sb"]) / (par["T_sb"])))
        now["fstream_ref_"*hm] = (par["fstream_max_"*hm] - par["fstream_min_"*hm]) * now["alpha_"*hm] + par["fstream_min_"*hm]  # the reference value of the streaming fraction is linear with alpha (thus temperature) -- jas

        # change in the stream fraction 
        now["fstreamdot_"*hm] = (now["fstream_ref_"*hm] - now["fstream_"*hm]) / tau_kin

        # fstream
        now["fstream_"*hm] = max(now["fstream_"*hm] + now["fstreamdot_"*hm] * ctl["dt"], 0.0)
    end
    return now
end

"""
    calc_U(now, par)
calculates total velocity

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_U(now, par)
    for hm in par["hemisphere"]
        now["U_"*hm] = now["U_d_"*hm] + now["fstream_"*hm] * now["U_b_"*hm]
    end
    return now
end


#############################
# Time derivatives
#############################
"""
    update_Z(now, par)
updates ice surface elevation

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function update_Z(now::OrderedDict, par::OrderedDict)
    for hm in par["hemisphere"]
        if now["H_"*hm] == 0.0
            now["Z_"*hm] = 0.0
        elseif par["active_iso"]
            now["Z_"*hm] = max(now["H_"*hm] + now["B_"*hm],
                0.0 + (1 - (rhoi / rhow) * now["H_"*hm]))     # Pattyn 2017, Robinson 2020
        else
            now["Z_"*hm] = now["H_"*hm] + par["B_eq_"*hm]
        end
    end
    return now
end

"""
    calc_Hdot(now, par)
calculates ice thickness derivative

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_Hdot(now, par)
    for hm in par["hemisphere"]
        if par["active_ice"]
            now["Hdot_"*hm] = now["TMB_"*hm] - now["U_"*hm] * now["H_"*hm] / par["L"]
        else
            now["Hdot_"*hm] = now["TMB_"*hm] + 0.0
        end
    end
    return now
end

"""
    calc_Hseddot(now, par)
calculates sediments thickness derivative

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_Hseddot(now, par)
    for hm in par["hemisphere"]
        now["Hseddot_"*hm] = -par["f_1"] * now["U_"*hm] + par["f_2"] * now["M_"*hm] #/ ctl["dt"]
    end
    return now
end

"""
    calc_Bdot(now, par)
calculates bedrock elevation derivative

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters

## Return
updated `now` dictionary
"""
function calc_Bdot(now, par)
    for hm in par["hemisphere"]
        if par["active_iso"]
            now["Bdot_"*hm] = -(now["B_"*hm] - par["B_eq_"*hm] + now["H_"*hm] * rhoi / rhom) / par["tau_bed_"*hm] # needs further improvement -- spm 2022.11.17
        else
            now["Bdot_"*hm] = 0.0
        end
    end
    return now
end