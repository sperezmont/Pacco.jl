# =============================
#     Program: amod_dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================
@doc """
    calc_tau_d:
        calculates total driving stress through
        taud = rhoi * g * H * Z / L 
"""
function calc_tau_d(now_d, par_d)
    for hm in par_d["hemisphere"]
        now_d["tau_d_"*hm] = rhoi * g * now_d["H_"*hm] * now_d["Z_"*hm] / par_d["L"]
    end
    return now_d
end

@doc """
    calc_tau_b:
        calculates total basal stress through
        taub = rhoi * g * H * Z / L 
        It is assumed that tau_d = tau_b as in the SIA.     -- jas
        SIA --> The basal shear stress balances out completely the driving stress

"""
function calc_tau_b(now_d, par_d)
    for hm in par_d["hemisphere"]
        now_d["tau_b_"*hm] = rhoi * g * now_d["H_"*hm] * now_d["Z_"*hm] / par_d["L"]
    end
    return now_d
end

@doc """
    calc_U_d:
        calculates the driving velocity through different options
            --> glen case 
                    U_d = 2.0 * A * H * tau_d^glen_n / 5.0

"""
function calc_U_d(now_d, par_d)
    for hm in par_d["hemisphere"]
        if (par_d["ud_case"] == "sia")
            now_d["U_d_"*hm] = (2.0 * now_d["A"] * now_d["H_"*hm] * now_d["tau_d_"*hm]^par_d["glen_n"]) / (par_d["glen_n"] + 2)
        else
            error("ERROR, velocity option not recognized")
        end
    end
    return now_d
end

@doc """
    calc_U_b: 
        calculates the basal velocity taking 
            --> case weertman
                    Weertman sliding law
            
        sliding treatment --> I assume temperate base in streaming areas, so no dependence on temperature just on sediments -- jas
"""
function calc_U_b(now_d, par_d)
    for hm in par_d["hemisphere"]
        if par_d["ub_case"] == "weertmanq"      # weertman quadratic 
            beta = max(0.0, min(1.0, now_d["Hsed_"*hm]))       # sediments Csprime weight
            Csprime = beta * par_d["C_s"]                  # sliding par calculation based on amount of sediments
            now_d["U_b_"*hm] = Csprime * now_d["tau_b_"*hm]^2.0            # Pollard and DeConto (2012) take m = 2 so tau_b^(m-1) * tau_b = tau_b^2  
        else
            error("ERROR, basal velocity option not recognized")
        end
    end
    return now_d
end

@doc """
    calc_fstream:
        calculates stream fraction 
"""
function calc_fstream(now_d, par_d, ctl_d)
    tau_kin = par_d["L"] / par_d["v_kin"]   # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)

    for hm in par_d["hemisphere"]
        # streaming inland propagation
        now_d["alpha_"*hm] = max(0.0, min(1.0, (now_d["T_ice_"*hm] + par_d["T_sb"]) / (par_d["T_sb"])))
        now_d["fstream_ref_"*hm] = (par_d["fstream_max_"*hm] - par_d["fstream_min_"*hm]) * now_d["alpha_"*hm] + par_d["fstream_min_"*hm]  # the reference value of the streaming fraction is linear with alpha (thus temperature) -- jas

        # change in the stream fraction 
        now_d["fstreamdot_"*hm] = (now_d["fstream_ref_"*hm] - now_d["fstream_"*hm]) / tau_kin

        # fstream
        now_d["fstream_"*hm] = max(now_d["fstream_"*hm] + now_d["fstreamdot_"*hm] * ctl_d["dt"], 0.0)
    end
    return now_d
end

@doc """
    calc_U:
        calculates total velocity
"""
function calc_U(now_d, par_d)
    for hm in par_d["hemisphere"]
        now_d["U_"*hm] = now_d["U_d_"*hm] + now_d["fstream_"*hm] * now_d["U_b_"*hm]
    end
    return now_d
end