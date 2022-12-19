# =============================
#     Program: amod_dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================
@doc """
    calc_taud: calculates total driving stress through
        taud = rhoi * g * H * Z / L 
"""
function calc_taud(now_d, par_d)
    for hm in par_d["hemisphere"]
        now_d["tau_d"*hm] = rhoi * g * now_d["H"*hm] * now_d["Z"*hm] / par_d["L"]
    end
    return now_d
end

@doc """
    calc_taub: calculates total basal stress through
        taub = rhoi * g * H * Z / L 
        It is assumed that tau_d = tau_b as in the SIA.     -- jas
        SIA --> The basal shear stress balances out completely the driving stress

"""
function calc_taub(now_d, par_d)
    for hm in par_d["hemisphere"]
        now_d["tau_b"*hm] = rhoi * g * now_d["H"*hm] * now_d["Z"*hm] / par_d["L"]
    end
    return now_d
end

@doc """
    calc_Ud: calculates the driving velocity through different options
        --> glen case 
                U_d = 2.0 * A * H * tau_d^glen_n / 5.0

"""
function calc_Ud(now_d, par_d)
    for hm in par_d["hemisphere"]
        if (par_d["ud_case"] == "sia")
            now_d["U_d"*hm] = (2.0 * now_d["A"*hm] * now_d["H"*hm] * now_d["tau_d"*hm]^par_d["glen_n"]) / (par_d["glen_n"] + 2)
        else
            error("ERROR, velocity option not recognized")
        end
    end
    return now_d
end

@doc """
    calc_Ub: calculates the basal velocity taking 
        --> case weertman
                Weertman sliding law
        
    sliding treatment --> I assume temperate base in streaming areas, so no dependence on temperature just on sediments -- jas
"""
function calc_Ub(now_d, par_d)
    for hm in par_d["hemisphere"]
        if par_d["ub_case"] == "weertmanq"      # weertman quadratic 
            beta = max(0.0, min(1.0, now_d["Hsed"*hm]))       # sediments Csprime weight
            Csprime = beta * par_d["C_s"]                  # sliding par calculation based on amount of sediments
            now_d["U_b"*hm] = Csprime * now_d["tau_b"*hm]^2.0            # Pollard and DeConto (2012) take m = 2 so tau_b^(m-1) * tau_b = tau_b^2  
        else
            error("ERROR, basal velocity option not recognized")
        end
    end
    return now_d
end

@doc """
    calc_fstream: calculates stream fraction 
"""
function calc_fstream(now_d, par_d, ctl_d)
    tau_kin = par["L"] / par["v_kin"]   # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)

    for hm in par_d["hemisphere"]
        # streaming inland propagation
        now_d["alpha"*hm] = max(0.0, min(1.0, (now_d["T_ice"*hm] + par_d["T_sb"]) / (par_d["T_sb"])))
        now_d["fstream_ref"*hm] = (par_d["fstream_max"*hm] - par_d["fstream_min"*hm]) * now_d["alpha"*hm] + par_d["fstream_min"*hm]  # the reference value of the streaming fraction is linear with alpha (thus temperature) -- jas

        # change in the stream fraction 
        now_d["fstreamdot"*hm] = (now_d["fstream_ref"*hm] - now_d["fstream"*hm]) / tau_kin
        
        # fstream
        now_d["fstream"*hm] = max(now_d["fstream"*hm] + now_d["fstreamdot"*hm] * ctl_d["dt"], 0.0)
    end
    return now_d
end

@doc """
    calc_U: calculates total velocity
"""
function calc_U(now_d, par_d)
    for hm in par_d["hemisphere"]
        now["U_d"*hm] + now["fstream"*hm] * now["U_b"*hm]
    end
    return now_d
end