# =============================
#     Program: amod_dynamics.jl
#     Aim: This program contains functions to calculate dynamics
# =============================

@doc """
    calc_taud: calculates total driving stress through
        taud = rho * g * H * S / L 
"""
function calc_taud(now_d, par_d)
    return rho * g * now_d["H"] * now_d["S"] / par_d["L"]
end

@doc """
    calc_taub: calculates total basal stress through
        taub = rho * g * H * S / L 
        It is assumed that tau_d = tau_b as in the SIA.     -- jas
        SIA --> The basal shear stress balances out completely the driving stress

"""
function calc_taub(now_d, par_d)
    return rho * g * now_d["H"] * now_d["S"] / par_d["L"]
end

@doc """
    calc_Ud: calculates the driving velocity through different options
        --> glen case 
                U_d = 2.0 * A * H * tau_d^glen_n / 5.0

"""
function calc_Ud(now_d, par_d)
    if (par_d["ud_case"] == "sia")
        return 2.0 * now_d["A"] * now_d["H"] * now_d["tau_d"]^par_d["glen_n"] / (par_d["glen_n"] + 2)
    else
        error("ERROR, velocity option not recognized")
    end
end

@doc """
    calc_Ub: calculates the basal velocity taking 
        --> case weertman
                Weertman sliding law
        
    sliding treatment --> I assume temperate base in streaming areas, so no dependence on temperature just on sediments -- jas
"""
function calc_Ub(now_d, par_d)
    if par_d["ub_case"] == "weertmanq"      # weertman quadratic 
        beta = max(0.0, min(1.0, now_d["Hsed"]))       # sediments Csprime weight
        Csprime = beta * par_d["C_s"]                  # sliding par calculation based on amount of sediments
        return Csprime * now_d["tau_b"]^2.0            # Pollard and DeConto (2012) take m = 2 so tau_b^(m-1) * tau_b = tau_b^2  
    else
        error("ERROR, basal velocity option not recognized")
    end
end

@doc """
    calc_fstreamdot: calculates stream fraction change 
"""
function calc_fstreamdot(now_d, par_d, tau_kin)
    # streaming inland propagation
    now_d["alpha"] = max(0.0, min(1.0, (now_d["T"] + par_d["t_sb"] + degK) / (par_d["t_sb"] + degK)))
    now_d["fstream_ref"] = (par_d["fstream_max"] - par_d["fstream_min"]) * now_d["alpha"] + par_d["fstream_min"]  # the reference value of the streaming fraction is linear with alpha (thus temperature) -- jas

    # change in the stream fraction 
    return (now_d["fstream_ref"] - now_d["fstream"]) / tau_kin
end

