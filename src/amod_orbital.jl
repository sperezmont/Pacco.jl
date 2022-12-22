# =============================
#     Program: amod_orbital.jl
#     Aim: This program contains functions to calculate orbital parameters
# =============================
@doc """
    calc_artificial_insolation: Compute daily average insolation through different parameterizations
"""
function calc_artificial_insolation(now_o, par_o)
    # Calculate reference and amplitude
    ins_ref = (par_o["ins_max"] + par_o["ins_min"]) / 2
    A_ins = (par_o["ins_max"] - par_o["ins_min"]) / 2

    # Return artificial insolation -- I have to discuss this with jas
    ins = ins_ref + A_ins * (
        par_o["P_obl"] * cos(2.0 * pi * now_o["time"] / par_o["tau_obl"]) +
        par_o["P_pre"] * cos(2.0 * pi * now_o["time"] / par_o["tau_pre"]) +
        par_o["P_exc"] * cos(2.0 * pi * now_o["time"] / par_o["tau_exc"]))
    for hm in hemisphere
        now_o["ins_"*hm] = ins
    end
    return now_o
end

@doc """
    calc_laskar_insolation: Compute daily average insolation given latitude and time of year
"""
function calc_laskar_insolation(now_o, par_o)
    date::DateTime = DateTime(now_o["time"], par_o["ins_month"], par_o["ins_day"])
    for hm in hemisphere
        zen_ang, es_dist = daily_zenith_angle(date, par_o["ins_lat_"*hm], param_set)
        now_o["ins_"*hm] = insolation(zen_ang, es_dist, param_set)
    end
    return now_o
end

@doc """
    calc_insol_day: Compute daily average insolation
"""
function calc_ins(now_o, par_o)
    if par_o["ins_case"] == "artificial"
        now_o = calc_artificial_insolation(now_o, par_o)
    elseif par_o["ins_case"] == "laskar"
        now_o = calc_laskar_insolation(now_o, par_o)
    else
        error("ERROR, insolation option not recognized")
    end
    return now_o
end

