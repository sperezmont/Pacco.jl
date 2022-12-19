# =============================
#     Program: amod_diag.jl
#     Aim: This program contains functions to calculate diagnostic variables
# =============================
@doc """
    calc_E: calculates ice extension
"""
function calc_E(now_p, par_p)
    # Northern Hemisphere
    now_p["En"] = par_p["En_ref"] * (now_p["Tn"] - now_p["Tn_ref"]) / par_p["An_te"]
    # Southern Hemisphere
    if par_p["active_antarctica"]
        now_p["Es"] = par_p["Es_ref"] * (now_p["Ts"] - now_p["Ts_ref"]) / par_p["As_te"]
    end
    return now_p
end

@doc """
    calc_V: calculates ice volume from predicted mean H and diagnosed E
"""
function calc_V(now_p, par_p)
    # Northern Hemisphere
    now_p["Vn"] = now_p["En"] * now_p["Hn"] / par_p["L"]
    # Southern Hemisphere
    if par_p["active_antarctica"]
        now_p["Vs"] = now_p["Es"] * now_p["Hs"] / par_p["L"]
    end
    return now_p
end





