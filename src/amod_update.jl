# =============================
#     Program: amod_update.jl
#     Aim: functions to update amod variables
# =============================
@doc """
    update_amod_out: updates output variable (OUT)
"""
function update_amod_out(outf::OrderedDict, vals::OrderedDict)
    for (key, value) in outf
        push!(outf[key], vals[key])
    end
    return outf
end

@doc """
    update_forward: takes vars2update and calculates their time evolution
                    vars2update is a vector with the names of the variables to update
                    vars2update is computed in run_amod() function 
"""
function update_forward(now_u, par_u, ctl_u, vars2update)
    for hm in par_u["hemisphere"], v in vars2update
        # -- calculate time evolution
        variab, vardot = v * "_" * hm, v * "dot_" * hm
        now_u[variab] += now_u[vardot] * ctl_u["dt"]    # now = now + dnow/dt * dt

        # -- modify if desired
        if variab in ["H_n", "H_s"]
            now_u[variab] = max(now_u[variab], 0.0)
        elseif variab in ["Hsed_n", "Hsed_s"]
            now_u[variab] = min(max(now_u[variab], 0.0), 1.0) # sediments go from 0 to 1
        elseif variab in ["B_n", "B_s"]
            (~par_u["active_iso"]) && (now_u[variab] = par_u["B_eq_"*hm]) # reupdate to equilibrium value
        elseif variab in ["T_ice_n", "T_ice_s"]
            now_u[variab] = min(now_u[variab], degK)
        elseif variab in ["co2_n", "co2_s"]
            now_u[variab] = max(now_u[variab], 0.0)
        end
    end
    return now_u
end

@doc """
    update_Z: updates ice surface elevation
"""
function update_Z(now_u, par_u)
    for hm in par_u["hemisphere"]
        if par_u["active_iso"]
            now_u["Z_"*hm] = max(now_u["H_"*hm] + now_u["B_"*hm],
                0.0 + (1 - (rhoi / rhow) * now_u["H_"*hm]))     # Pattyn 2017, Robinson 2020
        else
            now_u["Z_"*hm] = now_u["H_"*hm] + par_u["B_eq_"*hm]
        end
    end
    return now_u
end
