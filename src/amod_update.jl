# =============================
#     Program: amod_update.jl
#     Aim: functions to update amod variables
# =============================
@doc """
    update_amod_out: updates output variable (OUT)
"""
function update_amod_out(d::OrderedDict, vals::OrderedDict)
    for (key, value) in vals
        push!(d[key], vals[key])
    end
    return d
end

@doc """
    update_forward: takes vars2update and calculates their time evolution
                    vars2update is a vector with the names of the variables to update
                    vars2update is computed in run_amod() function 
"""
function update_forward(now_u, par_u, ctl_u, vars2update)
    for hm in par_u["hemisphere"], v in vars2update
        # -- calculate time evolution
        variab, vardot = v*hm, v*"dot"*hm
        now_u[variab] += now_u[vardot] * ctl_u["dt"]    # now = now + dnow/dt * dt

        # -- modify if desired
        if variab in ["Hn", "Hs"]
            now_u[variab] = max(now_u[variab], 0.0)
        elseif variab in ["Hsedn", "Hseds"]
            now_u[variab] = min(max(now_u[variab], 0.0), 1.0) # sediments go from 0 to 1
        elseif variab in ["Bn", "Bs"]
            (~par_u["active_iso"]) && (now_u[variab] = par_u["B_eq"*hm]) # reupdate to equilibrium value
        elseif variab in ["T_icen", "T_ices"]
            now_u[variab] = min(now_u[variab], degK)
        end
    end 
    return now_u
end

@doc """
    update_Z: updates ice surface elevation
"""
function update_Z(now_u, par_u)
    for hm in par_u["hemisphere"], v in vars2update
        if par_u["active_iso"]
            now_u["Z"*hm] = max(now_u["H"*hm] + now_u["B"*hm], 
                        0.0 + (1 - (rhoi / rhow) * now_u["H"*hm]))     # Pattyn 2017, Robinson 2020
        else
            now_u["Z"*hm] = now_u["H"*hm] + par_u["B_eq"*hm]
        end
    end
    return now_u
end
