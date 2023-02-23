# =============================
#     Program: amod_update.jl
#     Aim: functions to update amod output variables and check runs
# =============================
@doc """
    update_amod_out:
        updates output variable (OUT)
"""
function update_amod_out(outf::OrderedDict, vals::OrderedDict)
    for (key, value) in outf
        push!(outf[key], vals[key])
    end
    return outf
end

@doc """
    check_run:
        checks if run is doing well
"""
function check_run(run_to_check::OrderedDict)
    counter, break_iteration = 0, false
    for (key, value) in run_to_check
        if isnan(run_to_check[key])
            printstyled("NaN value found in $(key) \n", color=:red)
            counter += 1
        end
    end
    if counter > 0
        printstyled("Run stopped at $(run_to_check["time"]) yr \n", color=:red)
        break_iteration = true
    end
    return break_iteration
end




