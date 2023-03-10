# =============================
#     Program: pacco_update.jl
#     Aim: functions to update pacco output variables and check runs
# =============================
"""
    update_pacco_out(now, out)
updates output variable

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `out` Dictionary with stored values of the model variables (output)

## Return
updated `out` variable
"""
function update_pacco_out(now::OrderedDict, out::OrderedDict)
    for (key, value) in out
        push!(out[key], now[key])
    end
    return out
end

"""
    check_run(run_to_check)
checks if run is doing well

## Arguments
* `run_to_check` Dictionary with values to check

## Return
`true`/`false` 
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




