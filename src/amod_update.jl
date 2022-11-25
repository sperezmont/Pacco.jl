# =============================
#     Program: amod_update.jl
#     Aim: functions to update amod variables
# =============================
function update_amod_out(d::OrderedDict, vals::OrderedDict)
    for (key, value) in vals
        push!(d[key], vals[key])
    end
    return d
end


