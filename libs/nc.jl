# =============================
#     Program: nc.jl
#     Aim: several functions to handle netCDF files
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Load dependencies
using NCDatasets

function genout_nc(outputpath::AbstractString, filename::AbstractString, outvalues, parameters::Params{T}, attr::Dict; prec=Float64) where {T<:AbstractFloat}
    # Note: here outvalues.u includes (u, v)
    ds = NCDataset(outputpath * filename, "c")
    lt = length(outvalues.t)

    # Define time dimension
    defDim(ds, "time", lt)
    defVar(ds, "time", prec, ("time",))
    ds["time"][:] = outvalues.t

    # Define and assign the variables using global variables states_u, states_v and states_comp
    new_out_matrix = reduce(hcat, outvalues.u)' # each column has the entire time series of each variable
    for (index, variable) in enumerate(states_u)
        defVar(ds, variable, prec, ("time",), attrib=attr[variable])
        ds[variable][:] .= new_out_matrix[:, index]
    end

    # -- composite variables
    for variable in states_comp # THIS IS HARCODED; CAREFUL WITH INDEXES IF A CHANGE IS MADE IN pacco_defs.jl
        defVar(ds, variable, prec, ("time",), attrib=attr[variable])
        if variable == "Inorm"
            ds[variable][:] .= 2.0 .* (new_out_matrix[:, I_idx] .- parameters.insol_min) ./ (parameters.insol_max - parameters.insol_min) .- 1.0
        elseif variable == "Ianom"
            ds[variable][:] .= new_out_matrix[:, I_idx] .- parameters.insol_ref
        elseif variable == "m"
            ds[variable][:] .= new_out_matrix[:, s_idx] .- new_out_matrix[:, a_idx]
        elseif variable == "q"
            ds[variable][:] .= new_out_matrix[:, v_idx] .* new_out_matrix[:, H_idx] ./ parameters.L
        elseif variable == "sw"
            ds[variable][:] .= max.(0.0, parameters.kI .* (1 .- new_out_matrix[:, albedo_eff_idx]) .* (new_out_matrix[:, I_idx] .- parameters.insol_ref))
        elseif variable == "lw"
            ds[variable][:] .= max.(0.0, parameters.lambda .* (new_out_matrix[:, T_idx] .- parameters.Tref))
        elseif variable == "RCO2"
            ds[variable][:] .= 5.35 .* log.(new_out_matrix[:, C_idx] ./ 280.0)
        end
    end
    close(ds)
end

function load_nc(filename::AbstractString, vars::Vector; time_name="time")
    if length(vars) == 1
        return NCDataset(filename)[vars[1]], NCDataset(filename)[time_name]
    else
        return [NCDataset(filename)[v] for v in vars], NCDataset(filename)[time_name]
    end
end

