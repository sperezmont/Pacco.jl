# =============================
#     Program: nc.jl
#     Aim: several functions to handle netCDF files
#     Author: Sergio Pérez-Montero, 2022.11.11
# =============================

# Load dependencies
using NCDatasets
using DataStructures

function genout_nc(out::AbstractString, filename::AbstractString, d::OrderedDict, precc, attr::OrderedDict)
    ds = NCDataset(out * filename, "c")

    # define time dimension
    defDim(ds, "time", Inf)

    # define groups # NOT IMPLEMENTED YET -- spm 2022.10.27
    # for i in eachindex(out_groups)
    #     defGroup(ds, out_groups[i], attrib = Dict("title" => out_groups[i]))
    # end

    # define the variables
    for (key, val) in d
        defVar(ds, key, precc, ("time",), attrib=attr[key])
    end

    # assign values
    for (key, val) in d
        try # only assign value if it is defined in amod_defs (out_attr)
            ds[key][:] = d[key]
        catch
            continue
        end
    end

    close(ds)
end

function load_nc(filename::AbstractString, vars::Vector)
    if length(vars) == 1
        return NCDataset(filename)[vars[1]], NCDataset(filename)["time"]
    else
        return [NCDataset(filename)[v] for v in vars], NCDataset(filename)["time"]
    end
end

