# =============================
#     Program: analysis.jl
#     Aim: analyze results from AMOD simulations
#     Author: Sergio Pérez-Montero, 2023.01.26
# =============================


@doc """
    interp_2series:
        Takes d2 and interpoles it to d2 grid to make them comparable
        Both time series should start and end at similar times!!
"""
function interp_2series(d1, d2)
    ld = length(d2)
    # create grid and interpolant
    grid = 1:ld
    d_interp = linear_interpolation(grid, d2)

    # calculate new grid
    new_grid = range(1, stop=ld, length=length(d1))

    # return new time series
    new_d2 = d_interp.(new_grid)
    return new_d2
end

@doc """
    compute_comp_stats:
        Takes vectors d1 and d2 and compute some comparison-oriented statistics
"""
function compute_comp_stats(d1::Vector, d2::Vector)
    return Dict("dif" => d2 .- d1,
        "corr" => cor(d1, d2),
        "spect_dif"=> nothing)
end

@doc """
    compute_stats:
        Takes a vector a compute some statistics
"""
function compute_stats(d::Vector)
    return Dict("min" => minimum(skipmissing(d)),
        "max" => maximum(skipmissing(d)),
        "mean" => mean(skipmissing(d)),
        "vari" => var(skipmissing(d)))
end

@doc """
    analyze_amod:
        Analyzes a simulation or ensemble against available proxies
        experiment      --> experiment name to analyze
        isens           --> is it an ensemble?
"""
function analyze_amod(; experiment::String="test_default_ens", isens::Bool=true)
    # Define some local variables
    locdir = pwd() * "/output/" * experiment * "/"

    # Check if experiment is an ensemble and extract subdirectories
    elements = []
    for element in readdir(locdir)
        if element[end-3:end] != ".png"
            push!(elements, element)
        end
    end

    # Check available proxies
    data_list, proxy_data, vars2compare = readdir(pwd() * "/data/"), Dict(), []
    for i in 1:length(data_list)
        element_name = data_list[i]
        if element_name[end-2:end] == ".nc" # -- filt netCDF files
            df = NCDataset(pwd() * "/data/" * element_name) # -- load data frame
            varnames = keys(df)

            proxy_data[element_name[1:end-3]] = Dict(v => df[v][:] for v in varnames)   # -- store info in dictionary

            varidx = 2
            for letter in eachindex(element_name[1:end-3])
                if element_name[1:end-3][letter] == '_'
                    varidx = letter - 1 
                    break
                end
            end
            if element_name[1:varidx] in vars2compare
                continue
            else
                push!(vars2compare, element_name[1:varidx])
            end
        end
    end

    # Load certain results from the ensemble (comparable variables to proxy info)
    amod_data = Dict()
    vars2load = ["time"; vars2compare .* "_n"]
    vars2compare_aux = ["time"; vars2compare]
    for i in eachindex(elements)
        df = NCDataset(locdir * elements[i] * "/amod.nc")
        amod_data[elements[i]] = Dict(vars2compare_aux[v] => df[vars2load[v]][:] for v in eachindex(vars2load))
    end

    # Cut time series, assume all amod runs have the same time length
    first_ts = [abs(proxy_data[k]["time"][1]) for k in keys(proxy_data)] 
    proxy_min_t, amod_min_t = findmin(first_ts), abs(amod_data[elements[1]]["time"][1])
    min_t = min(proxy_min_t[1], amod_min_t)
    if min_t == proxy_min_t[1]
        min_name = collect(keys(proxy_data))[proxy_min_t[2]]
        new_t1, new_tend = proxy_data[min_name]["time"][1], proxy_data[min_name]["time"][end]
    else
        min_name = elements[1]
        new_t1, new_tend = amod_data[elements[1]]["time"][1], amod_data[elements[1]]["time"][end]  
    end
    
    for k in keys(proxy_data)
        if k != min_name
            idx = findmin(abs.(proxy_data[k]["time"] .- new_t1))[2]
            for k2 in keys(proxy_data[k])
                proxy_data[k][k2] = proxy_data[k][k2][idx:end]
            end 
        end
    end

    for k in keys(amod_data)
        if k != min_name
            idx = findmin(abs.(amod_data[k]["time"] .- new_t1))[2]
            for k2 in keys(amod_data[k])
                amod_data[k][k2] = amod_data[k][k2][idx:end]
            end 
        end
    end

    # Compute Statistics
    # -- proxy
    proxy_stats = Dict(v => Dict() for v in keys(proxy_data))
    for i in keys(proxy_data)
        stats_i = Dict(v => Dict() for v in keys(proxy_data[i]))
        for j in keys(proxy_data[i])
            if j != "time"
                stats_ij = compute_stats(proxy_data[i][j])
                stats_i[j] = stats_ij
            end
        end
        delete!(stats_i, "time")
        proxy_stats[i] = stats_i
    end

    # -- amod ensemble
    # ---- ensemble stats
    amod_stats = Dict(v => Dict() for v in keys(amod_data))
    for i in keys(amod_data)
        stats_i = Dict(v => Dict() for v in keys(amod_data[i]))
        for j in keys(amod_data[i])
            if j != "time"
                if j == "T"
                    amod_data[i][j] = amod_data[i][j] .- 273.15 # t_ref
                end
                stats_ij = compute_stats(amod_data[i][j])
                stats_i[j] = stats_ij
            end
        end
        delete!(stats_i, "time")
        amod_stats[i] = stats_i
    end

    # ---- ensemble vs proxy
    comp_data = Dict()
    for kp in keys(proxy_data)
        stats_kp = Dict()
        for v in vars2compare
            stats_kpv = Dict()
            if v in keys(proxy_data[kp])
                for ka in keys(amod_data)
                    d1, d2 = proxy_data[kp][v], amod_data[ka][v]
                    new_d2 = interp_2series(d1, d2)
                    stats_kpv[ka] = compute_comp_stats(d1, new_d2)
                end
                stats_kp[v] = stats_kpv
            end
        end
        comp_data[kp] = stats_kp
    end

    # Plot results, 1 file for each variable
    stat_vars = ["min", "mean", "max", "vari"]
    # -- convert to array stats data
    d = Array{Real}(undef, length(elements), length(vars2compare), length(stat_vars))
    for i in eachindex(elements)
        for j in eachindex(vars2compare)
            for k in eachindex(stat_vars)
                d[i, j, k] = amod_stats[elements[i]][vars2compare[j]][stat_vars[k]]
            end
        end
    end

    # -- Figure 1, time series and stats
    fig = Figure(resolution=(1000, 500))
    prx_clr = [:navy, :firebrick, :olive, :indigo]
    for i in eachindex(vars2compare)
        # ---- time series
        ax = Axis(fig[i, 1], ylabel=vars2compare[i])
        for j in eachindex(elements)
            d2p = amod_data[elements[j]][vars2compare[i]]
            lines!(ax, amod_data[elements[j]]["time"], d2p, color="grey")
        end
        k = 1
        for p in 1:length(keys(proxy_data))
            prxnm = collect(keys(proxy_data))[p]
            if vars2compare[i] in keys(proxy_data[prxnm])
                lines!(ax, proxy_data[prxnm]["time"], proxy_data[prxnm][vars2compare[i]], color=prx_clr[k])
                k += 1
            end
        end
        xlims!(ax, amod_data[elements[1]]["time"][1], amod_data[elements[1]]["time"][end])

        # ---- ensemble stats
        ax_stats = Axis(fig[i, 2])
        for j in eachindex(stat_vars[1:3])  # only min, mean and max
            boxplot!(ax_stats, j .* ones(length(elements)), d[:, i, j], color="grey")
            k = 1
            for p in 1:length(keys(proxy_stats))
                prxnm = collect(keys(proxy_stats))[p]
                if vars2compare[i] in keys(proxy_stats[prxnm])
                    if j == 1
                        scatter!(ax_stats, j, proxy_stats[prxnm][vars2compare[i]][stat_vars[j]], color=prx_clr[k])
                    else
                        scatter!(ax_stats, j, proxy_stats[prxnm][vars2compare[i]][stat_vars[j]], color=prx_clr[k])
                    end
                    k += 1
                end
            end
        end
        linkyaxes!(ax, ax_stats)
        ax_stats.xtickformat = k -> stat_vars[1:3]
        ax_stats.yticklabelsvisible = false

        # ---- variance
        ax_vari = Axis(fig[i, 3])
        boxplot!(ax_vari, ones(length(elements)), d[:, i, end], color="grey", label="AMOD")
        k = 1
        for p in 1:length(keys(proxy_stats))
            prxnm = collect(keys(proxy_stats))[p]
            if vars2compare[i] in keys(proxy_stats[prxnm])
                scatter!(ax_vari, 1, proxy_stats[prxnm][vars2compare[i]]["vari"], color=prx_clr[k], label=prxnm[length(vars2compare[i]*"_")+1:end])
                k += 1
            end
        end

        fig[i, 4] = Legend(fig, ax_vari, framevisible=false, labelsize=15)

        ax_vari.xticks = ([1.0], ["variance"])

        if i != length(vars2compare)
            ax.xticklabelsvisible = false
            ax_stats.xticklabelsvisible = false
            ax_vari.xticklabelsvisible = false
        end
    end

    colsize!(fig.layout, 2, Relative(2/10))
    colsize!(fig.layout, 3, Relative(2/10))
    colsize!(fig.layout, 4, Relative(2/10))
    save(locdir * "time-series_ens_stats.png", fig)

    # -- Figure 2, comparative stats
    fig = Figure(resolution=(1000, 800))
    prx_clr = [:navy, :firebrick, :olive, :indigo]
    for i in eachindex(collect(keys(comp_data)))
        proxy_label = collect(keys(comp_data))[i]
        ax = Axis(fig[i, 1], ylabel=proxy_label)
        idx = 1
        for l in 1:length(proxy_label)
            if proxy_label[l] == '_'
                idx = l
                break
            end
        end
        var2plot = proxy_label[1:idx-1]
        for j in eachindex(elements)
            lines!(ax, comp_data[proxy_label][var2plot][elements[j]]["dif"], color="grey")
        end

        # ---- ensemble correlation boxplot
        ax_cor = Axis(fig[i, 2])
        cor_ar = [comp_data[proxy_label][var2plot][e]["corr"] for e in elements]
        max_cor = findmax(cor_ar)   
        lines!(ax, comp_data[proxy_label][var2plot][elements[max_cor[2]]]["dif"], color="green", label=elements[max_cor[2]])
        boxplot!(ax_cor, ones(length(elements)), cor_ar, color="grey")

        ax_cor.xticks = ([1.0], ["correlation"])
        axislegend(ax, framevisible=false, position=:rb)

        if i != length(collect(keys(comp_data)))
            ax.xticklabelsvisible = false
            ax_cor.xticklabelsvisible = false
        end

    end
    colsize!(fig.layout, 2, Relative(1/10))
    save(locdir * "comp-stats_ens.png", fig)
end

