# =============================
#     Program: analysis.jl
#     Aim: analyze results from AMOD simulations
#     Author: Sergio Pérez-Montero, 2023.01.26
# =============================

@doc """


"""
function calc_spectrum(d, fs)
    # -- spectrum blackman tuckey
    N = length(d)
    Nmax = Int(ceil(N / 2))
    P = periodogram(d, fs=fs, window=blackman(N))
    G, freq = P.power, P.freq
    return G, freq
end

@doc """
    calc_wavelet:
        Generate an array with the values of the wavelet applied to d
"""
function calc_wavelet(d, fs; sigma=π)
    wvt = ContinuousWavelets.wavelet(Morlet(sigma), s=8, boundary=ZPBoundary(), averagingType=NoAve(), β=1)
    S = ContinuousWavelets.cwt(d, wvt)                                           
    freq= getMeanFreq(ContinuousWavelets.computeWavelets(length(d), wvt)[1], fs)  
    S = abs.(S) .^ 2 
    S = S ./ sum(S, dims=2)

    # # -- cone of influence
    # s = 8
    # taus = sqrt(2) * s
    # coi = 1 ./ exp(t ./ taus)
    return S, freq
end

@doc """
    interp_2series:
        Takes d2 and interpoles it to d2 grid to make them comparable
        Both time series should start and end at similar times!!
"""
function interp_2series(d1::Vector, d2::Vector)
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
    calc_spect_dif:
        Compute the difference between the spectra of d1 and d2
"""
function calc_spect_dif(d1::Vector, d2::Vector, fs::Real)
    # compute spectrum
    G1, freq = calc_spectrum(d1, fs)
    G2, freq = calc_spectrum(d2, fs)

    # filt very low frequencies
    ftocut = 2e5
    G1, G2 = G1[freq .> 1 / ftocut], G2[freq .> 1 / ftocut]
    freq = freq[freq .> 1 / ftocut]

    # Normalize
    G1_norm, G2_norm = G1 ./ sum(G1), G2 ./ sum(G2)

    # Compute periods
    period = 1 ./ freq

    # Return difference and periods
    difG = G2_norm .- G1_norm
    return (difG, period)
end

@doc """
    compute_comp_stats:
        Takes vectors d1 and d2 and compute some comparison-oriented statistics
"""
function compute_comp_stats(d1::Vector, d2::Vector, fs::Real)
    return Dict("dif" => d2 .- d1,
        "corr" => cor(d1, d2),
        "spect_dif"=> calc_spect_dif(collect(skipmissing(d1)), d2, fs))
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
    readdir_store_sims:
"""
function readdir_store_sims(path_to_experiment)
    elements = []
    for element in readdir(path_to_experiment)
        if element[end-3:end] != ".png"
            push!(elements, element)
        end
    end
    return elements
end


@doc """
    nans_detector:
        This function will scream if any nan or missing is found in the ensemble runs
"""
function nans_detector(; experiment::String="test_default_ens", variable="H_n")
    # Define some local variables
    locdir = pwd() * "/output/" * experiment * "/"

    # Read directory
    elements = readdir_store_sims(locdir)

    # Look for nans and scream their names, times and variables!
    k = 0
    for e in elements
        d = NCDataset(locdir * "/" * e * "/amod.nc")[variable]
        isnan_list, ismissing_list = isnan.(d), ismissing.(d)
        if sum(isnan_list) > 0
            printstyled("NaN found in $e", color=:red)
            k += 1
        end
        if sum(ismissing_list) > 0
            printstyled("missing found in $e", color=:red)
            k += 1
        end
    end
    if k == 0
        printstyled("no NaN or missing found", color=:green)
    end
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
                    fs = 1 / (proxy_data[kp]["time"][2] - proxy_data[kp]["time"][1])
                    new_d2 = interp_2series(d1, d2)
                    stats_kpv[ka] = compute_comp_stats(d1, new_d2, fs)
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
        if i == length(vars2compare)
            ax = Axis(fig[i, 1], xlabel="Time (kyr)", ylabel=vars2compare[i])
        else
            ax = Axis(fig[i, 1], ylabel=vars2compare[i])
        end
        for j in eachindex(elements)
            d2p = amod_data[elements[j]][vars2compare[i]]
            lines!(ax, amod_data[elements[j]]["time"] ./ 1000, d2p, color="grey")
        end
        k = 1
        for p in 1:length(keys(proxy_data))
            prxnm = collect(keys(proxy_data))[p]
            if vars2compare[i] in keys(proxy_data[prxnm])
                scatter!(ax, proxy_data[prxnm]["time"] ./ 1000, proxy_data[prxnm][vars2compare[i]], color=prx_clr[k], markersize=5)
                k += 1
            end
        end
        xlims!(ax, amod_data[elements[1]]["time"][1]/1000, amod_data[elements[1]]["time"][end]/1000)

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
    fig = Figure(resolution=(1600, 900))
    prx_clr = [:navy, :firebrick, :olive, :indigo]
    for i in eachindex(collect(keys(comp_data)))
        # ---- preamble
        proxy_label = collect(keys(comp_data))[i]
        idx = 1
        for l in 1:length(proxy_label)
            if proxy_label[l] == '_'
                idx = l
                break
            end
        end
        var2plot = proxy_label[1:idx-1]

        cor_ar = [comp_data[proxy_label][var2plot][e]["corr"] for e in elements]
        max_cor = findmax(cor_ar)   
        best_cor = comp_data[proxy_label][var2plot][elements[max_cor[2]]]["dif"]

        # ---- time series
        if i == 1
            ax = Axis(fig[i, 1], title="Time series", ylabel=proxy_label)
        elseif i == length(collect(keys(comp_data)))
            ax = Axis(fig[i, 1], ylabel=proxy_label, xlabel="Time (kyr)")
        else
            ax = Axis(fig[i, 1], ylabel=proxy_label)
        end
        for j in eachindex(elements)
            d2p = amod_data[elements[j]][var2plot]
            lines!(ax, amod_data[elements[j]]["time"] ./ 1000, d2p, color="grey")
        end
        lines!(ax, amod_data[elements[max_cor[2]]]["time"] ./ 1000, amod_data[elements[max_cor[2]]][var2plot], color="green")
        lines!(ax, proxy_data[proxy_label]["time"] ./ 1000, proxy_data[proxy_label][var2plot], color="red")

        # ---- differences
        if i == 1
            ax = Axis(fig[i, 2], title="Difference, Simulated - Proxy")
        elseif i == length(collect(keys(comp_data)))
            ax = Axis(fig[i, 2], xlabel="Time (kyr)")
        else
            ax = Axis(fig[i, 2])
        end
        idx = 1
        for l in 1:length(proxy_label)
            if proxy_label[l] == '_'
                idx = l
                break
            end
        end
        for j in eachindex(elements)
            dij = comp_data[proxy_label][var2plot][elements[j]]["dif"]
            lines!(ax, range(new_t1, new_tend, length=length(dij)) ./1000, dij, color="grey")
        end

        # ---- ensemble correlation boxplot
        if i == 1
            ax_cor = Axis(fig[i, 3], title="Correlation")
        else
            ax_cor = Axis(fig[i, 3])
        end
        lines!(ax, range(new_t1, new_tend, length=length(best_cor)) ./ 1000, best_cor, color="green", label=elements[max_cor[2]])
        boxplot!(ax_cor, ones(length(elements)), cor_ar, color="grey")
        scatter!(ax_cor, 1.0, max_cor[1], color="green", markersize=5)

        ax_cor.xticks = ([1.0], [""])
        fig[i, 6] = Legend(fig, ax, framevisible=false, labelsize=15)

        if i != length(collect(keys(comp_data)))
            ax.xticklabelsvisible = false
            ax_cor.xticklabelsvisible = false
        end

        # ---- proxy percentile values vs simulation
        if i == 1
            ax_vs = Axis(fig[i, 4], aspect=AxisAspect(1.0), title="Percentile")
        elseif i == length(collect(keys(comp_data)))
            ax_vs = Axis(fig[i, 4], aspect=AxisAspect(1.0), xlabel="Simulated")
        else
            ax_vs = Axis(fig[i, 4], aspect=AxisAspect(1.0))
        end
        
        d1 = proxy_data[proxy_label][var2plot]
        p1 = [percentile(d1, p) for p in 1:100]
        list_of_mins, list_of_maxs = [], []
        for j in eachindex(elements)
            d2 = amod_data[elements[j]][var2plot]
            new_d2 = interp_2series(d1, d2)
            p2 = [percentile(new_d2, p) for p in 1:100]
            scatter!(ax_vs, p2, p1, color="grey", markersize=2)
            push!(list_of_mins, minimum(p2))
            push!(list_of_maxs, maximum(p2))
        end
        d2_best_cor = interp_2series(d1, amod_data[elements[max_cor[2]]][var2plot])
        best_p2 = [percentile(d2_best_cor, p) for p in 1:100]
        scatter!(ax_vs, best_p2, p1, color="green", markersize=4)

        min_proxy, max_proxy = minimum(p1), maximum(p1)
        min_amod, max_amod = minimum(list_of_mins), maximum(list_of_maxs)
        min_ax, max_ax = min(min_proxy, min_amod), max(max_proxy, max_amod)
        xlims!(ax_vs, min_ax, max_ax)
        ylims!(ax_vs, min_ax, max_ax)

        # ---- difference in spectrum
        if i == 1
            ax_spect = Axis(fig[i, 5], title="Difference in PSD")
        elseif i == length(collect(keys(comp_data)))
            ax_spect = Axis(fig[i, 5], xlabel="Period (kyr)")
        else
            ax_spect = Axis(fig[i, 5])
        end
        
        vlines!(ax_spect, [21, 41, 100], color="red", linestyle=:dash)
        for j in eachindex(elements)
            dij = comp_data[proxy_label][var2plot][elements[j]]["spect_dif"][1]
            dinx = comp_data[proxy_label][var2plot][elements[j]]["spect_dif"][2]
            lines!(ax_spect, dinx ./ 1000, dij, color="grey")
        end
        best_cor_spect = comp_data[proxy_label][var2plot][elements[max_cor[2]]]["spect_dif"]
        lines!(ax_spect, best_cor_spect[2] ./ 1000, best_cor_spect[1], color="green")
    end
    colsize!(fig.layout, 1, Relative(4/16))
    colsize!(fig.layout, 2, Relative(4/16))
    colsize!(fig.layout, 3, Relative(1/16))
    colsize!(fig.layout, 4, Relative(2/16))
    colsize!(fig.layout, 5, Relative(4/16))
    colsize!(fig.layout, 6, Relative(1/16))
    save(locdir * "comp-stats_ens.png", fig)

    # -- Figure 3
    fig = Figure()
    ax_heatmap = Axis(fig[1, 1], aspect=AxisAspect(1.0), title="Correlation", xlabel="Ensemble")
    cor_arr = Array{Real}(undef, length(elements), length(collect(keys(proxy_data))))
    lbls = []
    for i in 1:size(cor_arr)[2]
        proxy_label = collect(keys(proxy_data))[i]
        idx = 1
        for l in 1:length(proxy_label)
            if proxy_label[l] == '_'
                idx = l
                break
            end
        end
        var2plot = proxy_label[1:idx-1]
        k = 1
        for j in elements
            cor_arr[k, i] = comp_data[proxy_label][var2plot][j]["corr"]
            k += 1
        end
        push!(lbls, proxy_label)
    end
    max_corr = findmax(cor_arr, dims=1)
    c = heatmap!(ax_heatmap, cor_arr, colormap=:speed)
    for point in max_corr[2]
        scatter!(ax_heatmap, point[1], point[2], color="red")
    end
    ax_heatmap.yticks = (1:size(cor_arr)[2], lbls)
    ax_heatmap.xticks = (1:Int(size(cor_arr)[1]/10):size(cor_arr)[1], elements[1:Int(size(cor_arr)[1]/10):end])
    Colorbar(fig[1, 2], c, height=Relative(1 / 3), width=20)
    save(locdir * "corr_heatmap.png", fig)  
end

