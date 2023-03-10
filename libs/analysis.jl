# =============================
#     Program: analysis.jl
#     Aim: analyze results from PACCO simulations
#     Author: Sergio Pérez-Montero, 2023.01.26
# =============================

function write_filtered_ensemble(outpath, runs)
    isfile(outpath * "/good_runs.txt") && rm(outpath * "/good_runs.txt")
    f = open(outpath * "/good_runs.txt", "w")
    for run in runs
        write(f, run * "\n")
    end
    close(f)
end

"""
    detect_deglaciation(x, dt)
detects complete deglaciations in a time series

## Arguments
* `x` time series (> 0) to analyze (e.g ice thickness or ice volume)
* `t` time dimension

## Return
`number_of_deglaciations` and `timing`
"""
function detect_deglaciation(x::Vector, t::Vector)
    dx = x[2:end] .- x[1:end-1]
    dt = t[2:end] .- t[1:end-1]

    dxdt = (dx ./ dt)

    is_deglaciation = zeros(length(dxdt))
    for i in eachindex(dxdt)
        if x[i+1] == 0.0    # no ice
            if sign(dxdt[i]) < 0   # glacial termination found
                is_deglaciation[i] = 1.0
            end
        else
            is_deglaciation[i] = 0.0
        end
    end
    number_of_deglaciations = sum(is_deglaciation)
    timing = vcat(0.0, is_deglaciation) # in this way, we extend to the same length as original data
    return number_of_deglaciations, timing
end

function cut_time_series(d1::Dict, d2::Dict)
    elements_d1, elements_d2 = collect(keys(d1)), collect(keys(d2))

    first_ts = [abs(d1[k]["time"][1]) for k in elements_d1]
    proxy_min_t, pacco_min_t = findmin(first_ts), abs(d2[elements_d2[1]]["time"][1])
    min_t = min(proxy_min_t[1], pacco_min_t)
    if min_t == proxy_min_t[1]
        min_name = collect(elements_d1)[proxy_min_t[2]]
        new_t1, new_tend = d1[min_name]["time"][1], d1[min_name]["time"][end]
    else
        min_name = elements[1]
        new_t1, new_tend = d2[elements[1]]["time"][1], d2[elements[1]]["time"][end]
    end

    for k in elements_d1
        if k != min_name
            idx = findmin(abs.(d1[k]["time"] .- new_t1))[2]
            for k2 in keys(d1[k])
                d1[k][k2] = d1[k][k2][idx:end]
            end
        end
    end

    for k in elements_d2
        if k != min_name
            idx = findmin(abs.(d2[k]["time"] .- new_t1))[2]
            for k2 in keys(d2[k])
                d2[k][k2] = d2[k][k2][idx:end]
            end
        end
    end
    return d1, d2, new_t1, new_tend
end


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
    freq = getMeanFreq(ContinuousWavelets.computeWavelets(length(d), wvt)[1], fs)
    S = abs.(S) .^ 2
    S = S ./ sum(S, dims=2)

    return S, freq
end

@doc """
    calc_coi:
        calculates cone of influence from a wavelet analysis
        adapted from Torrence and Compo, 1998, BAMS
        https://es.mathworks.com/help/wavelet/ug/boundary-effects-and-the-cone-of-influence.html
"""
function calc_coi(t, f, cf)
    predtimes = sqrt(2) .* cf ./ f
    tmax, tmin = maximum(t), minimum(t)
    t1 = tmin .+ predtimes
    t2 = tmax .- predtimes
    t_samples = vcat(t1, t2)
    p_samples = vcat(1 ./ f, 1 ./ f)
    return t_samples, p_samples
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
    G1, G2 = G1[freq.>1/ftocut], G2[freq.>1/ftocut]
    freq = freq[freq.>1/ftocut]

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
        "spect_dif" => calc_spect_dif(collect(skipmissing(d1)), d2, fs))
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
    gen_stats_dict:
        generates the dictionary of stats
        dimensions of d --> d[element][variable][time]
"""
function gen_stats_dict(d::Dict)
    stats_dict = Dict(v => Dict() for v in keys(d))
    for i in keys(d)
        stats_i = Dict(v => Dict() for v in keys(d[i]))
        for j in keys(d[i])
            if j != "time"
                stats_ij = compute_stats(d[i][j])
                stats_i[j] = stats_ij
            end
        end
        delete!(stats_i, "time")
        stats_dict[i] = stats_i
    end
    return stats_dict
end

sentence_rule1 = "only considers runs with at least 1 complete deglaciation in ice thickness and with main periodicity between 80-120 kyr"
"""
    rule1(runs)
$(sentence_rule1)
"""
function rule1(runs; sentence=sentence_rule1)
    var_to_load = "H_n" 

    good, bad, stats = [], [], []
    for run in runs
        df = NCDataset(run, "r")
        x, t = df[var_to_load][:], df["time"][:]
        number_of_deglaciations, timing = detect_deglaciation(x, t)
        G, f = calc_spectrum(x, 1 / (t[2] - t[1]))
        G, f = G[f.>1/150e3], f[f.>1/150e3]   # filtering
        G = G ./ sum(G)
        periods = 1 ./ f
        maxG = findmax(G)   # (maxval, index)
        main_period = periods[maxG[2]]

        stats_run = Dict("number_of_deglaciations" => number_of_deglaciations, 
                         "mean" => mean(x), "std" => std(x), "min" => minimum(x), "max" => maximum(x),
                         "main_period" => main_period)
        if number_of_deglaciations > 1
            if stats_run["main_period"] in 80e3..120e3
                if stats_run["mean"] > 1000
                    push!(good, run)
                end
            end
        else
            push!(bad, run)
        end
        push!(stats, stats_run)
        close(df)
    end

    return good, bad, stats, sentence
end


"""
    analyze_runs(;experiment = "test_ensemble_default", experiments = [], rules = 1)
analyzes the run/ensemble `experiment` or the runs in `experiments` based on the rules dictated by `rule`

## Arguments
* `experiment` Name of the run or ensemble to analyze
* `experiments` Vector with the names of the specific runs to analyze
* `rule` rules to apply (check ?rule# to obtain info about the different options)

## Return
* `good_runs` Vector with the names of the runs that fullfill the rules
* `bad_runs` Vector with the names of the runs that do not fullfill the rules
* `stats` Dictionary with the stats associated to `good_runs`
"""
function analyze_runs(; experiment="test_ensemble_default", experiments=[], rule::Integer=1)
    isensemble, data_to_load = is_experiment_or_experiments(experiment, experiments)

    if rule == 1
        good_runs, bad_runs, stats, rule_sentence = rule1(data_to_load)
    else
        error("Rule option not recognized")
    end
        
    if experiments != []
        println("In the experiments selected I have found $(string(length(good_runs))) / $(string(length(data_to_load))) good runs")
    else
        println("In the experiment $(experiment) I have found $(string(length(good_runs))) / $(string(length(data_to_load))) good runs")
    end

    if experiments != []
        write_filtered_ensemble(pwd() * "/output/$(experiments[1])/", good_runs)
    else
        if isensemble
            write_filtered_ensemble(pwd() * "/output/$(experiment)/results/", good_runs)
        else
            write_filtered_ensemble(pwd() * "/output/$(experiment)/", good_runs)
        end
    end
    println("following rule $(string(rule)): $(rule_sentence)")

    return good_runs, bad_runs, stats
end
