# =============================
#     Program: analysis.jl
#     Aim: analyze results from PACCO simulations
#     Author: Sergio Pérez-Montero, 2023.01.26
# =============================

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
sentence_rule0 = "considers all runs"
"""
    rule0(runs)
$(sentence_rule0)
"""
function rule0(runs; sentence=sentence_rule0)
    mask = []
    for run in runs
        push!(mask, true)
    end

    return mask, sentence
end

sentence_rule1 = "only considers runs \n 
                    * with at least 1 complete deglaciation in ice thickness \n "
"""
    rule1(runs)
$(sentence_rule1)
"""
function rule1(runs; sentence=sentence_rule1)
    var_to_load = "H_n"

    mask = []
    for run in runs
        if run[1:3] == "NS_"
            push!(mask, false)
            continue
        end

        df = NCDataset(run, "r")
        x, t = df[var_to_load][:], df["time"][:]
        number_of_deglaciations, timing = detect_deglaciation(x, t)

        if number_of_deglaciations > 1
            push!(mask, true)
        else
            push!(mask, false)
        end
        close(df)
    end

    return mask, sentence
end

sentence_rule2 = "only considers runs \n 
                    * with main periodicities between 80-120 kyr \n"
"""
    rule2(runs)
$(sentence_rule2)
"""
function rule2(runs; sentence=sentence_rule2)
    var_to_load = "H_n"
    mask = []
    for run in runs
        if run[1:3] == "NS_"
            push!(mask, false)
            continue
        end

        df = NCDataset(run, "r")
        x, t = df[var_to_load][:], df["time"][:]
        G, f = calc_spectrum(x, 1 / (t[2] - t[1]))
        time2cut = t[Int(ceil(length(t) / 2))]
        G, f = G[f.>-1/time2cut], f[f.>-1/time2cut]   # filtering
        G = G ./ sum(G)
        periods = 1 ./ f
        maxG = findmax(G)   # (maxval, index)
        maxperiod = periods[maxG[2]]

        #peaks = [idx[1] for idx in findlocalmaxima(G)]  # find main periodicities 
        #peak_is_valid = 0.0
        #for peak in peaks
        #     if periods[peak] in 80e3 .. 120e3
        #         peak_is_valid += 1.0
        #     end
        # end

        # if peak_is_valid > 0
        #     push!(mask, true)
        # else
        #     push!(mask, false)
        # end
        if maxperiod in 80e3 .. 120e3
            push!(mask, true)
        else
            push!(mask, false)
        end

        close(df)
    end

    return mask, sentence
end

sentence_rule3 = "only considers runs with\n 
                    *  mean ice thickness greater than 500 m\n"
"""
    rule3(runs)

$(sentence_rule3)
"""
function rule3(runs; sentence=sentence_rule3)
    var_to_load = "H_n"

    mask = []
    for run in runs
        if run[1:3] == "NS_"
            push!(mask, false)
            continue
        end

        df = NCDataset(run, "r")
        x = df[var_to_load][:]
        mean_x = mean(x)
        if mean_x > 500
            push!(mask, true)
        else
            push!(mask, false)
        end
        close(df)
    end

    return mask, sentence
end

sentence_rule4 = "only considers runs with\n 
                    *  different values of derivative in the two sides of each peak\n"
"""
    rule4(runs)
$(sentence_rule4)
"""
function rule4(runs; sentence=sentence_rule4, window=10)
    printstyled("This rule is work in progress, use carefully\n", color=:red)
    mask = []
    for run in runs
        if run[1:3] == "NS_"
            push!(mask, false)
            continue
        end

        df = NCDataset(run, "r")
        x, dxdt = df["H_n"][:], df["Hdot_n"][:]

        peaks = [idx[1] for idx in findlocalmaxima(x)]  # find glacial maxima 
        number_of_true_glacial_cycles = 0
        for idx in peaks[1:end-1]
            if (peaks[1] > window) && (idx < (length(dxdt) - window))   # takes a window of `window * dt_out` size
                if abs(abs(mean(dxdt[idx-window:idx-1])) - abs(mean(dxdt[idx+1:idx+window]))) < 0.01
                    number_of_true_glacial_cycles += 1
                end
                # if abs(dxdt[idx-window]) < abs(dxdt[idx+window])
                #     number_of_true_glacial_cycles += 1
                # end
            end
        end
        if number_of_true_glacial_cycles > 1
            push!(mask, true)
        else
            push!(mask, false)
        end
        close(df)
    end

    return mask, sentence
end

sentence_rule5 = "only considers runs \n 
                    *  that end with mean ice thickness in a window `window` less than 1% of the average ice thickness\n"
"""
    rule5(runs)

$(sentence_rule5)
"""
function rule5(runs; sentence=sentence_rule5, window=10)
    var_to_load = "H_n"

    mask = []
    for run in runs
        if run[1:3] == "NS_"
            push!(mask, false)
            continue
        end

        df = NCDataset(run, "r")
        x = df[var_to_load][:]
        mean_x = mean(x[end-window:end])

        if mean_x < 0.01 * mean(x)
            push!(mask, true)
        else
            push!(mask, false)
        end
        close(df)
    end

    return mask, sentence
end

sentence_rule6 = "only considers runs with\n 
                    *  mean ice sheet velocity major than 10 m/yr\n"
"""
    rule6(runs)

$(sentence_rule6)
"""
function rule6(runs; sentence=sentence_rule6)
    var_to_load = "U_n"

    mask = []
    for run in runs
        if run[1:3] == "NS_"
            push!(mask, false)
            continue
        end

        df = NCDataset(run, "r")
        x = df[var_to_load][:]
        if mean(x) > 10
            push!(mask, true)
        else
            push!(mask, false)
        end
        close(df)
    end

    return mask, sentence
end

"""
    analyze_runs(;experiment = "test_ensemble_default", experiments = [], rules = 1)
analyzes the run/ensemble `experiment` or the runs in `experiments` based on the rules dictated by `rule`

## Arguments
* `experiment` Name of the run or ensemble to analyze
* `experiments` Vector with the names of the specific runs to analyze
* `rule` rules to apply (check ?rule# to obtain info about the different options)
* `reanalyze` if != [] reads and analyze a the provided pre-analyzed ensemble

## Return
* `good_runs` Vector with the names of the runs that fullfill the rules
* `bad_runs` Vector with the names of the runs that do not fullfill the rules
* `stats` Dictionary with the stats associated to `good_runs`
"""
function analyze_runs(; experiment="test_ensemble_default", experiments=[], rule::Integer=1, reanalyze=[])
    if experiments != []
        path_to_results = get_path_to_results_or_runs(experiments[1], "results")
    else
        path_to_results = get_path_to_results_or_runs(experiment, "results")
    end

    if reanalyze == []
        gr_files = read_dir_good_runs(path_to_results)
        if gr_files != []
            rm.(path_to_results .* "/" .* gr_files)
        end
        isensemble, data_to_load = is_experiment_or_experiments(experiment, experiments)
    else
        isensemble = true
        data_to_load = readlines(path_to_results * "/$(reanalyze)") .* "/pacco.nc"
    end

    if rule === 0
        mask, rule_sentence = rule0(data_to_load)
    elseif rule == 1
        mask, rule_sentence = rule1(data_to_load)
    elseif rule == 2
        mask, rule_sentence = rule2(data_to_load)
    elseif rule == 3
        mask, rule_sentence = rule3(data_to_load)
    elseif rule == 4
        mask, rule_sentence = rule4(data_to_load)
    elseif rule == 5
        mask, rule_sentence = rule5(data_to_load)
    elseif rule == 6
        mask, rule_sentence = rule6(data_to_load)
    else
        error("Rule option not recognized/implemented")
    end

    if experiments != []
        println("In the experiments selected I have found $(sum(mask)) / $(string(length(data_to_load))) good runs")
    else
        println("In the experiment $(experiment) I have found $(sum(mask)) / $(string(length(data_to_load))) good runs")
    end

    write_filtered_ensemble(path_to_results, mask, rule, reanalyze)
    println("following rule $(string(rule)): $(rule_sentence)")
end

"""
    see_rules()
shows documentation of each rule defined
"""
function see_rules()
    println("rule 0: $(sentence_rule0)")
    println("rule 1: $(sentence_rule1)")
    println("rule 2: $(sentence_rule2)")
    println("rule 3: $(sentence_rule3)")
    println("rule 4: $(sentence_rule4)")
    println("rule 5: $(sentence_rule5)")
    println("rule 6: $(sentence_rule6)")
end