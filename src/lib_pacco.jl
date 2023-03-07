# =============================
#     Program: lib_pacco.jl
#     Aim: main functions of PACCO
# =============================

"""
    update_forward(now, par, ctl, vars2update) 

takes vars2update and calculates their time evolution vars2update is a vector with the names of the variables to update vars2update is computed in run_pacco() function 

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters
* `ctl` Run control parameters
* `vars2update` Vector of variables to update in time as var = var + vardot * dt
## Return
updated `now` dictionary
"""
function update_forward(now::OrderedDict, par::OrderedDict, ctl::OrderedDict, vars2update::Vector)
    for hm in par["hemisphere"], v in vars2update
        # -- calculate time evolution
        variab, vardot = v * "_" * hm, v * "dot_" * hm
        now[variab] += now[vardot] * ctl["dt"]    # now = now + dnow/dt * dt

        # -- modify if desired
        if variab in ["H_n", "H_s"]
            now[variab] = max(now[variab], 0.0)
        elseif variab in ["Hsed_n", "Hsed_s"]
            now[variab] = min(max(now[variab], 0.0), 1.0) # sediments go from 0 to 1
        elseif variab in ["B_n", "B_s"]
            (~par["active_iso"]) && (now[variab] = par["B_eq_"*hm]) # reupdate to equilibrium value
        elseif variab in ["T_ice_n", "T_ice_s"]
            now[variab] = min(now[variab], degK)
        elseif variab in ["co2_n", "co2_s"]
            now[variab] = max(now[variab], 1.0)
        end
    end
    return now
end

"""
     pacco(now, par, ctl, vars2update)

PACCO model 

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters
* `ctl` Run control parameters
* `vars2update` Vector of variables to update in time as var = var + vardot * dt
## Return
updated `now` dictionary
"""
function pacco(now::OrderedDict, par::OrderedDict, ctl::OrderedDict, vars2update::Vector)
    # Update variables
    now = update_forward(now, par, ctl, vars2update)
    now = update_Z(now, par)

    now = calc_E(now, par)
    now = calc_V(now, par)

    # -- calculate reference value for albedo
    (par["active_climate"]) && (now = calc_albedo_ref(now, par))

    # -- thermomechanical coupling (For now, the rate factor A is constant -- jas)

    ## Ice Dynamics
    if par["active_ice"]
        now = calc_tau_d(now, par)           # -- driving stress
        now = calc_tau_b(now, par)           # -- basal shear stress

        now = calc_U_d(now, par)             # -- mean driving velocity
        now = calc_U_b(now, par)             # -- mean basal velocity

        now = calc_fstream(now, par, ctl)   # -- stream fraction

        now = calc_U(now, par)              # -- total velocity
    end

    # -- surface temperature
    now = calc_T_surf(now, par)

    # -- bedrock temperature (currently prescribed -- jas)

    # -- mass balance
    now = calc_SMB(now, par)
    now = calc_TMB(now, par)

    ## Ice thermodynamics
    if par["active_ice"]
        now = calc_Qdif(now, par)       # -- diffusion
        now = calc_Qdrag(now, par)      # -- drag
    end

    # Calculate time evolution
    # -- ice thickness
    now = calc_Hdot(now, par)

    # -- sediment thickness
    if par["active_ice"]
        if par["active_sed"]
            now = calc_Hseddot(now, par)
        end
        # -- ice temperature
        now = calc_T_icedot(now, par)
    end

    # -- bedrock elevation
    now = calc_Bdot(now, par)

    # -- regional temperature
    (par["active_climate"]) && (now = calc_Tdot(now, par))

    # -- albedo
    now = calc_albedodot(now, par)

    # -- co2
    now = calc_co2dot(now, par)

    return now
end

"""
    pacco_loop(now, out, par, ctl, vars2update, file)
calculates the time loop of the model

## Attributes
* `now` Dictionary with values of the model variables at current timestep
* `out` Dictionary with stored values of the model variables (output)
* `par` Dictionary with run parameters
* `ctl` Run control parameters
* `vars2update` Vector of variables to update in time as var = var + vardot * dt
* `file` out.out file  
## Return
updated `out` dictionary
"""
function pacco_loop(now::OrderedDict, out::OrderedDict, par::OrderedDict, ctl::OrderedDict, vars2update::Vector, file)
    time_length = ceil((ctl["time_end"] - ctl["time_init"]) / ctl["dt"])
    for n in 1:time_length
        # update simulation time
        now["time"] = ctl["time_init"] + n * ctl["dt"]

        # update contour variables
        if par["active_climate"]
            now = calc_T_rf(now, par)     # -- rf due to coupled climate    
        else
            now = calc_T_sl(now, par)    # -- T_sl from insolation
        end

        # run PACCO
        now = pacco(now, par, ctl, vars2update)

        # only update output variable at desired frequency
        if mod(now["time"], ctl["dt_out"]) == 0
            out = update_pacco_out(now, out)
            if par["active_outout"]
                write(file, "time = " * string(now["time"]) * "\n")
            end
        end

        # Check if run presents instabilities
        break_iter = check_run(now)
        (break_iter) && (break)
    end
    return out
end

"""
    run_pacco(; experiment, par_file, par2change)
        
main function of PACCO (all arguments are optional)

## Attributes 
* `experiment`  Name of the experiment to perform
* `par_file`    Name of the parameter file to use
* `pars2change` Dictionary with the (parameter => new value) to change in the original `par_file` parameter file
## Return
nothing
"""
function run_pacco(; experiment::String="test_default", par_file::String="pacco_default.jl", par2change=[])
    ## Now, load arguments
    output_path = load_out(pacco_path, experiment)
    load_parf(pacco_path, output_path, par_file)

    # -- check if some parameters need to be changed
    if par2change != []
        change_namelist(pacco_path * "/output/" * experiment, "namelist.jl", par2change)
    end

    # Assign parameters
    # -- assign parameters
    ctl, now, par, out, out_prec, out_attr = load_defs(output_path * "namelist.jl")

    ## Open out.out
    # if out.out exists remove it
    isfile(output_path * "/out.out") && rm(output_path * "/out.out")
    f = open(output_path * "/out.out", "w")
    if par["active_outout"]
        write(f, "**** Starting PACCO " * experiment * " ... ****\n")
    end

    ## Initialize
    out = update_pacco_out(now, out) # update output

    if par["active_outout"]
        write(f, "time = " * string(now["time"]) * "\n")
    end

    ## Let's run!
    # -- define updatable (in time) variables
    vars2update = ["H", "B"]
    (par["active_ice"]) && (vars2update = vcat(vars2update, ["Hsed", "T_ice"]))
    (par["active_climate"]) && (vars2update = vcat(vars2update, ["T", "co2", "albedo", "ice_time"]))

    # -- run PACCO
    out = pacco_loop(now, out, par, ctl, vars2update, f)

    ## Create output nc file
    # if outfile exists remove it
    isfile(output_path * "/pacco.nc") && rm(output_path * "/pacco.nc")
    genout_nc(output_path, "pacco.nc", out, out_prec, out_attr)

    if par["active_outout"]
        write(f, "**** PACCO " * experiment * " done! ****" * "\n")
    end
    close(f)
    (par["active_outout"] == false) && rm(output_path * "/out.out")

    ## Done!
end

"""
    run_pacco_ensemble(par2per; experiment, par_file)
        
Takes an (ordered) dictionary of parameters to be exchanged in order to run an ensemble and runs PACCO for each combination
Deprecated, use `run_pacco_lhs()` instead

## Attributes 
* `par2per`     Dictionary with parameters to be exchanged
## Optional attributes
* `experiment`  Name of the experiment to perform
* `par_file`    Name of the parameter file to use
"""
function run_pacco_ensemble(par2per::OrderedDict; experiment::String="test_default_ens", par_file::String="pacco_default.jl")
    # First, obtain simulations
    perm = calc_permutations(par2per)

    # Second, get number of simulations
    nperm = length(perm)

    # Third, create ensemble directory and run each permutation in it
    isdir(pacco_path * "/output/" * experiment) || mkdir(pacco_path * "/output/" * experiment)
    for i in 1:nperm
        if i < 10   # create subdirectory name
            experimenti = "/s0$i" * "/"
        else
            experimenti = "/s$i" * "/"
        end
        run_pacco(experiment=experiment * experimenti, par_file=par_file, par2change=perm[i])
    end

    # Done!
end

"""
    run_pacco_lhs(par2per, nsim; experiment, par_file)
        
Takes a dictionary of parameters and create LHS, key => (min, max) to run an ensemble and runs PACCO for each combination

## Attributes 
* `par2per`     Dictionary with parameters to be exchanged, key => (min, max)
* `nsim`        Number of simulations to perform

## Optional attributes
* `experiment`  Name of the experiment to perform
* `par_file`    Name of the parameter file to use
## Return
nothing
"""
function run_pacco_lhs(par2per::Dict, nsim::Int; experiment::String="test_default_ens", par_file::String="pacco_default.jl")
    parnames = collect(keys(par2per))
    # First, calculate LHS
    permutations, permutations_dict = gen_lhs(par2per, nsim; pars_type=2)

    # Now, create ensemble directory
    if isdir(pacco_path * "/output/" * experiment)
        rm(pacco_path * "/output/" * experiment, recursive=true)
    end
    mkdir(pacco_path * "/output/" * experiment)

    # Plot (if possible) the LHS
    mkdir(pacco_path * "/output/" * experiment * "/results")
    if length(par2per) == 2
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2])
        scatter!(ax, permutations)
        save(pacco_path * "/output/" * experiment * "/results/lhs.png", fig)
    elseif length(par2per) == 3
        fig = Figure()
        ax = Axis3(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2], zlabel=parnames[3])
        scatter!(ax, permutations)
        save(pacco_path * "/output/" * experiment * "/results/lhs.png", fig)
    end

    # Run each permutation
    nperms = length(permutations_dict)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in ProgressBar(1:nperms)
        experimenti = "/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        run_pacco(experiment=experiment * experimenti, par_file=par_file, par2change=permutations_dict[i])
    end

    # Done!
end


# Some shortcuts
"""
    run_tests(; test1a, test1b, test1c, test2, test3)
        
Shortcut to run some general tests (dev), all attributes are optional

## Attributes (default value is `false`)
* `test1a`  Just ice dynamics and artificial insolation
* `test1b`  Just ice dynamics, artificial insolation and linear Clausius-clapeyron
* `test1c`  Just ice dynamics, laskar insolation and linear Clausius-clapeyron
* `test2`   Just climate
* `test3`   Default mode, ice dynamics + coupled climate
## Return
nothing
"""
function run_tests(; test1a=false, test1b=false, test1c=false, test2=false, test3=false)

    # Test 1a. Just ice dynamics and artificial insolation
    if test1a
        expname = "test1a_ice-artif_def"
        pars = Dict("time_init" => -2e6, "time_end" => 2e6,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "ins", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 1.5e-8,
            "pr_ref" => 0.7, "lambda" => 0.07)
        run_pacco(experiment=expname, par_file="pacco_default.jl", par2change=pars)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        println("test 1a $(expname) completed")
    end

    # Test 1b. Just ice dynamics, artificial insolation and linear Clausius-clapeyron
    if test1b
        # -- only paleo
        pars = Dict("time_init" => -2e6, "time_end" => 2e6,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "linear", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 0.3e-7,
            "ka" => 0.008,
            "lambda" => 0.07, "Acc_ref_n" => 0.4)
        expname = "test1b_ice-artif-lin_def"
        run_pacco(experiment=expname, par_file="pacco_default.jl", par2change=pars)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)

        # -- no iso
        pars["active_iso"] = false
        expname2 = "test1b_ice-artif-lin-noiso_def"
        run_pacco(experiment=expname2, par_file="pacco_default.jl", par2change=pars)
        plot_pacco(experiment=expname2, vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname2, MPT=true, fs=1 / 1000, sigma=π)

        plot_pacco(experiments=[expname, expname2], vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true)

        println("test 1b $(expname) completed")
    end

    # Test 1c. Just ice dynamics, laskar insolation and linear Clausius-clapeyron
    if test1c
        # -- only paleo
        expname = "test1c_ice_def"
        pars = Dict("time_init" => -2e6, "time_end" => 2e6,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "laskar", "ac_case" => "linear", "sm_case" => "PDD",
            "A_t" => 35.0,
            "f_1" => 0.3e-7,
            "ka" => 0.003,
            "lambda" => 0.1, "Acc_ref_n" => 0.4)
        run_pacco(experiment=expname, par_file="pacco_default.jl", par2change=pars)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "B_n", "Hsed_n"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        println("test 1c $(expname) completed \n Note: not fully calibrated")
    end

    # Test 2. Just climate
    if test2
        expname = "test2_clim_def"
        pars = Dict("time_init" => -5e5, "time_end" => 1e6,
            "height_temp" => "useZ",
            "active_iso" => true, "active_sed" => false, "active_climate" => true, "active_ice" => false,
            "ins_case" => "laskar", "ac_case" => "linear", "sm_case" => "ITM",
            "csi" => 0.07, "cs" => 0.65, "csz" => 0.0065, "cco2" => 2.0,
            "ka" => 0.008, "ki" => 0.0095, "ktco2" => 7.0, "melt_offset" => -5.0,
            "lambda" => 0.01, "Acc_ref_n" => 0.1)
        run_pacco(experiment=expname, par_file="pacco_default.jl", par2change=pars)
        plot_pacco(experiment=expname, vars2plot=["ins_anom_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        println("test 2 $(expname) completed \n Note: needs further improvement in calibration")
    end

    if test3 # default mode, ice + coupled climate 
        expname = "test3_ice-clim_def"
        pars = Dict("time_init" => -8e5, "time_end" => 2e5,
            "height_temp" => "useZ",
            "active_iso" => true, "active_sed" => false, "active_climate" => true, "active_ice" => true,
            "ins_case" => "laskar", "ac_case" => "linear", "sm_case" => "ITM",
            "csi" => 0.08, "cs" => 0.65, "csz" => 0.0065, "cco2" => 2.0,
            "ka" => 0.008, "ki" => 0.008, "ktco2" => 7.0, "t_threshold" => -5.0,
            "lambda" => 0.1, "Acc_ref_n" => 0.4,
            "f_1" => 1e-6, "C_s" => 1e-7,
            "fstream_min_n" => 0.4, "fstream_max_n" => 0.4)
        vrs2plt = ["ins_anom_n", "H_n", "T_ref_n", "M_n", "SMB_n", "U_n", "co2_n"]
        
        pars["Hsed_init_n"] = 0.0
        run_pacco(experiment=expname, par_file="pacco_default.jl", par2change=pars)
        plot_pacco(experiment=expname, vars2plot=vrs2plt, plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=3.45)

        println("test 3 $(expname) completed \n Note: Work in progress")
    end
end