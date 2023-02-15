# =============================
#     Program: lib_amod.jl
#     Aim: main functions of AMOD
# =============================
@doc """
        ===========================================================
            Model: AMOD (adimensional ice-sheet-sediment model)
                 by Jorge Alvarez-Solas (Fortran, 2017)
        ===========================================================
            Adapted to Julia by Sergio Pérez-Montero (2022)
"""
function amod(now::OrderedDict, par::OrderedDict, ctl::OrderedDict, vart::Vector)
    # Update variables
    now = update_forward(now, par, ctl, vart)
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
            now = calc_Hseddot(now, par, ctl)
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

@doc """
    amod_loop: contains the time loop of the model
"""
function amod_loop(now::OrderedDict, out::OrderedDict, par::OrderedDict, ctl::OrderedDict, vart::Vector, file)
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

        # run AMOD
        now = amod(now, par, ctl, vart)

        # only update output variable at desired frequency
        if mod(now["time"], ctl["dt_out"]) == 0
            out = update_amod_out(out, now)
            if par["active_outout"]
                write(file, "time = " * string(now["time"]) * "\n")
            end
        end
    end
    return out
end

@doc """
    run_amod: main function of AMOD
"""
function run_amod(; experiment::String="test_default", par_file::String="amod_default.jl", par2change=[])
    ## Now, load arguments
    output_path = load_out(amod_path, experiment)
    load_parf(amod_path, output_path, par_file)

    # -- check if some parameters need to be changed
    if par2change != []
        change_namelist(amod_path * "/output/" * experiment, "namelist.jl", par2change)
    end

    # Assign parameters
    # -- assign parameters, CTL, INCOND, PAR, amod_INCOND
    CTL, amod_INCOND, PAR, OUT, out_precc, out_attr = load_defs(output_path * "namelist.jl")

    ## Open out.out
    # if out.out exists remove it
    isfile(output_path * "/out.out") && rm(output_path * "/out.out")
    f = open(output_path * "/out.out", "w")
    if PAR["active_outout"]
        write(f, "**** Starting AMOD " * experiment * " ... ****\n")
    end

    ## Initialize
    NOW = copy(amod_INCOND)
    OUT = update_amod_out(OUT, NOW) # update output

    if PAR["active_outout"]
        write(f, "time = " * string(NOW["time"]) * "\n")
    end

    ## Let's run!
    # -- define updatable (in time) variables
    VARt = ["H", "B"]
    (PAR["active_ice"]) && (VARt = vcat(VARt, ["Hsed", "T_ice"]))
    (PAR["active_climate"]) && (VARt = vcat(VARt, ["T", "co2", "albedo", "ice_time"]))

    # -- run AMOD
    OUT = amod_loop(NOW, OUT, PAR, CTL, VARt, f)

    ## Create output nc file
    # if outfile exists remove it
    isfile(output_path * "/amod.nc") && rm(output_path * "/amod.nc")
    genout_nc(output_path, "amod.nc", OUT, out_precc, out_attr)

    if PAR["active_outout"]
        write(f, "**** AMOD " * experiment * " done! ****" * "\n")
    end
    close(f)
    (PAR["active_outout"] == false) && rm(output_path * "/out.out")

    ## Done!
end

@doc """
    run_amod_ensemble: 
        Takes an (ordered) dictionary of parameters to be exchanged
        in order to run an ensemble and runs AMOD for each combination
"""
function run_amod_ensemble(par2per::OrderedDict; experiment::String="test_default_ens", par_file::String="amod_default.jl")
    # First, obtain simulations
    perm = calc_permutations(par2per)

    # Second, get number of simulations
    nperm = length(perm)

    # Third, create ensemble directory and run each permutation in it
    isdir(amod_path * "/output/" * experiment) || mkdir(amod_path * "/output/" * experiment)
    for i in 1:nperm
        if i < 10   # create subdirectory name
            experimenti = "/s0$i" * "/"
        else
            experimenti = "/s$i" * "/"
        end
        run_amod(experiment=experiment * experimenti, par_file=par_file, par2change=perm[i])
    end

    # Done!
end

@doc """
    run_amod_lhs: 
        Takes a dictionary of parameters and create LHS, key => (min, max)
        to run an ensemble and runs AMOD for each combination (not valid for bool parameters)
"""
function run_amod_lhs(par2per::Dict, nsim::Int; experiment::String="test_default_ens", par_file::String="amod_default.jl")
    parnames = collect(keys(par2per))
    # First, calculate LHS
    permutations, permutations_dict = gen_lhs(par2per, nsim; pars_type=2)

    # Now, create ensemble directory
    if isdir(amod_path * "/output/" * experiment)
        rm(amod_path * "/output/" * experiment, recursive=true)
    end
    mkdir(amod_path * "/output/" * experiment)

    # Plot (if possible) the LHS
    mkdir(amod_path * "/output/" * experiment * "/results")
    if length(par2per) == 2
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2])
        scatter!(ax, permutations)
        save(amod_path * "/output/" * experiment * "/results/lhs.png", fig)
    elseif length(par2per) == 3
        fig = Figure()
        ax = Axis3(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2], zlabel=parnames[3])
        scatter!(ax, permutations)
        save(amod_path * "/output/" * experiment * "/results/lhs.png", fig)
    end

    # Run each permutation
    nperms = length(permutations_dict)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in ProgressBar(1:nperms)
        experimenti = "/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        run_amod(experiment=experiment * experimenti, par_file=par_file, par2change=permutations_dict[i])
    end

    # Done!
end


# Some shortcuts

function run_tests(; test1a=true, test1b=true, test1c=true, test2=true, test3=true)

    # Test 1a. Just ice dynamics and artificial insolation
    if test1a
        # -- only paleo
        par_ice_artif = Dict("time_init" => -2e6, "time_end" => 0,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "ins", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 1.5e-8,
            "pr_ref" => 0.7, "lambda" => 0.07)
        run_amod(experiment="test_ice-artif_def", par_file="amod_default.jl", par2change=par_ice_artif)
        plot_amod(experiment="test_ice-artif_def", vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true)
        plot_wavelet(experiment="test_ice-artif_def", MPT=true, fs=1 / 1000, sigma=π)

        # -- paleo + anth
        par_ice_artif_anth = Dict("time_init" => -2e6, "time_end" => 2e6,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "ins", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 1.5e-8,
            "pr_ref" => 0.7, "lambda" => 0.07)
        run_amod(experiment="test_ice-artif-anth_def", par_file="amod_default.jl", par2change=par_ice_artif_anth)
        plot_amod(experiment="test_ice-artif-anth_def", vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true)
        plot_wavelet(experiment="test_ice-artif-anth_def", MPT=true, fs=1 / 1000, sigma=π)
        println("test 1a completed")
    end

    # Test 1b. Just ice dynamics, artificial insolation and linear Clausius-clapeyron
    if test1b
        # -- only paleo
        par_ice_artif_lin = Dict("time_init" => -2e6, "time_end" => 0,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "linear", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 0.3e-7,
            "ka" => 0.008,
            "lambda" => 0.07, "Acc_ref_n" => 0.4)
        run_amod(experiment="test_ice-artif-lin_def", par_file="amod_default.jl", par2change=par_ice_artif_lin)
        plot_amod(experiment="test_ice-artif-lin_def", vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true)
        plot_wavelet(experiment="test_ice-artif-lin_def", MPT=true, fs=1 / 1000, sigma=π)

        # -- no iso
        par_ice_artif_lin_noiso = Dict("time_init" => -2e6, "time_end" => 0,
            "active_iso" => false, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "linear", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 0.3e-7,
            "ka" => 0.008,
            "lambda" => 0.07, "Acc_ref_n" => 0.4)
        run_amod(experiment="test_ice-artif-lin-noiso_def", par_file="amod_default.jl", par2change=par_ice_artif_lin_noiso)
        plot_amod(experiment="test_ice-artif-lin-noiso_def", vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true)
        plot_wavelet(experiment="test_ice-artif-lin-noiso_def", MPT=true, fs=1 / 1000, sigma=π)

        plot_amod_runs(experiment=["test_ice-artif-lin_def", "test_ice-artif-lin-noiso_def"], vars2plot=["ins_n", "H_n", "Hsed_n"], MPT=true)

        # -- paleo + anth
        par_ice_artif_lin_anth = Dict("time_init" => -2e6, "time_end" => 2e6,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "artificial", "ac_case" => "linear", "sm_case" => "PDD",
            "A_t" => 20.0,
            "f_1" => 0.3e-7,
            "ka" => 0.008,
            "lambda" => 0.07, "Acc_ref_n" => 0.4)
        run_amod(experiment="test_ice-artif-lin-anth_def", par_file="amod_default.jl", par2change=par_ice_artif_lin_anth)
        plot_amod(experiment="test_ice-artif-lin-anth_def", vars2plot=["ins_n", "H_n", "Hsed_n"], plot_MPT=true)
        plot_wavelet(experiment="test_ice-artif-lin-anth_def", MPT=true, fs=1 / 1000, sigma=π)
        println("test 1b completed")
    end

    # Test 1c. Just ice dynamics, laskar insolation and linear Clausius-clapeyron
    if test1c
        # -- only paleo
        par_ice = Dict("time_init" => -2e6, "time_end" => 0,
            "active_iso" => true, "active_sed" => true, "active_climate" => false, "active_ice" => true,
            "ins_case" => "laskar", "ac_case" => "linear", "sm_case" => "PDD",
            "A_t" => 35.0,
            "f_1" => 0.3e-7,
            "ka" => 0.003,
            "lambda" => 0.1, "Acc_ref_n" => 0.4)
        run_amod(experiment="test_ice_def", par_file="amod_default.jl", par2change=par_ice)
        plot_amod(experiment="test_ice_def", vars2plot=["ins_n", "H_n", "B_n", "Hsed_n"], plot_MPT=true)
        plot_wavelet(experiment="test_ice_def", MPT=true, fs=1 / 1000, sigma=π)
        println("test 1c completed")
    end

    # Test 2. Just climate
    if test2
        # -- only paleo
        par_clim = Dict("time_init" => -1e6, "time_end" => 0,
            "height_temp" => "useH",
            "active_iso" => true, "active_sed" => false, "active_climate" => true, "active_ice" => false,
            "ins_case" => "laskar", "ac_case" => "linear", "sm_case" => "ITM", "csi" => 0.11, "cs" => 0.65,
            "lambda" => 0.05, "Acc_ref_n" => 0.1)
        run_amod(experiment="test_clim_useH", par_file="amod_default.jl", par2change=par_clim)
        plot_amod(experiment="test_clim_useH", vars2plot=["ins_n", "H_n", "Z_n", "B_n", "T_n", "co2_n", "V_n"], plot_MPT=true)
        plot_wavelet(experiment="test_clim_useH", MPT=true, fs=1 / 1000, sigma=π)

        par_clim = Dict("time_init" => -1e6, "time_end" => 0,
            "height_temp" => "useZ",
            "active_iso" => true, "active_sed" => false, "active_climate" => true, "active_ice" => false,
            "ins_case" => "laskar", "ac_case" => "linear", "sm_case" => "ITM", "csi" => 0.11, "cs" => 0.65,
            "lambda" => 0.05, "Acc_ref_n" => 0.1)
        run_amod(experiment="test_clim_useZ", par_file="amod_default.jl", par2change=par_clim)
        plot_amod(experiment="test_clim_useZ", vars2plot=["ins_n", "H_n", "Z_n", "B_n", "T_n", "co2_n", "V_n"], plot_MPT=true)
        plot_wavelet(experiment="test_clim_useZ", MPT=true, fs=1 / 1000, sigma=π)

        println("test 2 completed")
    end

end