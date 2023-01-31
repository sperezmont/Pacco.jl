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
function amod(now, par, ctl, vart)
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
function amod_loop(now, out, par, ctl, vart, file)
    time_length = ceil((ctl["time_end"] - ctl["time_init"]) / ctl["dt"])
    for n in 1:time_length
        # update simulation time
        now["time"] = ctl["time_init"] + n * ctl["dt"]

        # update contour variables
        if par["active_climate"]
            now = calc_rf(now, par)     # -- rf due to coupled climate    
        else
            now = calc_Tsl(now, par)    # -- T_sl from insolation
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
function run_amod(; out_name="test_default", par_file="amod_default.jl", par2change=[])
    ## Now, load arguments
    output_path = load_out(amod_path, out_name)
    load_parf(amod_path, output_path, par_file)

    # -- check if some parameters need to be changed
    if par2change != []
        change_namelist(amod_path * "/output/" * out_name, "namelist.jl", par2change)
    end

    # Assign parameters
    # -- assign parameters, CTL, INCOND, PAR, amod_INCOND
    CTL, amod_INCOND, PAR, OUT, out_precc, out_attr = load_defs(output_path * "namelist.jl")

    ## Open out.out
    # if out.out exists remove it
    isfile(output_path * "/out.out") && rm(output_path * "/out.out")
    f = open(output_path * "/out.out", "w")
    if PAR["active_outout"]
        write(f, "**** Starting AMOD " * out_name * " ... ****\n")
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
        write(f, "**** AMOD " * out_name * " done! ****" * "\n")
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
function run_amod_ensemble(par2per::OrderedDict; out_name="test_default_ens", par_file="amod_default.jl")
    # First, obtain simulations
    perm = calc_permutations(par2per)

    # Second, get number of simulations
    nperm = length(perm)

    # Third, create ensemble directory and run each permutation in it
    isdir(amod_path * "/output/" * out_name) || mkdir(amod_path * "/output/" * out_name)
    for i in 1:nperm
        if i < 10   # create subdirectory name
            out_namei = "/s0$i" * "/"
        else
            out_namei = "/s$i" * "/"
        end
        run_amod(out_name=out_name * out_namei, par_file=par_file, par2change=perm[i])
    end

    # Done!
end

@doc """
    run_amod_lhs: 
        Takes a dictionary of parameters and create LHS, key => (min, max)
        to run an ensemble and runs AMOD for each combination (not valid for bool parameters)
"""
function run_amod_lhs(par2per::Dict, nsim::Int; out_name="test_default_ens", par_file="amod_default.jl")
    parnames = collect(keys(par2per))
    # First, calculate LHS
    permutations, permutations_dict = gen_lhs(par2per, nsim; pars_type=2)

    # Now, create ensemble directory
    isdir(amod_path * "/output/" * out_name) || mkdir(amod_path * "/output/" * out_name)

    # Plot (if possible) the LHS
    if length(par2per) == 2
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2])
        scatter!(ax, permutations)
        save(amod_path * "/output/" * out_name * "/lhs.png", fig)
    elseif length(par2per) == 3
        fig = Figure()
        ax = Axis3(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2], zlabel=parnames[3])
        scatter!(ax, permutations)
        save(amod_path * "/output/" * out_name * "/lhs.png", fig)
    end

    # Run each permutation
    nperms = length(permutations_dict)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in 1:nperms
        out_namei = "/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        run_amod(out_name=out_name * out_namei, par_file=par_file, par2change=permutations_dict[i])
    end

    # Done!
end

