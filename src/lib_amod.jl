# =============================
#     Program: lib_amod.jl
#     Aim: main functions of AMOD
# =============================
@doc """
        ===========================================================
            Function: amod.jl
            Model: AMOD (adimensional ice-sheet-sediment model)
                 by Jorge Alvarez-Solas (Fortran, 2017)
        ===========================================================
            Adapted to Julia by Sergio Pérez-Montero (2022)
"""
function amod(now, par, ctl, vart)    
    # Update variables
    now = update_forward(now, par, ctl, vart)
    now = update_Z(now, par)

    # -- calculate reference value for albedo
    now = calc_albedo_ref(now, par)

    # -- thermomechanical coupling (For now, the rate factor A is constant -- jas)

    ## Ice Dynamics
    if par["active_dynamics"]
        now = calc_taud(now, par)           # -- driving stress
        now = calc_taub(now, par)           # -- basal shear stress

        now = calc_Ud(now, par)             # -- mean driving velocity
        now = calc_Ub(now, par)             # -- mean basal velocity

        now = calc_fstream(now, par, ctl)   # -- stream fraction

        now = calc_U(now, par)              # -- total velocity
    end

    # -- surface temperature
    now = calc_T_surf(now, par)

    # -- bedrock temperature (currently prescribed -- jas)

    # -- total pressure
    now = calc_P(now, par)

    # -- mass balance
    now = calc_SMB(now, par)
    now = calc_TMB(now, par)

    ## Ice thermodynamics
    now = calc_Qdif(now, par)       # -- diffusion
    now = calc_Qdrag(now, par)      # -- drag

    # Calculate time evolution
    # -- ice thickness
    now = calc_Hdot(now, par)

    # -- sediment thickness
    (par["active_sed"]) && (now = calc_Hseddot(now, par, ctl))

    # -- bedrock elevation
    now = calc_Bdot(now, par)

    # -- ice temperature
    now = calc_Ticedot(now, par)

    # -- regional temperature
    now = calc_Tdot(now, par)

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
        now = calc_Tsl(now, par)

        # run AMOD
        now = amod(now, par, ctl, vart)

        # only update output variable at desired frequency
        if mod(now["time"], ctl["dt_out"]) == 0
            out = update_amod_out(out, now)
            if par["active_outout"]
                write(file, "time = " * string(now["time"]) * " --> " * "ins = " * string(now["ins"]) * " --> " * "T_sl = " * string(now["T_sl"]) * " --> " * "H = " * string(now["H"]) * "\n")
            end
        end
    end
    return out
end

@doc """
    run_amod: main function of AMOD
"""
function run_amod(out_name="test_default", par_file="amod_default.jl", par2change=[])
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
        write(f, "time = " * string(NOW["time"]) * " --> " * "ins = " * string(NOW["ins"]) * " --> " * "T_sl = " * string(NOW["T_sl"]) * " --> " * "H = " * string(NOW["H"]) * "\n")
    end

    ## Let's run!
    # -- define updatable (in time) variables
    VARt = ["H", "Hsed", "B", "E", "V", "T", "T_ice", "co2", "albedo", "ice_time"]
    
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
    run_amod_ens: 
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
            out_namei = "/s0$i" * "_" * join(perm[i].keys .* string.(perm[i].vals), "-") * "/"
        else
            out_namei = "/s$i" * "_" * join(perm[i].keys .* string.(perm[i].vals), "-") * "/"
        end
        run_amod(out_name * out_namei, par_file, perm[i])
    end

    # Done!
end




