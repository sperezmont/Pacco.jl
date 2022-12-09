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
function amod(now, par, ctl)
    # Define some local variables and parameteres
    kt_ann = par["kt"] * sec_year
    #qgeo_ann = par["Q_geo"] * sec_year * 1e-3

    tau_kin = par["L"] / par["v_kin"]   # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)

    # Update variables
    # -- thicknesses
    now["H"] = max(now["H"] + now["Hdot"] * ctl["dt"], 0.0)
    if active_sed
        now["Hsed"] = min(max(now["Hsed"] + now["Hseddot"] * ctl["dt"], 0.0), 1.0)  # sediments go from 0 to 1
    end

    if par["active_iso"]
        # -- bedrock
        now["B"] = now["B"] + now["Bdot"] * ctl["dt"]

        # -- ice surface
        now["Z"] = max(now["H"] + now["B"],
            0.0 + (1 - (rhoi / rhow) * now["H"]))     # Pattyn 2017, Robinson 2020
    else
        # -- bedrock
        now["B"] = par["B_eq"]

        # -- ice surface
        now["Z"] = now["H"] + par["B_eq"]
    end

    # -- ice temperature
    now["T"] = min(now["T"] + now["Tdot"] * ctl["dt"], degK)       # we do not allow ice above 0ºC

    # -- thermomechanical coupling (For now, the rate factor A is constant -- jas)

    ## Ice Dynamics
    # -- driving stress
    now["tau_d"] = calc_taud(now, par)

    # -- basal shear stress
    now["tau_b"] = calc_taub(now, par)

    # -- velocities
    now["U_d"] = calc_Ud(now, par)
    now["U_b"] = calc_Ub(now, par)

    # -- stream fraction
    now = calc_fstreamdot(now, par, tau_kin)
    now["fstream"] = max(now["fstream"] + now["fstreamdot"] * ctl["dt"], 0.0)

    # -- total velocity
    now["U"] = now["U_d"] + now["fstream"] * now["U_b"]

    # -- surface temperature
    now["T_surf"] = calc_T_surf(now, par)

    # -- bedrock temperature (currently prescribed -- jas)

    # -- total pressure
    now["P"] = P_sl * exp((-now["Z"] * g) / (Rd * now["T_surf"]))   # http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-1/

    # -- surface mass balance
    now = calc_SMB(now, par)
    now["TMB"] = now["SMB"] + 0.0   # for now, TMB = SMB -- 2022.11.17 spm

    ## Ice thermodynamics
    if now["H"] < 10.0  # -- check if there is no ice
        now["Q_dif"] = 0.0
        now["Q_drag"] = 0.0
        now["T"] = now["T_sl"]
    else
        # -- update total diffusion
        now = calc_Qdif(now, par, kt_ann)

        # -- update basal drag heat
        now = calc_Qdrag(now, par)

        # -- update advective heat ??
    end

    # Calculate time evolution
    # -- ice thickness
    now["Hdot"] = now["TMB"] - now["U"] * now["H"] / par["L"]

    # -- sediment thickness
    now["Hseddot"] = -par["f_1"] * now["U"] + par["f_2"] * now["M"] / ctl["dt"]

    # -- bedrock elevation
    if par["active_iso"]
        now["Bdot"] = -(now["B"] - par["B_eq"] + now["H"] * rhoi / rhom) / par["tau_bed"] # needs further improvement -- spm 2022.11.17
    else
        now["Bdot"] = 0.0
    end

    # -- temperature
    now["Tdot"] = now["Q_dif"] + now["Q_drag"]

    return now
end

@doc """
    amod_loop: contains the time loop of the model
"""
function amod_loop(now, out, par, ctl, file)
    time_length = ceil((ctl["time_end"] - ctl["time_init"]) / ctl["dt"])
    for n in 1:time_length
        # update simulation time
        now["time"] = ctl["time_init"] + n * ctl["dt"]

        # update contour variables
        now = calc_Tsl(now, par)

        # run AMOD
        now = amod(now, par, ctl)

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
    OUT = amod_loop(NOW, OUT, PAR, CTL, f)

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

# data = OrderedDict("a"=>[1000,2000,3000], "b"=>[1,2], "c"=>[10,20,30,50])
# data_values, data_keys = collect(values(data)), collect(keys(data)) 
# nperm = 24 
# maxlen, minlen, nvars = 4, 2, 3

# perm = copy(data_values[1])
# for i in 2:nvars
#     perm = Iterators.product(perm, data_values[i])
# end
# perm = collect(perm)



