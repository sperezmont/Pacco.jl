# Load src modules
include("./amod_dynamics.jl")
include("./amod_thermodynamics.jl")

@doc """
        ===========================================================
            Function: run_amod
            Model: AMOD (adimensional ice-sheet-sediment model)
                 by Jorge Alvarez-Solas (Fortran, 2017)
        ===========================================================
            Adapted to Julia by Sergio Pérez-Montero (2022)
"""
function run_amod(now, par, ctl)
    # Define some local variables and parameteres
    kt_ann = par["k"] * sec_year
    #qgeo_ann = par["Q_geo"] * sec_year * 1e-3

    tau_kin = par["L"] / par["v_kin"]   # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)

    # Update variables
    # -- thicknesses
    now["H"] = max(now["H"] + now["Hdot"] * ctl["dt"], 0.0)
    if active_sed
        now["Hsed"] = min(max(now["Hsed"] + now["Hseddot"] * ctl["dt"], 0.0), 1.0)  # sediments go from 0 to 1
    end


    # -- bedrock
    now["B"] = now["B"] + now["Bdot"] * ctl["dt"]

    # -- ice surface
    if par["active_iso"]
        now["Z"] = max(now["H"] + now["B"], 
                        0.0 + (1 - (rhoi/rhow) * now["H"]))     # Pattyn 2017, Robinson 2020
    else
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
    now["Bdot"] = -(now["B"] - par["B_eq"] + now["H"] * rhoi / rhom) / par["tau_bed"] # needs further improvement -- spm 2022.11.17

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
        now = run_amod(now, par, ctl)

        # only update output variable at desired frequency
        if mod(now["time"], ctl["dt_out"]) == 0
            out = update_amod_out(out, now)
            write(file, "time = " * string(now["time"]) * " --> " * "ins = " * string(now["ins"]) * " --> " * "T_sl = " * string(now["T_sl"]) * " --> " * "H = " * string(now["H"]) * "\n")
        end
    end
    return out
end