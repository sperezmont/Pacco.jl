# Load src modules
include("./amod_dynamics.jl")
include("./amod_thermodynamics.jl")

function run_amod(ctl, par, now)
    # =============================
    #     Module: run_amod
    #     Model: AMOD (adimensional ice-sheet-sediment model)
    #            by Jorge Alvarez-Solas (Fortran, 2017)
    # =============================
    #     Adapted to Julia by Sergio Pérez-Montero

    # Define some local variables and parameteres
    kt_ann = par["k"] * sec_year
    qgeo_ann = par["Q_geo"] * sec_year * 1e-3

    tau_kin = par["L"] / par["v_kin"]   # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)

    # Update variables
    # -- thicknesses
    now["H"] = max(now["H"] + now["Hdot"] * ctl["dt"], 0.0)
    now["Hsed"] = max(now["Hsed"] + now["Hseddot"] * ctl["dt"], 0.0)

    # -- bedrock
    now["B"] = now["B"] + now["Bdot"] * ctl["dt"]

    # -- surface
    if par["active_iso"]
        now["Z"] = now["H"] + now["B"]
    else
        now["Z"] = now["H"] + par["B_eq"]
    end

    # -- ice temperature
    now["T"] = now["T"] + now["Tdot"] * ctl["dt"]
    now["T"] = min(now["T"], 0.0)       # we do not allow ice above 0ºC

    # -- thermomechanical coupling (For now, the rate factor A is constant -- jas)

    # -- driving stress
    now["tau_d"] = calc_taud(now, par)

    # -- basal shear stress
    now["tau_b"] = calc_taub(now, par)

    # -- velocities
    now["U_d"] = calc_Ud(now, par)
    now["U_b"] = calc_Ub(now, par)

    # -- stream fraction
    now["fstreamdot"] = calc_fstreamdot(now, par, tau_kin)
    now["fstream"] = max(now["fstream"] + now["fstreamdot"] * ctl["dt"], 0.0)

    # -- total velocity
    now["U"] = now["U_d"] + now["fstream"] * now["U_b"]

    # -- surface temperature
    now["T_surf"] = calc_T_surf(now, par)

    # -- bedrock temperature (currently prescribed -- jas)

    # -- total pressure
    now["P"] = P_sl * exp((-now["Z"] * g) / (Rd * now["T_surf"]))   # http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-1/

    # -- surface mass balance
    now["SMB"] = calc_SMB(now, par)

    # Now, compute thermodynamics
    # -- check if there is no ice
    if now["H"] < 10.0
        now["Q_dif"] = 0.0
        now["Q_drag"] = 0.0
        now["T"] = now["T_sl"]
    end

    # -- update total diffusion
    calc_Qdif(now, PAR, kt_ann)

    # -- update basal drag heat
    now["Q_drag"] = now["fstream"] * now["tau_b"] * now["U_b"] / (par["c"] * rho)

    # -- update advective heat ??

    # Calculate time evolution
    # -- ice thickness
    now["Hdot"] = now["SMB"] - now["U"] * now["H"] / par["L"]

    # -- sediment thickness
    now["Hseddot"] = par["f_1"] * now["U"] + par["f_2"] * now["M"] / ctl["dt"]

    # -- bedrock elevation
    now["Bdot"] = -(now["B"] - par["B_eq"] + now["H"] / 3) / par["tau_bed"] # -- ajr: improve 1/3 to use densities, etc...

    # -- temperature
    now["Tdot"] = now["Q_dif"] + now["Q_drag"]

    return
end