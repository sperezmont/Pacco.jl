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

    # Update thicknesses
    now["H"] = max(now["H"] + now["Hdot"] * ctl["dt"], 0.0)
    now["Hsed"] = max(now["Hsed"] + now["Hseddot"] * ctl["dt"], 0.0)

    # Update bedrock
    now["B"] = now["B"] + now["Bdot"] * ctl["dt"]

    # Update surface
    if par["active_iso"]
        now["S"] = now["H"] + now["B"]
    else
        now["S"] = now["H"] + par["B_eq"]
    end

    # Update ice temperature
    now["T"] = now["T"] + now["Tdot"] * ctl["dt"]
    now["T"] = min(now["T"], 0.0)       # we do not allow ice above 0ºC

    # Thermomechanical coupling (For now, the rate factor A is constant -- jas)

    # Update driving stress: (rho*g*H*dS/dx)
    now["tau_d"] = calc_taud(now, par)

    # Update basal shear stress
    now["tau_b"] = calc_taub(now, par)

    # Update velocities
    now["U_d"] = calc_Ud(now, par)
    now["U_b"] = calc_Ub(now, par)

    # Update stream fraction
    now["fstream"] = calc_fstream(ctl, now, par, tau_kin)

    # Update total velocity
    now["U"] = now["U_d"] + now["fstream"] * now["U_b"]

    # Update surface temperature
    now["T_surf"] = calc_T_surf(now, par)

    # Update bedrock temperature (currently prescribed -- jas)

    # Update surface mass balance
    now["Acc"] = calc_Acc(now, par)

end