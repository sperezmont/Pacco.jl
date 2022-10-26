# =============================
#     Program: amod_initialization.jl
#     Aim: functions to set initial conditions
# =============================

function assign_parameters()    # assign the parameters of the namelist
    ctl_assign = RunSettings(time_init, time_end, dt, dt_out)
    inicond_assign = InitialConditions(H_init, Hsed_init, T_init, A_init)
    radpar_assign = RadiativeParameters(ins_min, ins_max, ins_prei, co2_prei, T_0, T_ref)
    orbpar_assign = OrbitalParameters(orbital_mode, P_obl, tau_obl, P_pre, tau_pre, P_exc, tau_exc)
    geopar_assign = GeoParameters(active_iso, B_eq, tau_bed, T_mantle, Q_geo)
    icepar_assign = IceParameters(v_kin, k_1, f_1, f_2, fstream_min, fstream_max, L, C_s, k, c, A_m, A_t, lambda)
    return ctl_assign, inicond_assign, radpar_assign, orbpar_assign, geopar_assign, icepar_assign 
end

function amod_init()    # initialize the model
    now_init = NOW(ctl.time_init,        # time
                inicond.H_init,     # H
                inicond.Hsed_init,  # H_sed
                inicond.T_init,     # T
                inicond.A_init,     # A
                radpar.T_ref,       # T_sl
                0.0,                # TMB                
                0.0,                # S
                0.0,                # B
                0.0,                # M
                0.0,                # Acc
                0.0,                # SMB
                0.0,                # U_d
                0.0,                # U_b
                0.0,                # U
                0.0,                # T_surf
                0.0,                # tau_b
                0.0,                # tau_d
                0.0,                # Q_dif
                0.0,                # Q_difup
                0.0,                # Q_difdown
                0.0,                # Q_drag
                0.0,                # alpha
                0.0,                # Q_adv
                0.0,                # fstream
                0.0,                # fstream_ref
                0.0,                # Hdot
                0.0,                # Hseddot
                0.0,                # Bdot
                0.0,                # Tdot
                0.0,                # fstreamdot
                radpar.ins_prei,    # ins
                radpar.co2_prei)    # co2
    return now_init
end
