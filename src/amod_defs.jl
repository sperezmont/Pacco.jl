# =============================
#     Program: amod_defs.jl
#     Aim: definition of parameters as dictionaries
# =============================
function load_defs(par_path)
    # Include parameters file
    include(par_path)

    # We define run parameters as dictionaries and model variables as a vector 
    # -- run control settings
    CTL = OrderedDict(
        "time_init" => time_init::Real,
        "time_end" => time_end::Real,
        "dt" => dt::Real,
        "dt_out" => dt_out::Real,
    )

    # -- initial conditions
    INCOND = OrderedDict(
        "H_init_n" => H_init_n,
        "H_init_s" => H_init_s,
        "Hsed_init_n" => Hsed_init_n,
        "Hsed_init_s" => Hsed_init_s,
        "B_init_n" => B_init_n,
        "B_init_s" => B_init_s,
        "T_init_n" => t_init_n + degK,
        "T_init_s" => t_init_s + degK,
        "T_ice_init_n" => t_ice_init_n + degK,
        "T_ice_init_s" => t_ice_init_s + degK,
        "A_init" => A_init
    )

    # -- run parameters
    PAR = OrderedDict(
        # -- dev par (this should be removed for official release)
        "height_temp" => height_temp,
        # -- ctl
        "hemisphere" => hemisphere,
        # -- Switches
        "active_outout" => active_outout,
        "active_iso" => active_iso,
        "active_sed" => active_sed,
        "active_climate" => active_climate,
        "active_ice" => active_ice,
        # -- Cases
        "ins_case" => ins_case,
        "ud_case" => ud_case,
        "ub_case" => ub_case,
        "tsurf_case" => tsurf_case,
        "ac_case" => ac_case,
        "sm_case" => sm_case,
        # -- Orbital forcing parameters 
        "P_obl" => P_obl,
        "tau_obl" => tau_obl,
        "P_pre" => P_pre,
        "tau_pre" => tau_pre,
        "P_exc" => P_exc,
        "tau_exc" => tau_exc,
        "ins_day" => ins_day,
        "ins_lat_n" => ins_lat_n, "ins_lat_s" => ins_lat_s,
        # -- Radiative forcing parameters 
        "ins_min" => ins_min,
        "ins_max" => ins_max,
        "ins_ref_n" => ins_ref_n, "ins_ref_s" => ins_ref_s,
        "ins_prei" => ins_prei,
        "co2_prei" => co2_prei,
        "co2_ref" => co2_ref,
        "tau_co2" => tau_co2,
        "cco2" => cco2,
        "ktco2" => ktco2,
        "A_t" => A_t,
        "time_anth" => time_anth,
        "co2_anth" => co2_anth,
        "At_anth" => At_anth,
        "Ac_anth" => Ac_anth,
        "tau_anth" => tau_anth,
        "albedo_land" => albedo_land,
        "albedo_oldice" => albedo_oldice,
        "albedo_newice" => albedo_newice,
        "albedo_slope" => albedo_slope,
        "albedo_quad" => albedo_quad,
        "tau_albedo" => tau_albedo,
        "csi" => csi,
        "cs" => cs,
        "csz" => csz,
        "T_ref_n" => t_ref_n + degK, "T_ref_s" => t_ref_s + degK,
        "tau_rf_n" => tau_rf_n, "tau_rf_s" => tau_rf_s,
        # -- Geophysical parameters
        "B_eq_n" => B_eq_n, "B_eq_s" => B_eq_s,
        "tau_bed_n" => tau_bed_n, "tau_bed_s" => tau_bed_s,
        "T_mantle_n" => t_mantle_n + degK, "T_mantle_s" => t_mantle_s + degK,
        "H_mantle_n" => H_mantle_n, "H_mantle_s" => H_mantle_s,
        "Q_geo" => Q_geo,
        # -- Geometry
        "L" => L,
        "E_ref_n" => E_ref_n,
        "E_ref_s" => E_ref_s,
        # -- Dynamics
        "v_kin" => v_kin,
        "f_1" => f_1,
        "f_2" => f_2,
        "fstream_min_n" => fstream_min_n, "fstream_min_s" => fstream_min_s,
        "fstream_max_n" => fstream_max_n, "fstream_max_s" => fstream_max_s,
        "C_s" => C_s,
        "glen_n" => glen_n,
        # -- Thermodynamics
        "T_sb" => t_sb + degK,
        "kt" => kt,
        "pr_ref" => pr_ref,
        "A_pr" => A_pr,
        "T_snow" => t_snow + degK,
        "T_rain" => t_rain + degK,
        "lambda" => lambda,
        "melt_offset" => melt_offset + degK,
        "c" => c,
        "km" => km,
        "ki" => ki,
        "ka" => ka,
        "Acc_ref_n" => Acc_ref_n, "Acc_ref_s" => Acc_ref_s,
        "A_te_n" => A_te_n, "A_te_s" => A_te_s
    )

    # Assign initial conditions (first NOW step)
    amod_INCOND = OrderedDict(
        # -- ctl
        "time" => CTL["time_init"],
        # -- forcing
        "ins_n" => PAR["ins_prei"], "ins_s" => PAR["ins_prei"],
        "ins_norm_n" => 0.0, "ins_norm_s" => 0.0,
        "ins_anom_n" => 0.0, "ins_anom_s" => 0.0,
        "T_rf_n" => 0.0, "T_rf_s" => 0.0,
        "T_sl_n" => PAR["T_ref_n"], "T_sl_s" => PAR["T_ref_s"],      # -- amod variables
        # ---- time-updatable
        "H_n" => INCOND["H_init_n"], "H_s" => INCOND["H_init_s"],
        "B_n" => INCOND["B_init_n"], "B_s" => INCOND["B_init_s"],
        "Hsed_n" => INCOND["Hsed_init_n"], "Hsed_s" => INCOND["Hsed_init_s"],
        "E_n" => 0.0, "E_s" => 0.0,
        "V_n" => 0.0, "V_s" => 0.0,
        "T_ice_n" => INCOND["T_ice_init_n"], "T_ice_s" => INCOND["T_ice_init_s"],
        "T_n" => INCOND["T_init_n"], "T_s" => INCOND["T_init_s"],
        "co2_n" => PAR["co2_prei"], "co2_s" => PAR["co2_prei"],
        "albedo_n" => PAR["albedo_land"], "albedo_s" => PAR["albedo_land"],
        "ice_time_n" => 0.0, "ice_time_s" => 0.0,
        "Z_n" => 0.0, "Z_s" => 0.0,
        # ---- albedo
        "albedo_ref_n" => PAR["albedo_land"], "albedo_ref_s" => PAR["albedo_land"],        # ---- ice dynamics
        "tau_d_n" => 0.0, "tau_d_s" => 0.0,
        "tau_b_n" => 0.0, "tau_b_s" => 0.0,
        "U_d_n" => 0.0, "U_d_s" => 0.0,
        "U_b_n" => 0.0, "U_b_s" => 0.0,
        "U_n" => 0.0, "U_s" => 0.0,
        "alpha_n" => 0.0, "alpha_s" => 0.0,
        "fstream_ref_n" => fstream_min_n, "fstream_ref_s" => fstream_min_s,
        "fstreamdot_n" => 0.0, "fstreamdot_s" => 0.0,
        "fstream_n" => fstream_min_n, "fstream_s" => fstream_min_s,
        # ---- thermodynamics
        "A" => INCOND["A_init"],
        "T_surf_n" => degK, "T_surf_s" => degK,
        "Acc_n" => 0.0, "Acc_s" => 0.0,
        "M_n" => 0.0, "M_s" => 0.0,
        "SMB_n" => 0.0, "SMB_s" => 0.0,
        "TMB_n" => 0.0, "TMB_s" => 0.0,
        "Q_difup_n" => 0.0, "Q_difup_s" => 0.0,
        "Q_difdown_n" => 0.0, "Q_difdown_s" => 0.0,
        "Q_dif_n" => 0.0, "Q_dif_s" => 0.0,
        "Q_drag_n" => 0.0, "Q_drag_s" => 0.0,
        # ---- time derivatives
        "Hdot_n" => 0.0, "Hdot_s" => 0.0,
        "Hseddot_n" => 0.0, "Hseddot_s" => 0.0,
        "T_icedot_n" => 0.0, "T_icedot_s" => 0.0,
        "Bdot_n" => 0.0, "Bdot_s" => 0.0,
        "Tdot_n" => 0.0, "Tdot_s" => 0.0,
        "albedodot_n" => 0.0, "albedodot_s" => 0.0,
        "ice_timedot_n" => 1.0, "ice_timedot_n" => 1.0, # this value is not updated, "dumb" variable (ice_time += dt) -- spm 2023.01.05
        "co2dot_n" => 0.0, "co2dot_s" => 0.0
    )

    # -- model variables
    OUT = OrderedDict(
        # -- ctl
        "time" => [],
        # -- forcing
        "ins_n" => [], "ins_s" => [],
        "ins_norm_n" => [], "ins_norm_s" => [],
        "ins_anom_n" => [], "ins_anom_s" => [],
        "T_rf_n" => [], "T_rf_s" => [],
        "T_sl_n" => [], "T_sl_s" => [],
        # -- amod variables
        # ---- time-updatable
        "H_n" => [], "H_s" => [],
        "B_n" => [], "B_s" => [],
        "Hsed_n" => [], "Hsed_s" => [],
        "E_n" => [], "E_s" => [],
        "V_n" => [], "V_s" => [],
        "T_ice_n" => [], "T_ice_s" => [],
        "T_n" => [], "T_s" => [],
        "co2_n" => [], "co2_s" => [],
        "albedo_ref_n" => [], "albedo_ref_s" => [],
        "albedo_n" => [], "albedo_s" => [],
        "ice_time_n" => [], "ice_time_s" => [],
        "Z_n" => [], "Z_s" => [],
        # ---- ice dynamics
        "tau_d_n" => [], "tau_d_s" => [],
        "U_d_n" => [], "U_d_s" => [],
        "U_b_n" => [], "U_b_s" => [],
        "U_n" => [], "U_s" => [],
        "fstream_n" => [], "fstream_s" => [],
        # ---- thermodynamics
        "T_surf_n" => [], "T_surf_s" => [],
        "Acc_n" => [], "Acc_s" => [],
        "M_n" => [], "M_s" => [],
        "SMB_n" => [], "SMB_s" => [],
        "Q_dif_n" => [], "Q_dif_s" => [],
        "Q_drag_n" => [], "Q_drag_s" => [],
        # ---- time derivatives
        "Hdot_n" => [], "Hdot_s" => [],
        "Hseddot_n" => [], "Hseddot_s" => [],
        "T_icedot_n" => [], "T_icedot_s" => [],
        "Bdot_n" => [], "Bdot_s" => [],
        "Tdot_n" => [], "Tdot_s" => [],
        "albedodot_n" => [], "albedodot_s" => [],
        "co2dot_n" => [], "co2dot_s" => []
    )

    # Output file settings
    out_precc = Float64
    out_attr = OrderedDict(
        # -- ctl
        "time" => Dict("units" => "yr", "long_name" => "Simulation Time", "group" => "Time"),
        # -- forcing
        "ins_n" => Dict("units" => "W/m²", "long_name" => "Insolation", "group" => "Forcing"),
        "ins_s" => Dict("units" => "W/m²", "long_name" => "Insolation", "group" => "Forcing"),
        "ins_norm_n" => Dict("units" => "W/m²", "long_name" => "Normalized Insolation", "group" => "Forcing"),
        "ins_norm_s" => Dict("units" => "W/m²", "long_name" => "Normalized Insolation", "group" => "Forcing"),
        "ins_anom_n" => Dict("units" => "W/m²", "long_name" => "Insolation Anomaly", "group" => "Forcing"),
        "ins_anom_s" => Dict("units" => "W/m²", "long_name" => "Insolation Anomaly", "group" => "Forcing"),
        "T_rf_n" => Dict("units" => "K/yr", "long_name" => "Radiative Forcing Temp.", "group" => "Forcing"),
        "T_rf_s" => Dict("units" => "K/yr", "long_name" => "Radiative Forcing Temp.", "group" => "Forcing"),
        "T_sl_n" => Dict("units" => "K", "long_name" => "Sea-level Temperature", "group" => "Forcing"),
        "T_sl_s" => Dict("units" => "K", "long_name" => "Sea-level Temperature", "group" => "Forcing"),
        # -- amod variables
        # ---- time-updatable
        "H_n" => Dict("units" => "m", "long_name" => "Ice Thickness", "group" => "Geometry"),
        "H_s" => Dict("units" => "m", "long_name" => "Ice Thickness", "group" => "Geometry"),
        "B_n" => Dict("units" => "m", "long_name" => "Bedrock Elevation", "group" => "Geometry"),
        "B_s" => Dict("units" => "m", "long_name" => "Bedrock Elevation", "group" => "Geometry"),
        "Hsed_n" => Dict("units" => "m", "long_name" => "Sediment Thickness", "group" => "Geometry"),
        "Hsed_s" => Dict("units" => "m", "long_name" => "Sediment Thickness", "group" => "Geometry"),
        "E_n" => Dict("units" => "km²", "long_name" => "Ice Extent", "group" => "Geometry"),
        "E_s" => Dict("units" => "km²", "long_name" => "Ice Extent", "group" => "Geometry"),
        "V_n" => Dict("units" => "m SLE", "long_name" => "Ice Volume", "group" => "Geometry"),
        "V_s" => Dict("units" => "m SLE", "long_name" => "Ice Volume", "group" => "Geometry"),
        "T_ice_n" => Dict("units" => "K", "long_name" => "Ice Temperature", "group" => "Thermodynamics"),
        "T_ice_s" => Dict("units" => "K", "long_name" => "Ice Temperature", "group" => "Thermodynamics"),
        "T_n" => Dict("units" => "K", "long_name" => "Regional Temperature", "group" => "Thermodynamics"),
        "T_s" => Dict("units" => "K", "long_name" => "Regional Temperature", "group" => "Thermodynamics"),
        "co2_n" => Dict("units" => "ppm", "long_name" => "co2 Concentration", "group" => "Radiative"),
        "co2_s" => Dict("units" => "ppm", "long_name" => "co2 Concentration", "group" => "Radiative"),
        "albedo_n" => Dict("units" => "--", "long_name" => "System Albedo", "group" => "Radiative"),
        "albedo_s" => Dict("units" => "--", "long_name" => "System Albedo", "group" => "Radiative"),
        "ice_time_n" => Dict("units" => "yr", "long_name" => "Ice Age", "group" => "Radiative"),
        "ice_time_s" => Dict("units" => "yr", "long_name" => "Ice Age", "group" => "Radiative"),
        "Z_n" => Dict("units" => "m", "long_name" => "Ice Surface Elevation", "group" => "Geometry"),
        "Z_s" => Dict("units" => "m", "long_name" => "Ice Surface Elevation", "group" => "Geometry"),
        # ---- albedo
        "albedo_ref_n" => Dict("units" => "--", "long_name" => "System Reference Albedo", "group" => "Radiative"),
        "albedo_ref_s" => Dict("units" => "--", "long_name" => "System Reference Albedo", "group" => "Radiative"),
        # ---- ice dynamics
        "tau_d_n" => Dict("units" => "Pa", "long_name" => "Driving stress", "group" => "Dynamics"),
        "tau_d_s" => Dict("units" => "Pa", "long_name" => "Driving stress", "group" => "Dynamics"),
        "tau_b_n" => Dict("units" => "Pa", "long_name" => "Basal stress", "group" => "Dynamics"),
        "tau_b_s" => Dict("units" => "Pa", "long_name" => "Basal stress", "group" => "Dynamics"),
        "U_d_n" => Dict("units" => "m/a", "long_name" => "Driving Velocity", "group" => "Dynamics"),
        "U_d_s" => Dict("units" => "m/a", "long_name" => "Driving Velocity", "group" => "Dynamics"),
        "U_b_n" => Dict("units" => "m/a", "long_name" => "Basal Velocity", "group" => "Dynamics"),
        "U_b_s" => Dict("units" => "m/a", "long_name" => "Basal Velocity", "group" => "Dynamics"),
        "U_n" => Dict("units" => "m/a", "long_name" => "Total Velocity", "group" => "Dynamics"),
        "U_s" => Dict("units" => "m/a", "long_name" => "Total Velocity", "group" => "Dynamics"),
        "alpha_n" => Dict("units" => "--", "long_name" => "Temperate Smoothing Function", "group" => "Thermodynamics"),
        "alpha_s" => Dict("units" => "--", "long_name" => "Temperate Smoothing Function", "group" => "Thermodynamics"),
        "fstream_ref_n" => Dict("units" => "--", "long_name" => "Reference Stream Fraction", "group" => "Dynamics"),
        "fstream_ref_s" => Dict("units" => "--", "long_name" => "Reference Stream Fraction", "group" => "Dynamics"),
        "fstreamdot_n" => Dict("units" => "--", "long_name" => "Stream Fraction Change", "group" => "Derivatives"),
        "fstreamdot_s" => Dict("units" => "--", "long_name" => "Stream Fraction Change", "group" => "Derivatives"),
        "fstream_n" => Dict("units" => "--", "long_name" => "Stream Fraction", "group" => "Dynamics"),
        "fstream_s" => Dict("units" => "--", "long_name" => "Stream Fraction", "group" => "Dynamics"),
        # ---- thermodynamics
        "T_surf_n" => Dict("units" => "K", "long_name" => "Surface Temperature", "group" => "Thermodynamics"),
        "T_surf_s" => Dict("units" => "K", "long_name" => "Surface Temperature", "group" => "Thermodynamics"),
        "Acc_n" => Dict("units" => "m/a", "long_name" => "Accumulation rate", "group" => "Thermodynamics"),
        "Acc_s" => Dict("units" => "m/a", "long_name" => "Accumulation rate", "group" => "Thermodynamics"),
        "M_n" => Dict("units" => "m/a", "long_name" => "Surface Melting rate", "group" => "Thermodynamics"),
        "M_s" => Dict("units" => "m/a", "long_name" => "Surface Melting rate", "group" => "Thermodynamics"),
        "SMB_n" => Dict("units" => "m/a", "long_name" => "Surface Mass Balance", "group" => "Thermodynamics"),
        "SMB_s" => Dict("units" => "m/a", "long_name" => "Surface Mass Balance", "group" => "Thermodynamics"),
        "TMB_n" => Dict("units" => "m/a", "long_name" => "Total Mass Balance", "group" => "Thermodynamics"),
        "TMB_s" => Dict("units" => "m/a", "long_name" => "Total Mass Balance", "group" => "Thermodynamics"),
        "Q_difup_n" => Dict("units" => "K/a", "long_name" => "Diffusive heat flux (Atmosphere)", "group" => "Thermodynamics"),
        "Q_difup_s" => Dict("units" => "K/a", "long_name" => "Diffusive heat flux (Atmosphere)", "group" => "Thermodynamics"),
        "Q_difdown_n" => Dict("units" => "K/a", "long_name" => "Diffusive heat flux (at the base)", "group" => "Thermodynamics"),
        "Q_difdown_s" => Dict("units" => "K/a", "long_name" => "Diffusive heat flux (at the base)", "group" => "Thermodynamics"),
        "Q_dif_n" => Dict("units" => "K/a", "long_name" => "Diffusive heat flux", "group" => "Thermodynamics"),
        "Q_dif_s" => Dict("units" => "K/a", "long_name" => "Diffusive heat flux", "group" => "Thermodynamics"),
        "Q_drag_n" => Dict("units" => "K/a", "long_name" => "Basal Friction Heating", "group" => "Thermodynamics"),
        "Q_drag_s" => Dict("units" => "K/a", "long_name" => "Basal Friction Heating", "group" => "Thermodynamics"),
        # ---- time derivatives
        "Hdot_n" => Dict("units" => "m/a", "long_name" => "Ice Thickness Change", "group" => "Derivatives"),
        "Hdot_s" => Dict("units" => "m/a", "long_name" => "Ice Thickness Change", "group" => "Derivatives"),
        "Hseddot_n" => Dict("units" => "m/a", "long_name" => "Sediments Thickness Change", "group" => "Derivatives"),
        "Hseddot_s" => Dict("units" => "m/a", "long_name" => "Sediments Thickness Change", "group" => "Derivatives"),
        "T_icedot_n" => Dict("units" => "K/a", "long_name" => "Ice Temperature Change", "group" => "Derivatives"),
        "T_icedot_s" => Dict("units" => "K/a", "long_name" => "Ice Temperature Change", "group" => "Derivatives"),
        "Bdot_n" => Dict("units" => "m/a", "long_name" => "Bedrock Elevation Change", "group" => "Derivatives"),
        "Bdot_s" => Dict("units" => "m/a", "long_name" => "Bedrock Elevation Change", "group" => "Derivatives"),
        "Tdot_n" => Dict("units" => "K/a", "long_name" => "Regional Temperature Change", "group" => "Derivatives"),
        "Tdot_s" => Dict("units" => "K/a", "long_name" => "Regional Temperature Change", "group" => "Derivatives"),
        "albedodot_n" => Dict("units" => "1/a", "long_name" => "Albedo Change", "group" => "Derivatives"),
        "albedodot_s" => Dict("units" => "1/a", "long_name" => "Albedo Change", "group" => "Derivatives"),
        "co2dot_n" => Dict("units" => "ppm/a", "long_name" => "co2 Change", "group" => "Derivatives"),
        "co2dot_s" => Dict("units" => "ppm/a", "long_name" => "co2 Change", "group" => "Derivatives")
    )
    return CTL, amod_INCOND, PAR, OUT, out_precc, out_attr
end


