# =============================
#     Program: pacco_defs.jl
#     Aim: definition of states and parameters 
# =============================
# pronostic variables
prognostic = ["T", "C", "iceage", "albedo", "H", "Hsed", "B", "Tice", "fstr"]
prognostic_longnames = ["Regional air temperature", "Carbon dioxide concentration",
    "Ice age", "System albedo", "Ice thickness", "Sediment layer thickness",
    "Bedrock elevation", "Ice temperature", "Streaming fraction"]

diagnostic = ["I", "Tsl", "Tref", "Cref",  # radiative forcing and climate response
    "z", "Surf", "Vol",  # ice geometry
    "albedo_ref", "Tsurf",   # climate parameters
    "s", "a",   # ice-sheet mass balance
    "taud", "taub", "vd", "vb", "fstr_ref", "beta", "v", # ice dynamics
    "hcond", "hvadv", "hdrag", "hgeo",   # ice thermodynamics
    "L", "Hb", "Pe", "albedo_eff", "precipitation"]    # misc   
diagnostic_longnames = ["Insolation", "Sea-level temperature", "Reference air temperature", "Reference carbon dioxide",
    "Ice-sheet surface elevation", "Ice-sheet area", "Ice-sheet volume",
    "Reference albedo", "Ice surface temperature",
    "Snowfall", "ablation",
    "Driving stress", "Basal stress", "Deformational velocity", "Basal velocity", "Reference streaming fraction", "Basal velocity coefficient", "Total velocity",
    "Conductive heat flux", "Vertical advection heat flux", "Drag heat flux", "Geothermal heat flux",
    "Ice-sheet aspect ratio", "Temperate layer thickness", "Peclet number", "Effective albedo", "precipitation"]

global T_idx, C_idx, iceage_idx, albedo_idx, H_idx, Hsed_idx, B_idx, Tice_idx, fstr_idx = 1, 2, 3, 4, 5, 6, 7, 8, 9
global I_idx, Tsl_idx, Tref_idx, Cref_idx = 10, 11, 12, 13
global z_idx, Surf_idx, Vol_idx = 14, 15, 16
global albedo_ref_idx, Tsurf_idx = 17, 18
global s_idx, a_idx = 19, 20
global taud_idx, taub_idx, vd_idx, vb_idx, fstr_ref_idx, beta_idx, v_idx = 21, 22, 23, 24, 25, 26, 27
global hcond_idx, hvadv_idx, hdrag_idx, hgeo_idx = 28, 29, 30, 31
global L_idx, Hb_idx, Pe_idx = 32, 33, 34
global albedo_eff_idx, precipitation_idx = 35, 36

global states_u = vcat(prognostic, diagnostic)
states_names = vcat(prognostic_longnames, diagnostic_longnames)

global lsu = length(states_u)
global lprog = length(prognostic)
global ldiag = length(diagnostic)

global states_comp = ["Inorm", "Ianom", "m", "q", "sw", "lw", "RCO2"]

"""
    display_pacco_variables()
Display the variables, index ...
    
"""
function display_pacco_states()
    printstyled("All PACCO variables are stored in the vector `u = {prognostic, diagnostic}`:\n", color=:light_blue)
    for i in eachindex(states_u)
        println("   $(i)        ---->   $(states_u[i]), $(states_names[i])")
    end

    printstyled("Model/run parameters are stored in the struct `p = Params()`\n", color=:light_blue)
    p = Params()
    nms = fieldnames(typeof(p))
    println("$(nms)")
end

"""
    load_defs
Include parameter file selected in Julia and assign values to model variables

## Arguments
* `p` Parameters struct

## Return
* `u0` Vector with the initial states of the model
* `out_attr` Attributes of each output variable (netCDF file)
"""
function load_defs(p)
    # Check some parameters
    (p.time_init > p.time_end) && (error("time_init is greater than time_end"))

    # Assign initial conditions
    u0 = vcat([p.T0, p.C0, p.iceage0, p.albedo0, p.H0, p.Hsed0, p.B0, p.Tice0, p.fstr0], # prognostic states
        zeros(ldiag)) # diagnostic variables

    calc_diagnostic_variables!(u0, p, p.time_init - p.time_spinup)

    # Check units for Insolation/Energy (I)
    if p.insol_case in ["ISI", "caloric"]
        I_units, I_name = "J/m²", "Energy"
        Inorm_units, Inorm_name = "J/m²", "Normalized Energy"
        Ianom_units, Ianom_name = "J/m²", "Energy anomaly"
    elseif p.insol_case == "input"
        insol_file = split(p.insol_input, "/")[3]
        if (insol_file[1:3] == "ISI") || (insol_file[1:8] == "mean_ISI") || (insol_file[1:7] == "caloric") || (insol_file[1:12] == "mean_caloric")
            I_units, I_name = "J/m²", "Energy"
            Inorm_units, Inorm_name = "J/m²", "Normalized Energy"
            Ianom_units, Ianom_name = "J/m²", "Energy anomaly"
        else
            I_units, I_name = "W/m²", "Insolation"
            Inorm_units, Inorm_name = "W/m²", "Normalized Insolation"
            Ianom_units, Ianom_name = "W/m²", "Insolation anomaly"
        end
    else
        I_units, I_name = "W/m²", "Insolation"
        Inorm_units, Inorm_name = "W/m²", "Normalized Insolation"
        Ianom_units, Ianom_name = "W/m²", "Insolation anomaly"
    end

    # Output file settings
    out_attr = Dict(
        # -- ctl
        "time" => Dict("units" => "yr", "longame" => "Simulation Time", "group" => "Time"),
        # -- Prognostic variables
        "T" => Dict("units" => "K", "longame" => "Regional Temperature", "group" => "Thermodynamics"),
        "C" => Dict("units" => "ppm", "longame" => "C Concentration", "group" => "Climate"),
        "iceage" => Dict("units" => "yr", "longame" => "Ice Age", "group" => "Climate"),
        "albedo" => Dict("units" => "--", "longame" => "System albedo", "group" => "Climate"),
        "H" => Dict("units" => "m", "longame" => "Ice Thickness", "group" => "Dynamics"),
        "Hsed" => Dict("units" => "m", "longame" => "Sediment layer Thickness", "group" => "Dynamics"),
        "B" => Dict("units" => "m", "longame" => "Bedrock Elevation", "group" => "Dynamics"),
        "Tice" => Dict("units" => "K", "longame" => "Ice Temperature", "group" => "Thermodynamics"),
        "fstr" => Dict("units" => "--", "longame" => "Stream Fraction", "group" => "Dynamics"),
        # -- Diagnostic variables
        "I" => Dict("units" => I_units, "longame" => I_name, "group" => "Forcing"),
        "R" => Dict("units" => "K", "longame" => "Temperature anomaly due to Radiative Forcing", "group" => "Climate"),
        "Tsl" => Dict("units" => "K", "longame" => "Sea-level Temperature", "group" => "Climate"),
        "Tref" => Dict("units" => "K", "longame" => "Reference Air Temperature", "group" => "Climate"),
        "Cref" => Dict("units" => "ppm", "longame" => "Reference Carbon Dioxide", "group" => "Climate"),
        "z" => Dict("units" => "m", "longame" => "Ice Surface Elevation", "group" => "Geometry"),
        "Surf" => Dict("units" => "km²", "longame" => "Ice Extent", "group" => "Geometry"),
        "Vol" => Dict("units" => "m SLE", "longame" => "Ice Volume", "group" => "Geometry"),
        "albedo_ref" => Dict("units" => "--", "longame" => "System Reference albedo", "group" => "Climate"),
        "Tsurf" => Dict("units" => "K", "longame" => "Surface Temperature", "group" => "Climate"),
        "s" => Dict("units" => "m/yr", "longame" => "Accumulation rate", "group" => "Thermodynamics"),
        "a" => Dict("units" => "m/yr", "longame" => "Surface Ablation rate", "group" => "Thermodynamics"),
        "taud" => Dict("units" => "Pa", "longame" => "Driving stress", "group" => "Dynamics"),
        "taub" => Dict("units" => "Pa", "longame" => "Basal stress", "group" => "Dynamics"),
        "vd" => Dict("units" => "m/yr", "longame" => "Deformational Velocity", "group" => "Dynamics"),
        "vb" => Dict("units" => "m/yr", "longame" => "Basal Velocity", "group" => "Dynamics"),
        "fstr_ref" => Dict("units" => "--", "longame" => "Reference Stream Fraction", "group" => "Dynamics"),
        "beta" => Dict("units" => "--", "longame" => "Basal velocity coefficient", "group" => "Dynamics"),
        "v" => Dict("units" => "m/yr", "longame" => "Total Velocity", "group" => "Dynamics"),
        "hcond" => Dict("units" => "W/m²", "longame" => "Conductive heat flux", "group" => "Thermodynamics"),
        "hvadv" => Dict("units" => "W/m²", "longame" => "Vertical Heat advection flux", "group" => "Thermodynamics"),
        "hdrag" => Dict("units" => "W/m²", "longame" => "Drag Heat flux", "group" => "Thermodynamics"),
        "hgeo" => Dict("units" => "W/m²", "longame" => "Geothermal Heat flux", "group" => "Thermodynamics"),
        "L" => Dict("units" => "m", "longame" => "Ice-sheet aspect ratio", "group" => "Geometry"),
        "Hb" => Dict("units" => "m", "longame" => "Temperate layer thickness", "group" => "Geometry"),
        "Pe" => Dict("units" => "", "longame" => "Peclet number", "group" => "Thermodynamics"),
        "albedo_eff" => Dict("units" => "", "longame" => "Effective albedo", "group" => "Thermodynamics"),
        "precipitation" => Dict("units" => "m/yr", "longame" => "Precipitation", "group" => "Thermodynamics"),
        # -- Composite variables
        "Inorm" => Dict("units" => Inorm_units, "longame" => Inorm_name, "group" => "Forcing"),
        "Ianom" => Dict("units" => Ianom_units, "longame" => Ianom_name, "group" => "Forcing"),
        "m" => Dict("units" => "m/yr", "longame" => "Mass Balance", "group" => "Thermodynamics"),
        "q" => Dict("units" => "m/yr", "longame" => "Ice flux divergence", "group" => "Thermodynamics"),
        "sw" => Dict("units" => "m/yr", "longame" => "Shortwave contribution to ablation", "group" => "Thermodynamics"),
        "lw" => Dict("units" => "m/yr", "longame" => "Longwave contribution to ablation", "group" => "Thermodynamics"),
        "RCO2" => Dict("units" => "W/m²", "longame" => "C radiative effect", "group" => "Climate"),
    )
    return u0, out_attr
end


