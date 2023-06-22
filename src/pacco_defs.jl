# =============================
#     Program: pacco_defs.jl
#     Aim: definition of states and parameters 
# =============================
# pronostic variables
prognostic = ["T", "pCO2", "iceage", "alpha", "H", "Hsed", "zb", "Tice", "fstr"]
prognostic_longnames = ["Regional air temperature", "Carbon dioxide concentration",
    "Ice age", "System albedo", "Ice thickness", "Sediment layer thickness",
    "Bedrock elevation", "Ice temperature", "Streaming fraction"]

diagnostic = ["I", "Tsl", "Tref",  # radiative forcing and climate response
    "z", "A", "V",  # ice geometry
    "alpha_ref", "Tsurf",   # climate parameters
    "s", "a",   # ice-sheet mass balance
    "taud", "taub", "vd", "vb", "fstr_ref", "v", # ice dynamics
    "Qdif", "Qdrag"]   # ice thermodynamics
diagnostic_longnames = ["Insolation", "Sea-level temperature", "Reference air temperature",
    "Ice-sheet surface elevation", "Ice-sheet area", "Ice-sheet volume",
    "Reference albedo", "Ice surface temperature",
    "Snowfall", "ablation",
    "Driving stress", "Basal stress", "Deformational velocity", "Basal velocity", "Reference streaming fraction", "Total velocity",
    "Diffusive heat", "Dragging heat"]

global states_u = vcat(prognostic, diagnostic)
states_names = vcat(prognostic_longnames, diagnostic_longnames)

global lsu = length(states_u)
global lprog = length(prognostic)
global ldiag = length(diagnostic)

global states_comp = ["Inorm", "Ianom", "m", "RCO2"]

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
    u0 = vcat([p.T0, p.pCO20, p.iceage0, p.alpha0, p.H0, p.Hsed0, p.zb0, p.Tice0, p.fstr0], # prognostic states
        zeros(ldiag)) # diagnostic variables

    calc_diagnostic_variables!(u0, p, p.time_init - p.time_spinup)

    # Output file settings
    out_attr = Dict(
        # -- ctl
        "time" => Dict("units" => "yr", "longame" => "Simulation Time", "group" => "Time"),
        # -- Prognostic variables
        "T" => Dict("units" => "K", "longame" => "Regional Temperature", "group" => "Thermodynamics"),
        "pCO2" => Dict("units" => "ppm", "longame" => "pCO2 Concentration", "group" => "Climate"),
        "iceage" => Dict("units" => "yr", "longame" => "Ice Age", "group" => "Climate"),
        "alpha" => Dict("units" => "--", "longame" => "System albedo", "group" => "Climate"),
        "H" => Dict("units" => "m", "longame" => "Ice Thickness", "group" => "Dynamics"),
        "Hsed" => Dict("units" => "m", "longame" => "Sediment layer Thickness", "group" => "Dynamics"),
        "zb" => Dict("units" => "m", "longame" => "Bedrock Elevation", "group" => "Dynamics"),
        "Tice" => Dict("units" => "K", "longame" => "Ice Temperature", "group" => "Thermodynamics"),
        "fstr" => Dict("units" => "--", "longame" => "Stream Fraction", "group" => "Dynamics"),
        # -- Diagnostic variables
        "I" => Dict("units" => "W/m²", "longame" => "Insolation", "group" => "Forcing"),
        "R" => Dict("units" => "K", "longame" => "Temperature anomaly due to Radiative Forcing", "group" => "Climate"),
        "Tsl" => Dict("units" => "K", "longame" => "Sea-level Temperature", "group" => "Climate"),
        "Tref" => Dict("units" => "K", "longame" => "Reference Air Temperature", "group" => "Climate"),
        "z" => Dict("units" => "m", "longame" => "Ice Surface Elevation", "group" => "Geometry"),
        "A" => Dict("units" => "km²", "longame" => "Ice Extent", "group" => "Geometry"),
        "V" => Dict("units" => "m SLE", "longame" => "Ice Volume", "group" => "Geometry"),
        "alpha_ref" => Dict("units" => "--", "longame" => "System Reference albedo", "group" => "Climate"),
        "Tsurf" => Dict("units" => "K", "longame" => "Surface Temperature", "group" => "Climate"),
        "s" => Dict("units" => "m/a", "longame" => "Accumulation rate", "group" => "Thermodynamics"),
        "a" => Dict("units" => "m/a", "longame" => "Surface Ablation rate", "group" => "Thermodynamics"),
        "taud" => Dict("units" => "Pa", "longame" => "Driving stress", "group" => "Dynamics"),
        "taub" => Dict("units" => "Pa", "longame" => "Basal stress", "group" => "Dynamics"),
        "vd" => Dict("units" => "m/a", "longame" => "Deformational Velocity", "group" => "Dynamics"),
        "vb" => Dict("units" => "m/a", "longame" => "Basal Velocity", "group" => "Dynamics"),
        "fstr_ref" => Dict("units" => "--", "longame" => "Reference Stream Fraction", "group" => "Dynamics"),
        "v" => Dict("units" => "m/a", "longame" => "Total Velocity", "group" => "Dynamics"),
        "Qdif" => Dict("units" => "K/a", "longame" => "Diffusive heat flux", "group" => "Thermodynamics"),
        "Qdrag" => Dict("units" => "K/a", "longame" => "Basal Friction Heating", "group" => "Thermodynamics"),
        # -- Composite variables
        "Inorm" => Dict("units" => "W/m²", "longame" => "Normalized Insolation", "group" => "Forcing"),
        "Ianom" => Dict("units" => "W/m²", "longame" => "Insolation anomaly", "group" => "Forcing"),
        "m" => Dict("units" => "m/a", "longame" => "Mass Balance", "group" => "Thermodynamics"),
        "RCO2" => Dict("units" => "W/m²", "longame" => "pCO2 radiative effect", "group" => "Climate"),
    )
    return u0, out_attr
end


