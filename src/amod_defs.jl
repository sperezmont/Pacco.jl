# =============================
#     Program: amod_defs.jl
#     Aim: definition of parameters as dictionaries
# =============================

using DataStructures

# We define run parameters as dictionaries and model variables as a vector 
# -- Earth constants
global T_0 = t_0 + degK
global T_ref = t_ref + degK

# -- run control settings
CTL = OrderedDict(
    "time_init" => time_init::Real,
    "time_end" => time_end::Real,
    "dt" => dt::Real,
    "dt_out" => dt_out::Real
)

# -- initial conditions
INCOND = OrderedDict(
    "H_init" => H_init::Real,
    "Hsed_init" => Hsed_init::Real,
    "T_init" => (t_init + degK)::Real,
    "A_init" => A_init::Real
)

# -- run parameters
PAR = OrderedDict(
    "ins_case" => ins_case::String,     # Radiative Parameters
    "ins_day" => ins_day::Real,
    "ins_lat" => ins_lat::Real,
    "ins_min" => ins_min::Real,
    "ins_max" => ins_max::Real,
    "ins_prei" => ins_prei::Real,
    "active_radco2" => active_radco2::Bool,
    "co2_prei" => co2_prei::Real,
    "A_t" => A_t::Real,
    "time_ant" => time_ant::Real,
    "A_ant" => A_ant::Real,
    "tau_ant" => tau_ant::Real,
    "P_obl" => P_obl::Real,              # Orbital Parameters
    "tau_obl" => tau_obl::Real,
    "P_pre" => P_pre::Real,
    "tau_pre" => tau_pre::Real,
    "P_exc" => P_exc::Real,
    "tau_exc" => tau_exc::Real,
    "active_iso" => active_iso::Bool,      # Geophysical Parameters
    "B_eq" => B_eq::Real,
    "tau_bed" => tau_bed::Real,
    "T_mantle" => (t_mantle + degK)::Real,
    "H_mantle" => H_mantle::Real,
    "Q_geo" => Q_geo::Real,
    "v_kin" => v_kin::Real,           # Dynamics
    "active_sed" => active_sed::Bool,
    "f_1" => f_1::Real,
    "f_2" => f_2::Real,
    "fstream_min" => fstream_min::Real,
    "fstream_max" => fstream_max::Real,
    "L" => L::Real,
    "C_s" => C_s::Real,
    "c" => c::Real,
    "ud_case" => ud_case::String,
    "glen_n" => glen_n::Real,
    "ub_case" => ub_case::String,
    "pr_ref" => pr_ref::Real,             # Thermodynamics
    "A_pr" => A_pr::Real,
    "lambda" => lambda::Real,
    "T_sb" => (t_sb + degK)::Real,
    "melt_offset" => (melt_offset + degK)::Real,
    "tsurf_case" => tsurf_case::String,
    "k" => k::Real,
    "cc_case" => cc_case::String,
    "RH" => RH::Real,
    "k_pr" => k_pr::Real,
    "tau_w" => tau_w::Real,
    "T_snow" => (t_snow + degK)::Real,
    "T_rain" => (t_rain + degK)::Real,
    "sm_case" => sm_case::String
)

# -- simulated variables
OUT = OrderedDict(
    "time" => [],            # values of the entire simulation
    "H" => [],
    "Hsed" => [],
    "T" => [],
    "A" => [],
    "T_sl" => [],
    "TMB" => [],
    "SMB" => [],
    "Z" => [],
    "B" => [],
    "M" => [],
    "Acc" => [],
    "U_d" => [],
    "U_b" => [],
    "U" => [],
    "T_surf" => [],
    "tau_b" => [],
    "tau_d" => [],
    "Q_dif" => [],
    "Q_difup" => [],
    "Q_difdown" => [],
    "Q_drag" => [],
    "alpha" => [],
    "Q_adv" => [],
    "fstream" => [],
    "fstream_ref" => [],
    "Hdot" => [],
    "Hseddot" => [],
    "Bdot" => [],
    "Tdot" => [],
    "fstreamdot" => [],
    "ins" => [],
    "ins_norm" => [],
    "co2" => [],
    "P" => [],
    "exc" => [],
    "long_peri" => [],
    "obl" => []
)

# Assign initial conditions
amod_INCOND = OrderedDict(
    "time" => CTL["time_init"],
    "H" => INCOND["H_init"],     # Add 1e-6 just to avoid NaN when calculating Q_dif
    "Hsed" => INCOND["Hsed_init"],
    "T" => INCOND["T_init"],
    "A" => INCOND["A_init"],
    "T_sl" => t_ref + degK,
    "TMB" => 0.0,
    "SMB" => 0.0,
    "Z" => 0.0,
    "B" => 0.0,      # For now it is 0.0, but it should depend on the amount of ice thickness -- spm 2022.11.17
    "M" => 0.0,
    "Acc" => 0.0,
    "U_d" => 0.0,
    "U_b" => 0.0,
    "U" => 0.0,
    "T_surf" => t_ref + degK,     # this should depend on if there is or not ice
    "tau_b" => 0.0,
    "tau_d" => 0.0,
    "Q_dif" => 0.0,
    "Q_difup" => 0.0,
    "Q_difdown" => 0.0,
    "Q_drag" => 0.0,
    "alpha" => 0.0,
    "Q_adv" => 0.0,
    "fstream" => 0.0,
    "fstream_ref" => 0.0,
    "Hdot" => 0.0,
    "Hseddot" => 0.0,
    "Bdot" => 0.0,
    "Tdot" => 0.0,
    "fstreamdot" => 0.0,
    "ins" => PAR["ins_prei"],
    "ins_norm" => 0.0,
    "co2" => PAR["co2_prei"],
    "P" => P_sl,
    "exc" => 0.0,
    "long_peri" => 0.0,
    "obl" => 0.0
)
if PAR["ins_case"] == "artificial" # initialize orbital parameters
    amod_INCOND["long_peri"], amod_INCOND["obl"], amod_INCOND["exc"] = orbital_params(amod_INCOND["time"])
else
    amod_INCOND["long_peri"], amod_INCOND["obl"], amod_INCOND["exc"] = 0, 0, 0
end

# Output file settings
out_precc = Float64
out_attr = OrderedDict(
    "time" => Dict("units" => "yr", "long_name" => "Simulation Time", "group" => "Time"),                    # Time
    "ins" => Dict("units" => "W/m²", "long_name" => "Insolation", "group" => "Radiation"),                        # Radiation
    "ins_norm" => Dict("units" => "W/m²", "long_name" => "Normalized Insolation", "group" => "Radiation"),
    "co2" => Dict("units" => "ppm", "long_name" => "co2 concentration", "group" => "Radiation"),
    "Z" => Dict("units" => "m", "long_name" => "Surface Elevation", "group" => "Geometry"),                   # Geometry
    "B" => Dict("units" => "m", "long_name" => "Bedrock Elevation", "group" => "Geometry"),
    "H" => Dict("units" => "m", "long_name" => "Ice Thickness", "group" => "Geometry"),
    "Hsed" => Dict("units" => "m", "long_name" => "Sediment Thickness", "group" => "Geometry"),
    "U" => Dict("units" => "m/a", "long_name" => "Total Velocity", "group" => "Dynamics"),                       # Dynamics
    "U_d" => Dict("units" => "m/a", "long_name" => "Deformational Velocity", "group" => "Dynamics"),
    "U_b" => Dict("units" => "m/a", "long_name" => "Basal Velocity", "group" => "Dynamics"),
    "tau_d" => Dict("units" => "Pa", "long_name" => "Driving stress", "group" => "Dynamics"),
    "tau_b" => Dict("units" => "Pa", "long_name" => "Basal Stress (friction)", "group" => "Dynamics"),
    "fstream" => Dict("units" => "--", "long_name" => "Stream Fraction", "group" => "Dynamics"),
    "fstream_ref" => Dict("units" => "--", "long_name" => "Stream Fraction", "group" => "Dynamics"),
    "TMB" => Dict("units" => "m/a", "long_name" => "Total Mass Balance", "group" => "Thermodynamics"),            # Thermodynamics
    "SMB" => Dict("units" => "m/a", "long_name" => "Surface Mass Balance", "group" => "Thermodynamics"),
    "Acc" => Dict("units" => "m/a", "long_name" => "Surface Accumulation", "group" => "Thermodynamics"),
    "M" => Dict("units" => "m/a", "long_name" => "Surface Melt", "group" => "Thermodynamics"),
    "T" => Dict("units" => "K", "long_name" => "Ice Temperature", "group" => "Thermodynamics"),
    "Q_dif" => Dict("units" => "K/a", "long_name" => "Diffusive heat", "group" => "Thermodynamics"),
    "Q_difup" => Dict("units" => "K/a", "long_name" => "Atmosphere diffusive heat", "group" => "Thermodynamics"),
    "Q_difdown" => Dict("units" => "K/a", "long_name" => "Basal diffusive heat", "group" => "Thermodynamics"),
    "Q_adv" => Dict("units" => "K/a", "long_name" => "Advective Heat Dissipation", "group" => "Thermodynamics"),
    "Q_drag" => Dict("units" => "K/a", "long_name" => "Basal Friction Heating", "group" => "Thermodynamics"),
    "alpha" => Dict("units" => "--", "long_name" => "Temperate Smoothing Function", "group" => "Thermodynamics"),
    "T_surf" => Dict("units" => "K", "long_name" => "Surface Temperature", "group" => "Thermodynamics"),
    "T_sl" => Dict("units" => "K", "long_name" => "Sea-level Temperature", "group" => "Thermodynamics"),
    "P" => Dict("units" => "Pa", "long_name" => "Air Pressure", "group" => "Thermodynamics"),
    "Hdot" => Dict("units" => "m/a", "long_name" => "Ice Thickness Change", "group" => "Derivatives"),              # Derivatives
    "Hseddot" => Dict("units" => "m/a", "long_name" => "Sediment Thickness Change", "group" => "Derivatives"),
    "Bdot" => Dict("units" => "m/a", "long_name" => "Bedrock Elevation Change", "group" => "Derivatives"),
    "Tdot" => Dict("units" => "K/a", "long_name" => "Ice Temperature Change", "group" => "Derivatives"),
    "fstreamdot" => Dict("units" => "--", "long_name" => "Stream fraction change", "group" => "Derivatives")
)


