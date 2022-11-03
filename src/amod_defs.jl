# =============================
#     Program: amod_defs.jl
#     Aim: definition of parameters as dictionaries
# =============================

using DataStructures

# We define run parameters as dictionaries and model variables as a vector 
# -- run control settings
global CTL = OrderedDict("time_init"  => time_init::Real,      
                    "time_end"  => time_end::Real,       
                    "dt"        => dt::Real,             
                    "dt_out"    => dt_out::Real         
                    )

# -- initial conditions
global INCOND = OrderedDict("H_init"       => H_init::Real,        
               "Hsed_init"    => Hsed_init::Real,     
               "T_init"       => T_init::Real,        
               "A_init"       => A_init::Real    
                )

# -- run parameters
global PAR = OrderedDict( "ins_min"     => ins_min::Real,         # Radiative Parameters
            "ins_max"     => ins_max::Real,
            "ins_prei"    => ins_prei::Real,
            "co2_prei"    => co2_prei::Real,
            "ins_max"     => ins_max::Real,
            "T_0"         => T_0::Real,
            "T_ref"       => T_ref::Real,
            "orb_case"=> orb_case::String,  # Orbital Parameters
            "P_obl"       => P_obl::Real,
            "tau_obl"     => tau_obl::Real,
            "P_pre"       => P_pre::Real,
            "tau_pre"     => tau_pre::Real,
            "P_exc"       => P_exc::Real,
            "tau_exc"     => tau_exc::Real,
            "active_iso"  => active_iso::Bool,      # Geophysical Parameters
            "B_eq"        => B_eq::Real,
            "tau_bed"     => tau_bed::Real,
            "T_mantle"    => T_mantle::Real,
            "Q_geo"       => Q_geo::Real,
            "v_kin"       => v_kin::Real,           # Dynamics
            "k_1"         => k_1::Real,
            "f_1"         => f_1::Real,
            "f_2"         => f_2::Real,
            "fstream_min" => fstream_min::Real,
            "fstream_max" => fstream_max::Real,
            "L"           => L::Real,
            "C_s"         => C_s::Real,
            "c"           => c::Real,
            "vel_case"    => vel_case::String,
            "glen_n"      => glen_n::Real,           
            "ub_case"    => ub_case::String,    
            "A_m"         => A_m::Real,             # Thermodynamics
            "A_t"         => A_t::Real,
            "lambda"      => lambda::Real,
            "T_sb"        => T_sb::Real,
            "melt_offset" => melt_offset::Real,
            "tsurf_case"  => tsurf_case::String,
            "k"           => k::Real,
            "cc_case"     => cc_case::String
            )

# -- simulated variables
global OUT = OrderedDict( "time"      => [],            # values of the entire simulation
                "H"         => [],
                "Hsed"      => [],
                "T"         => [],
                "A"         => [],    
                "T_sl"      => [], 
                "TMB"       => [],       
                "S"         => [],       
                "B"         => [],           
                "M"         => [],    
                "Acc"       => [],   
                "SMB"       => [],
                "U_d"       => [],       
                "U_b"       => [],  
                "U"         => [],      
                "T_surf"    => [],
                "tau_b"     => [],
                "tau_d"     => [],     
                "Q_dif"     => [],  
                "Q_difup"   => [],
                "Q_difdown" => [],      
                "Q_drag"    => [],   
                "alpha"     => [],  
                "Q_adv"     => [],               
                "fstream"   => [],
                "fstream_ref"=> [],
                "Hdot"      => [],   
                "Hseddot"   => [],  
                "Bdot"      => [],  
                "Tdot"      => [],   
                "fstreamdot"=> [],
                "ins"       => [],
                "co2"       => []
            )

# Assign initial conditions
global amod_INCOND = OrderedDict(
    "time"    => CTL["time_init"],        
    "H"       => INCOND["H_init"],
    "Hsed"    => INCOND["Hsed_init"],
    "T"       => INCOND["T_init"],
    "A"       => INCOND["A_init"],
    "T_sl"    => PAR["T_ref"],
    "TMB"     => 0.0,
    "S"       => 0.0,
    "TMB"     => 0.0,
    "S"       => 0.0,                
    "B"       => 0.0,
    "M"       => 0.0,
    "Acc"     => 0.0,
    "SMB"     => 0.0,
    "U_d"     => 0.0,
    "U_b"     => 0.0,
    "U"       => 0.0,
    "T_surf"  => 0.0,
    "tau_b"   => 0.0,
    "tau_d"   => 0.0,
    "Q_dif"   => 0.0,
    "Q_difup" => 0.0,
    "Q_difdown" => 0.0,
    "Q_drag"  => 0.0,
    "alpha"   => 0.0,
    "Q_adv"   => 0.0,
    "fstream" => 0.0,
    "fstream_ref" => 0.0,
    "Hdot"    => 0.0,
    "Hseddot" => 0.0,
    "Bdot"    => 0.0,
    "Tdot"    => 0.0,
    "fstreamdot" => 0.0,
    "ins"     => PAR["ins_prei"],
    "co2"     => PAR["co2_prei"]
    )

# Output file settings
out_precc = Float64
out_attr = OrderedDict("time"=>Dict("units"=>"yr", "long_name"=>"Simulation Time", "group"=>"Time"),                    # Time
                "ins"=>Dict("units"=>"W/m²", "long_name"=>"Insolation", "group"=>"Radiation"),                        # Radiation
                "co2"=>Dict("units"=>"ppm", "long_name"=>"co2 concentration", "group"=>"Radiation"),
                "S"=>Dict("units"=>"m", "long_name"=>"Surface Elevation", "group"=>"Geometry"),                   # Geometry
                "B"=>Dict("units"=>"m", "long_name"=>"Bedrock Elevation", "group"=>"Geometry"),
                "H"=>Dict("units"=>"m", "long_name"=>"Ice Thickness", "group"=>"Geometry"),
                "Hsed"=>Dict("units"=>"m", "long_name"=>"Sediment Thickness", "group"=>"Geometry"),
                "U"=>Dict("units"=>"m/a", "long_name"=>"Total Velocity", "group"=>"Dynamics"),                       # Dynamics
                "U_d"=>Dict("units"=>"m/a", "long_name"=>"Deformational Velocity", "group"=>"Dynamics"),
                "U_b"=>Dict("units"=>"m/a", "long_name"=>"Basal Velocity", "group"=>"Dynamics"),
                "tau_d"=>Dict("units"=>"Pa", "long_name"=>"Driving stress", "group"=>"Dynamics"),
                "tau_b"=>Dict("units"=>"Pa", "long_name"=>"Basal Stress (friction)", "group"=>"Dynamics"),
                "fstream"=>Dict("units"=>"--", "long_name"=>"Stream Fraction", "group"=>"Dynamics"),
                "fstream_ref"=>Dict("units"=>"--", "long_name"=>"Stream Fraction", "group"=>"Dynamics"),
                "TMB"=>Dict("units"=>"m/a", "long_name"=>"Total Mass Balance", "group"=>"Thermodynamics"),            # Thermodynamics
                "SMB"=>Dict("units"=>"m/a", "long_name"=>"Surface Mass Balance", "group"=>"Thermodynamics"),               
                "Acc"=>Dict("units"=>"m/a", "long_name"=>"Surface Accumulation", "group"=>"Thermodynamics"),
                "M"=>Dict("units"=>"m/a", "long_name"=>"Surface Melt", "group"=>"Thermodynamics"),
                "T"=>Dict("units"=>"ºC", "long_name"=>"Ice Temperature"),
                "Q_dif"=>Dict("units"=>"K/a", "long_name"=>"Temperature Solver", "group"=>"Thermodynamics"),
                "Q_difup"=>Dict("units"=>"K/a", "long_name"=>"Temperature Solver Term", "group"=>"Thermodynamics"),
                "Q_difdown"=>Dict("units"=>"K/a", "long_name"=>"Temperature Solver Term", "group"=>"Thermodynamics"),
                "Q_adv"=>Dict("units"=>"K/a", "long_name"=>"Advective Heat Dissipation", "group"=>"Thermodynamics"),
                "Q_drag"=>Dict("units"=>"K/a", "long_name"=>"Basal Friction Heating", "group"=>"Thermodynamics"),
                "alpha"=>Dict("units"=>"--", "long_name"=>"Temperate Smoothing Function", "group"=>"Thermodynamics"),
                "T_surf"=>Dict("units"=>"ºC", "long_name"=>"Surface Temperature", "group"=>"Thermodynamics"),
                "T_sl"=>Dict("units"=>"ºC", "long_name"=>"Sea-level Temperature", "group"=>"Thermodynamics"),
                "Hdot"=>Dict("units"=>"m/a", "long_name"=>"Ice Thickness Change", "group"=>"Derivatives"),              # Derivatives
                "Hseddot"=>Dict("units"=>"m/a", "long_name"=>"Sediment Thickness Change", "group"=>"Derivatives"),
                "Bdot"=>Dict("units"=>"m/a", "long_name"=>"Bedrock Elevation Change", "group"=>"Derivatives"),
                "Tdot"=>Dict("units"=>"m/a", "long_name"=>"Ice Temperature Change", "group"=>"Derivatives"),
                "fstreamdot"=>Dict("units"=>"--", "long_name"=>"Stream fraction change", "group"=>"Derivatives")
                )


# OLD CODE

# Help:
#   -- "mutable struct" works similar to fortran classes
#   -- Defining as {T<:Real} make sure that the values are correct
#   -- In this script can be found the struct of parameters and NOW

# mutable struct RunSettings{T<:Real}
#     time_init :: T       # [yr] Starting time (model years)
#     time_end  :: T       # [yr] Ending time (model years)  
#     dt        :: T       # [yr] Loop timestep 
#     dt_out    :: T       # [yr] Frequency of writing
# end

# mutable struct InitialConditions{T<:Real}
#     H_init      :: T     # [m] Initial condition for ice thickness
#     Hsed_init   :: T     # [--] Initial condition for sediments thickness
#     T_init      :: T     # [degC] Initial condition for ice temperature 
#     A_init      :: T     # [] Initial condition for Glenns law coefficient
# end

# mutable struct Parameters{T<:Real}
#     ins_min       :: T          # [] Insolation minimum value (June 21 65N)
#     ins_max       :: T          # [] Insolation maximum value (June 21 65N)
#     ins_prei      :: T          # [] Default preindustrial 
#     co2_prei      :: T          # [ppm] Default preindustrial
#     T_0           :: T          # [degC] Reference freezing temperature
#     T_ref         :: T          # [degC] Reference climatic temperature (-18 works fine with defaults values)  

#     orbital_mode  :: String          # o, op, oe, pe, ope # Tsl formula 
#     P_obl         :: T            # Power of obliquity (normalised to At)
#     tau_obl       :: T            # [yr] Obliquity period
#     P_pre         :: T            # Power of precession (normalised to At)
#     tau_pre       :: T            # [yr] Precession period
#     P_exc         :: T            # Power of excentricity (normalised to At)
#     tau_exc       :: T            # [yr] Excentricity period

#     active_iso    :: Bool         # Switch: include active isostatic rebound?
#     B_eq          :: T            # [m] equilibrium altitude of the bedrock
#     tau_bed       :: T            # [years] relaxation time for the astenosphere
#     T_mantle      :: T            # [degC] Mantle temperature 
#     Q_geo         :: T            # [mW/m²] Geothermal heat flux

#     v_kin         :: T         # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
#     k_1           :: T         # proportionality between Hsed and U
#     f_1           :: T         # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr ==> f1 ~ 1e-6)
#     f_2           :: T         # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
#     fstream_min   :: T         # [--] Mininmal Fraction of the ice sheet considered streaming
#     fstream_max   :: T         # [--] Maximal Fraction of the ice sheet considered streaming
#     L             :: T         # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
#     C_s           :: T         # 1e-7 by default?? #raw sliding parameter (m yr^(-1) Pa^(-2), between 10^(-10) and 10^(-5) in Pollard and deConto 
#     k             :: T         # [J s^(-1) m^(-1) K^(-1)] Ice thermal conductivity, EISMINT value
#     c             :: T         # [J Kg^(-1) K^(-1)] Ice specific heat capacity, EISMINT value 

#     A_m           :: T            # Amplitude of the cosinus for M (surface mass balance) # 0.1 dfault
#     A_t           :: T            # Amplitude of the cosinus for T (surface temperatures)
#     lambda        :: T            # proportionality between positive temperatures and surface melt (m / yr /K )

# end

# mutable struct NOW{T<:Real}     # calculated variables for t = time (it is updated with each time step)
#     time        :: T       # simulation time
#     H           :: T       
#     Hsed        :: T       
#     T           :: T       
#     A           :: T       
#     T_sl        :: T      
#     TMB         :: T       
#     S           :: T       
#     B           :: T              
#     M           :: T       
#     Acc         :: T        
#     SMB         :: T       
#     U_d         :: T       
#     U_b         :: T    
#     U           :: T       
#     T_surf      :: T 
#     tau_b       :: T    
#     tau_d       :: T     
#     Q_dif       :: T   
#     Q_difup     :: T 
#     Q_difdown   :: T        
#     Q_drag      :: T   
#     alpha       :: T   
#     Q_adv       :: T               
#     fstream     :: T
#     fstream_ref :: T
#     Hdot        :: T   
#     Hseddot     :: T  
#     Bdot        :: T   
#     Tdot        :: T   
#     fstreamdot  :: T
#     ins         :: T
#     co2         :: T 
# end





