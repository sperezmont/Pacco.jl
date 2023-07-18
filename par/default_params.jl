Base.@kwdef struct Params
    # Run Settings
    time_init::Real = -1.5e6           # [yr] Starting time (model years)
    time_end::Real = 0.0               # [yr] Ending time (model years)
    dt::Real = 10.0                    # [yr] Time step if fixed time step is activated
    dt_out::Real = 1000.0              # [yr] Frequency of writing
    time_spinup::Real = 5e5            # [yr] Time length employed to spinup the model

    # Initial Conditions
    degK::Real = 273.15                # [K] kelvin to celsius change 
    T0::Real = 0.0 + degK              # [K] Initial condition for regional temperature 
    C0::Real = 280.0                   # [ppm] Initial condition for carbon dioxide, C
    iceage0::Real = 0.0                # [a] Initial condition for ice age
    albedo0::Real = 0.2                # [--] Initial condition for model's albedo
    H0::Real = 0.0                     # [m] Initial condition for ice thickness
    Hsed0::Real = 0.0                  # [m] Initial condition for sediments thickness (0-1)
    B0::Real = 500.0                   # [m] Inital condition for bed elevation
    Tice0::Real = -20.0 + degK         # [K] Initial condition for ice temperature 
    fstr0::Real = 0.4                  # [--] Initial condition for streaming fraction

    Tref0::Real = 0.0 + degK           # [K] Reference climatic temperature 

    # Run parameters 
    ## Switches
    active_iso::Bool = true            # Switch: include active isostatic rebound?
    active_sed::Bool = false           # Switch: include interactive sediments?
    active_climate::Bool = true        # Switch: include climate routines?
    active_ice::Bool = true            # Switch: include ice sheet dynamics?
    active_aging::Bool = true          # Switch: include ice aging?

    ## Cases
    dt_case::String = "adaptive"       # Time step mode: "adaptive", "fixed"
    insol_case::String = "laskar"      # Insolation case: "constant", "artificial", "laskar", "ISI", "caloric", "input"
    carbon_case::String = "dynamic"    # Carbon cycle case: "dynamic", "constant", "trended" 
    dyn_case::String = "SIA"           # Ice flow approximation: "sia"
    basal_case::String = "weertmanq"   # parameterization of basal velocity: "weertman"
    snowfall_case::String = "linear"   # Clausius-Clapeyron approximation: "ins", "linear"
    ablation_case::String = "ITM"      # Surface melting case: "PDD", "ITM"

    ## I, insol, Insolation
    insol_day::Real = 170.0            # Day in which calculate insolation
    insol_lat::Real = 65.0             # [ºN] Latitude in which calculate the insolation in the northern hemisphere 
    insol_input::String = "data/insolation/solstice_insolation_65N170_10yr_5MyrBP-0.jld2"    # File to read as input if insol_case = "input"
    insol_const::Real = 400.0          # [W/m²] Insolation value for constant insolation mode
    insol_min::Real = 425.0            # [W/m²] Insolation minimum value 
    insol_max::Real = 565.0            # [W/m²] Insolation maximum value
    insol_threshold = 300.0            # [W/m²] Insolation threshold value for integrated summer insolation, Huybers (2006) used 275 W/m²
    insol_ref::Real = 480.0            # [W/m²] Reference value for insolation. For present day: 480.0 (solstice),  367.0 (caloric), 214.0 (annual) 

    Pobl::Real = 0.9                   # Power of obliquity (normalised to At)
    tauobl::Real = 41e3                # [yr] Obliquity period
    Ppre::Real = 0.0                   # Power of precession (normalised to At)
    taupre::Real = 23e3                # [yr] Precession period
    Pexc::Real = 0.1                   # Power of excentricity (normalised to At)
    tauexc::Real = 100e3               # [yr] Excentricity period
    At::Real = 35.0                    # [ºC or K] Amplitude of temperature forcing (surface temperatures)

    ## T, Regional air temperature
    cI::Real = 0.1                     # [K/Wm²] climate sensitivity to insolation (laskar method)
    cISI::Real = 0.1                   # [K/Wm²] climate sensitivity to Integrated Summer Insolation 
    cCAL::Real = 0.5                   # [K/Wm²] climate sensitivity to Caloric season insolation
    cC::Real = 0.65                    # [K/Wm²] climate sensitivity to C (CO2, carbon dioxide)
    cZ::Real = 0.00685                 # [K/m³] climate sensitivity to ice sheet elevation
    tauT::Real = 900.0                 # [yr] Characteristic time for temperature evolution w.r.t radiative forcing for northern hemisphere

    ## C, CO2, Carbon dioxide
    Cref::Real = 280.0                 # [ppm] Reference value for Cdot
    tauC::Real = 10.0                  # [yr] Characteristic time for C evolution
    kCT::Real = 2.0                    # [K] sensitivity of reference temperature to C anthorpogenic input (conversion C to T) -- perhaps kCT and kTC should be the same? spm 
    kTC::Real = 7.0                    # [ppm/K] Proportionality between temperature and C forcing (conversion T to C)
    kC::Real = -1e-5                   # [ppm/yr] Carbon dioxide imposed slope

    time_anth::Real = 2000.0           # [yr] Year in which we take into account the anthropogenic forcing
    C_anth::Real = 3000.0              # [Gt] Anthropogenic amount of CO2 produced
    AT_anth::Real = 55.0               # [ºC or K] Amplitude of anthropogenic temperature forcing (sea level temperatures)
    AC_anth::Real = 20.0               # [ppm] Amplitude of anthropogenic radiative forcing
    tau_anth::Real = 200e3             # [yr] Relaxation time for anthropogenic forcing

    ## A, iceage, ice age
    # no parameters

    ## α, albedo, Albedo
    albedo_land::Real = 0.2            # [--] ground albedo
    albedo_oldice::Real = 0.25         # [--] old ice albedo
    albedo_newice::Real = 0.9          # [--] new ice albedo 
    k_albedo::Real = 5e-6              # [yr⁻¹] Slope of the albedo - ice age parameterization (1 per 100 kyrs)
    albedo_quad::Real = 1e-10          # [yr⁻²] Sensitivity of the albedo - ice age quadratic parameterization (1 per 100 kyrs)        
    tau_albedo::Real = 1e3             # [yr] Characteristic time of albedo evolution w.r.t reference value

    ## H, Ice thickness
    ### ṡ, snowfall
    pr_ref::Real = 1.0                 # [m/yr] Reference value for precipitation
    Apr::Real = 0.5                    # Amplitude of the cosinus for M (surface mass balance) # 0.1 default -- jas
    Tsnow::Real = -11.6 + degK         # [K] Air temperature below which we consider full snowfall (Bales et al. 2009 take -11.4ºC, Robinson et al. 2010 take -7ºC)
    Train::Real = 7.4 + degK           # [K] Air temperature above which we consider full rain (Bales et al. 2009 take 7.4ºC, Robinson et al. 2010 take 7ºC)
    sref::Real = 0.3                   # [m/yr] Reference Accumulation for northern hemisphere (sref = smean + ka * ΔTmean)
    
    ### ȧ, ablation
    lambda::Real = 0.064               # [m yr⁻¹ K⁻¹] Proportionality between positive temperatures and surface melt
    Tthreshold::Real = -5.0 + degK     # [K] Temperature threshold that allows melting (default::Real = -5.0ºC)
    km::Real = 0.0                     # [m/yr] offset melting in ITM-like calculation
    kI::Real = 0.025                   # [m/yr/Wm²] sensitivity parameter of insolation melting ! 0.006 the default?
    ks::Real = 0.02                    # [m/yr/K] sensitivity parameter of snowfall accumulation to temperature (Clasuius clapeyron like) ! 0.004 the default?

    ### Ice sheet dynamics
    Ath::Real = 20.0                   # [K] Thermal amplitude due to ice-sheet area (Northern Hemisphere)
    L::Real = 1e6                      # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
    Aflow::Real = 1e-16                # [Pa-3 a−1] Flow parameter of the Glen's flow law
    glen_n::Real = 3.0                 # [--] Glen's flow law exponent
    Cs::Real = 1e-7                    # [m yr⁻¹ Pa⁻²] Raw sliding parameter (between 10^(-10) and 10^(-5) in Pollard and deConto (2012) after Schoof (2007)) 1e-7 by default?? 

    ## Hsed, sediment layer thickness
    Hsed_max::Real = 100.0             # [m] Maximum amount of sediments
    Hsed_min::Real = 0.0               # [m] Minimum amount of sediments
    f1::Real = 1.0e-7                  # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr::Real ==> f1 ~ 1e-6)
    f2::Real = 1.0e-6                  # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr

    ## B, bedrock elevation
    Beq::Real = 500.0                  # [m] equilibrium altitude of the bedrock
    taubedrock::Real = 5000.0          # [years] relaxation time for the astenosphere

    ## Tice, ice temperature
    Tmantle::Real = 0.0 + degK         # [K] Mantle temperature 
    Hmantle::Real = 1000.0             # [m] Mantle thickness
    Qgeo::Real = 50.0                  # [mW/m²] Geothermal heat flux
    Tsb::Real = 5.0 + degK             # [K] Streaming boundary temperature, Represents the thermal state of the boundary between deformational and streaming part -- jas 

    ## fstr, streaming fraction
    fstrmin::Real = 0.4                # [--] Mininimal Fraction of the ice sheet considered streaming
    fstrmax::Real = 0.4                # [--] Maximal Fraction of the ice sheet considered streaming
    vkin::Real = 1000.0                # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
    taukin::Real = L / vkin            # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)    

    ### Earth constants
    year_len::Real = 365.2422          # year length in days
    sec_year::Real = 31556926.0        # [s] Seconds in a year, EISMINT value
    I0::Real = 1365.2                  # [Wm⁻²] Solar constant
    g::Real = 9.81                     # [m/s²] Gravitational acceleration
    Surfoc::Real = 3.618e8                # [km²] Oceanic surface
    rhoi::Real = 910.0                 # [kg/m³] Ice density 
    rhow::Real = 1000.0                # [kg/m³] Water density
    rhom::Real = 2700.0                # [kg/m³] Lithosphere density 
    Γ::Real = 0.0065          # [K/m] Atmospheric lapse rate (Γ)
    Tfreezing::Real = 0.0 + degK       # [K] Reference freezing temperature
    Lv::Real = 2.5e6                   # [J kg⁻¹] Latent heat of vaporization http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    Rd::Real = 287.0                   # [J K⁻¹kg⁻¹] Dry-air gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    Rv::Real = 461.0                   # [J K⁻¹kg⁻¹] Water-vapor gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    cice::Real = 2009.0                   # [J Kg⁻¹K⁻¹] Ice specific heat capacity, EISMINT value 
    kth::Real = 2.1                        # [J s⁻¹ m⁻¹ K⁻¹] Ice thermal conductivity, (2.1, EISMINT value from Huybrecths et al. (1996))

end

## Insolation (orbital) parameters
global OrbData = Insolation.OrbitalData()   # This global variable stores the Orbital Parameters of Insolation.jl