Base.@kwdef mutable struct Params{T<:AbstractFloat}
    # Run Settings
    time_init::T = -3.0e6           # [yr] Starting time (model years)
    time_end::T = 0.0               # [yr] Ending time (model years), t0 in insolation is J2000 epoch
    dt::T = 10.0                    # [yr] Time step if fixed time step is activated
    dt_out::T = 1000.0              # [yr] Frequency of writing
    time_spinup::T = 0.0            # [yr] Time length employed to spinup the model

    # Initial Conditions
    degK::T = 273.15                # [K] kelvin to celsius change 
    T0::T = 0.0 + degK              # [K] Initial condition for regional temperature 
    C0::T = 280.0                   # [ppm] Initial condition for carbon dioxide, C
    iceage0::T = 0.0                # [a] Initial condition for ice age
    albedo0::T = 0.2                # [--] Initial condition for model's albedo
    H0::T = 0.0                     # [m] Initial condition for ice thickness
    Hsed0::T = 30.0                 # [m] Initial condition for sediments thickness
    B0::T = 500.0                   # [m] Inital condition for bed elevation
    Tice0::T = -15.0 + degK         # [K] Initial condition for ice temperature 
    fstr0::T = 0.4                  # [--] Initial condition for streaming fraction

    # Run parameters 
    ## Switches
    active_iso::Bool = true            # Switch: include active isostatic rebound?
    active_sed::Bool = true            # Switch: include interactive sediments?
    active_climate::Bool = true        # Switch: include climate routines?
    active_ice::Bool = true            # Switch: include ice sheet dynamics?
    active_thermo::Bool = false        # Switch: include ice sheet thermodynamics?
    active_aging::Bool = true          # Switch: include ice aging?
    active_snow_on_ice::Bool = true    # Switch: does snow accumulation affect ITM ablation?

    ## Cases
    dt_case::String = "adaptive"       # Time step mode: "adaptive", "fixed"
    insol_case::String = "laskar"      # Insolation case: "constant", "artificial", "laskar", "ISI", "caloric", "input"
    regtemp_case::String = "dynamic"   # Time step mode: "dynamic", "constant", "trend", "slope", "comp"
    carbon_case::String = "dynamic"    # Carbon cycle case: "dynamic", "constant", "trend", "slope", "comp" 
    deformationalflow_case::String = "profile"    # Deformational ice flow approximation: "profile"
    plugflow_case::String = "profile"             # Plug flow approximation: "profile", "constant"    
    deformational_case::String = "glen58"      # parameterization of deformational velocity: "glen58"
    basal_case::String = "weertmanq"   # parameterization of basal velocity: "weertman", "mb"
    albedo_case::String = "prognostic" # Albedo case: "prognostic", "diagnostic"
    snowfall_case::String = "linear"   # Clausius-Clapeyron approximation: "ins", "linear"
    ablation_case::String = "ITM"      # Surface melting case: "PDD", "ITM", "PDD-LIN"
    diffusion_case::String = "2pts"    # Diffusion equation case: "2pts", "3pts"
    streaming_case::String = "fixed"  # Reference streaming case: "fixed", "dynamic"
    ref_streaming_case::String = "thermo"  # Reference streaming case: "theo", "thermo"
    ice_hrz_scale_case::String = "fixed"   # Ice-sheet aspect ratio case: "fixed", "dynamic"
    ice_discharge_case::String = "fixed"   # Ice-sheet discharge case: "fixed", "dynamic"
    temperate_layer_case::String = "fixed"   # Temperate layer thickness case: "fixed", "verb2018"

    ## I, insol, Insolation
    insol_day::T = 170.0            # Day in which calculate insolation
    insol_lat::T = 65.0             # [ºN] Latitude in which calculate the insolation in the northern hemisphere 
    insol_input::String = "data/insolation/solstice_insolation_65N170_10yr_5MyrBP-0.jld2"    # File to read as input if insol_case = "input"
    insol_const::T = 400.0          # [W/m²] Insolation value for constant insolation mode
    insol_min::T = 425.0            # [W/m²] Insolation minimum value 
    insol_max::T = 565.0            # [W/m²] Insolation maximum value
    insol_threshold::T = 300.0      # [W/m²] Insolation threshold value for integrated summer insolation, Huybers (2006) used 275 W/m²
    insol_ref::T = 480.0            # [W/m²] Reference value for insolation. For present day: 480.0 (solstice),  5.8e9 (caloric), 6.7e9 (annual) 

    Ppre::T = 0.0                   # Power of precession (normalised to At)
    taupre::T = 23e3                # [yr] Precession period
    Pobl::T = 0.9                   # Power of obliquity (normalised to At)
    tauobl::T = 41e3                # [yr] Obliquity period
    Pexc::T = 0.1                   # Power of excentricity (normalised to At)
    tauexc::T = 100e3               # [yr] Excentricity period
    At::T = 35.0                    # [ºC or K] Amplitude of temperature forcing (surface temperatures)

    ## T, Regional air temperature
    cI::T = 0.1                     # [Km²/W] climate sensitivity to insolation (laskar method)
    cISI::T = 1e-7                  # [Km²/J] climate sensitivity to Integrated Summer Insolation 
    cCAL::T = 0.17e-7               # [Km²/J] climate sensitivity to Caloric season insolation
    cC::T = 0.65                    # [Km²/W] climate sensitivity to C (CO2, carbon dioxide)
    cZ::T = 0.007                   # [K/m] climate sensitivity to ice sheet elevation
    tauT::T = 900.0                 # [yr] Characteristic time for temperature evolution w.r.t radiative forcing for northern hemisphere
    kT::T = -1.2e-6                 # [K/yr] Temperature imposed slope
    Tref::T = 0.0 + degK            # [K] Initial or constant reference climatic temperature 
    deltaT::T = 0.0                 # [K] Temperature change in a given period of time (regtemp_case=trend, associated with kT)

    ## C, CO2, Carbon dioxide
    Cref::T = 280.0                 # [ppm] Initial or reference value for C
    tauC::T = 10.0                  # [yr] Characteristic time for C evolution
    kCT::T = 2.0                    # [K] sensitivity of reference temperature to C anthorpogenic input (conversion C to T) -- perhaps kCT and kTC should be the same? spm 
    kTC::T = 7.0                    # [ppm/K] Proportionality between temperature and C forcing (conversion T to C)
    kC::T = -1e-5                   # [ppm/yr] Carbon dioxide imposed slope
    deltaC::T = 0.0                 # [ppm] Carbon dioxide change in a given period of time (carbon_case=trend, associated with kC)

    time_anth::T = 2000.0           # [yr] Year in which we take into account the anthropogenic forcing
    C_anth::T = 3000.0              # [Gt] Anthropogenic amount of CO2 produced
    AT_anth::T = 55.0               # [ºC or K] Amplitude of anthropogenic temperature forcing (sea level temperatures)
    AC_anth::T = 20.0               # [ppm] Amplitude of anthropogenic radiative forcing
    tau_anth::T = 200e3             # [yr] Relaxation time for anthropogenic forcing

    ## A, iceage, ice age
    # no parameters

    ## α, albedo, Albedo
    albedo_land::T = 0.2            # [--] ground albedo
    albedo_oldice::T = 0.4          # [--] old ice albedo
    albedo_newice::T = 0.9          # [--] new ice albedo 
    k_albedo::T = 5e-6              # [yr⁻¹] Slope of the albedo - ice age parameterization (1 per 100 kyrs)
    albedo_quad::T = 1e-10          # [yr⁻²] Sensitivity of the albedo - ice age quadratic parameterization (1 per 100 kyrs)        
    tau_albedo::T = 1e3             # [yr] Characteristic time of albedo evolution w.r.t reference value

    ## H, Ice thickness
    ### ṡ, snowfall
    pr_ref::T = 1.0                 # [m/yr] Reference value for precipitation
    Apr::T = 0.5                    # Amplitude of the cosinus for M (surface mass balance) # 0.1 default -- jas
    Tsnow::T = -11.6 + degK         # [K] Air temperature below which we consider full snowfall (Bales et al. 2009 take -11.4ºC, Robinson et al. 2010 take -7ºC)
    Train::T = 7.4 + degK           # [K] Air temperature above which we consider full rain (Bales et al. 2009 take 7.4ºC, Robinson et al. 2010 take 7ºC)
    sref::T = 0.3                   # [m/yr] Reference Accumulation for northern hemisphere (sref = smean + ka * ΔTmean)
    ks::T = 0.02                    # [m/yr/K] sensitivity parameter of snowfall accumulation to temperature (Clasuius clapeyron like) ! 0.004 the default?

    ### ȧ, ablation
    lambda::T = 0.1                 # [m yr⁻¹ K⁻¹] Proportionality between positive temperatures and surface melt
    Tthreshold::T = -5.0 + degK     # [K] Temperature threshold that allows melting (default::T = -5.0ºC)
    km::T = 0.0                     # [m/yr] offset melting in ITM-like calculation
    kI::T = 0.027                   # [m/yr/Wm²] sensitivity parameter of insolation melting ! 0.006 the default?

    ### Ice sheet dynamics
    Ath::T = 15.0                   # [K] Thermal amplitude due to ice-sheet area (Northern Hemisphere)
    L::T = 1e6                      # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
    hrz_scale_coeff::T = 0.1        # [m⁻¹] Coefficient for diagnostic L
    hrz_scale_exp::T = 2.0          # [--] Exponent for diagnostic L
    hrz_scale_lb::T = 200e3         # [m] Lower bound for dynamic L: ~√50e3 [km] (rounded ice cap) 
    hrz_scale_ub::T = 3e6           # [m] Upper bound for dynamic L: ~√(13-14e6) [km] (LIS/AIS), just check it does not saturate results!
    hrz_vel_scale::T = 300.0        # [m yr⁻¹] Typical scale of sliding velocities
    vrt_vel_scale::T = -0.01        # [m yr⁻¹] Typical scale of vertical velocities in the base
    Aflow::T = 1e-16                # [Pa-3 a−1] Flow parameter of the Glen's flow law
    glen_n::T = 3.0                 # [--] Glen's flow law exponent
    Cs::T = 5e-5                    # [m yr⁻¹ Pa⁻²] Raw sliding parameter (between 10^(-10) and 10^(-5) in Pollard and deConto (2012) after Schoof (2007)) 1e-7 by default?? 
    beta_min::T = 0.0               # [--] Minimal basal coefficient allowed
    sintheta::T = 1e-3              # [--] Constant ice surface slope for plug flow constant approximation

    ## Hsed, sediment layer thickness
    Hsed_max::T = 30.0              # [m] Maximum amount of sediments (pre Pleistocene reconstructions ~30 m, Clark et al., 2006)
    Hsed_min::T = 5.0               # [m] Minimum amount of sediments (the mode of PD distribution is ~5 m, inferred from Laske and Masters, 1997)
    quarrying_frac::T = 1.5e-7                  # [--] fv, fraction of the sediments that is removed because of v # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr::Real ==> fv ~ 1e-6)
    weathering_frac::T = 5.0e-6                  # [--] fp, fraction of the ablation (a) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr

    ## B, bedrock elevation
    Beq::T = 500.0                  # [m] equilibrium altitude of the bedrock
    taubedrock::T = 5000.0          # [years] relaxation time for the astenosphere

    ## Tice, ice temperature
    hgeo::T = 50.0e-3               # [W/m²] Geothermal heat flux
    Tstr::T = 5.0 + degK            # [K] Streaming boundary temperature, Represents the thermal state of the boundary between deformational and streaming part -- jas 
    cice::T = 2009.0                # [J Kg⁻¹K⁻¹] Ice specific heat capacity, EISMINT value 
    kthr::T = 2.1                   # [J s⁻¹ m⁻¹ K⁻¹] Ice thermal conductivity, (2.1, EISMINT value from Huybrecths et al. (1996))
    Tmp::T = degK                   # [K] melting point temperature
    basal_scale::T = 10.0           # [m] Vertical scale of the basal interface (Robel et al., 2013, JGR)
    Pe::T = 10.0                    # [--] Peclet number

    ## fstr, streaming fraction
    fstrmin::T = 0.4                # [--] Mininimal Fraction of the ice sheet considered streaming
    fstrmax::T = 0.4                # [--] Maximal Fraction of the ice sheet considered streaming
    vkin::T = 1000.0                # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (1500 m/yr, Payne 2004)
    taukin::T = L / vkin            # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)    

    ### Earth constants
    year_len::T = 365.2422          # year length in days
    sec_year::T = 31556926.0        # [s/yr] Seconds in a year, EISMINT value
    I0::T = 1365.2                  # [Wm⁻²] Solar constant
    g::T = 9.81                     # [m/s²] Gravitational acceleration
    Surfoc::T = 3.618e8             # [km²] Oceanic surface (Cogley et al., 2012)
    rhoice::T = 910.0               # [kg/m³] Ice density 
    muice0::T = 1.928e-5            # [Pa⋅s] Ice reference viscosity
    rhowater::T = 1000.0            # [kg/m³] Water density
    rhobed::T = 2700.0                # [kg/m³] Lithosphere density 
    Γ::T = 0.0065                   # [K/m] Atmospheric lapse rate (Γ)
    Tfreezing::T = 0.0 + degK       # [K] Reference freezing temperature
    Lv::T = 2.5e6                   # [J kg⁻¹] Latent heat of vaporization http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    Rd::T = 287.0                   # [J K⁻¹kg⁻¹] Dry-air gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    Rv::T = 461.0                   # [J K⁻¹kg⁻¹] Water-vapor gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/

    # Ice-sheet discretization constraints
    ice_exists_thr::T = 0.0            # [m] strict threshold where ice sheet does not exist
    ice_is_big_thr::T = 100.0           # [m] threshold for ice sheet to be considered big enough 
    ice_is_old_thr::T = 10.0           # [yr] threshold for ice sheet to be considered old enough 
end

## Insolation (orbital) parameters
global OrbData = Insolation.OrbitalData()   # This global variable stores the Orbital Parameters of Insolation.jl