Base.@kwdef struct Params
    ### Run Settings
    time_init::Real = -2.0e6          # [yr] Starting time (model years)
    time_end::Real = 0.0              # [yr] Ending time (model years)
    dt_out::Real = 1000.0             # [yr] Frequency of writing

    ### Earth constants
    year_len::Real = 365.2422          # year length in days
    sec_year::Real = 31556926.0        # [s] Seconds in a year, EISMINT value
    P_sl::Real = 101325                # [Pa] Pressure at sea level
    I0::Real = 1365.2                  # [Wm⁻²] Solar constant
    g::Real = 9.81                     # [m/s²] Gravitational acceleration
    degK::Real = 273.15                # [K] kelvin to celsius change 
    Aoc::Real = 3.618e8                # [km²] Oceanic surface
    rhoi::Real = 910.0                 # [kg/m³] Ice density 
    rhow::Real = 1000.0                # [kg/m³] Water density
    rhom::Real = 2700.0                # [kg/m³] Lithosphere density 
    lapse_rate::Real = 0.0065          # [K/m] Atmospheric lapse rate (Γ)
    Tfreezing::Real = 0.0 + degK       # [K] Reference freezing temperature
    Lv::Real = 2.5e6                   # [J kg⁻¹] Latent heat of vaporization http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    Rd::Real = 287.0                   # [J K⁻¹kg⁻¹] Dry-air gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
    Rv::Real = 461.0                   # [J K⁻¹kg⁻¹] Water-vapor gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/

    ### Initial Conditions
    T0::Real = 0.0 + degK              # [K] Initial condition for regional temperature 
    pCO20::Real = 280.0                 # [ppm] Initial condition for pCO2
    iceage0::Real = 0.0                # [a] Initial condition for ice age
    alpha0::Real = 0.2                 # [--] Initial condition for model's albedo
    H0::Real = 0.0                     # [m] Initial condition for ice thickness
    Hsed0::Real = 1.0                  # [m] Initial condition for sediments thickness (0-1)
    B0::Real = 500.0                   # [m] Inital condition for bed elevation
    Tice0::Real = -20.0 + degK         # [K] Initial condition for ice temperature 
    fstream0::Real = 0.4               # [--] Initial condition for streaming fraction

    ### Run parameters 
    # -- Switches
    active_iso::Bool = true               # Switch: include active isostatic rebound?
    active_sed::Bool = true               # Switch: include interactive sediments?
    active_climate::Bool = false          # Switch: include climate routines?
    active_ice::Bool = true               # Switch: include ice sheet dynamics?

    # -- Cases
    I_case::String = "artificial"           # Insolation case: "constant", "artificial", "laskar"
    dyn_case::String = "SIA"                # Ice flow approximation: "sia"
    basal_case::String = "weertmanq"        # parameterization of basal velocity: "weertman"
    A_case::String = "linear"               # Clausius-Clapeyron approximation: "ins", "linear"
    M_case::String = "PDD"                  # Surface melting case: "PDD", "ITM"

    # -- Orbital forcing parameters 
    P_obl::Real = 0.9                     # Power of obliquity (normalised to At)
    tau_obl::Real = 41e3                  # [yr] Obliquity period
    P_pre::Real = 0.0                     # Power of precession (normalised to At)
    tau_pre::Real = 23e3                  # [yr] Precession period
    P_exc::Real = 0.1                     # Power of excentricity (normalised to At)
    tau_exc::Real = 100e3                 # [yr] Excentricity period

    I_day::Real = 170.0                   # Day in which calculate insolation
    I_lat::Real = 65.0                    # [ºN] Latitude in which calculate the insolation in the northern hemisphere 

    # -- Radiative forcing parameters   
    I_const::Real = 400.0               # [W/m²] Insolation value for constant insolation mode
    I_min::Real = 425.0                 # [W/m²] Insolation minimum value 
    I_max::Real = 565.0                 # [W/m²] Insolation maximum value

    I_ref::Real = 480.0                   # [W/m²] Present reference value for insolation 
    pCO2_ref::Real = 280.0                 # [ppm] Reference value for pCO2dot
    tau_pCO2::Real = 10.0                  # [yr] Characteristic time for pCO2 evolution
    cpCO2::Real = 2.0                      # [K] sensitivity of reference temperature to pCO2 Concentration
    ktpCO2::Real = 7.0                     # [ppm/K] Proportionality between temperature and pCO2 forcing

    At::Real = 25.0                       # [ºC or K] Amplitude of temperature forcing (surface temperatures)

    time_anth::Real = 2000.0              # [yr] Year in which we take into account the anthropogenic forcing
    pCO2_anth::Real = 3000.0               # [Gt] Anthropogenic amount of pCO2 produced
    At_anth::Real = 55.0                  # [ºC or K] Amplitude of anthropogenic temperature forcing (sea level temperatures)
    Ac_anth::Real = 20.0                  # [ppm] Amplitude of anthropogenic radiative forcing
    tau_anth::Real = 200e3                # [yr] Relaxation time for anthropogenic forcing

    alpha_land::Real = 0.2               # [--] ground albedo
    alpha_oldice::Real = 0.25            # [--] old ice albedo
    alpha_newice::Real = 0.9             # [--] new ice albedo 
    alpha_slope::Real = 5e-6             # [yr⁻¹] Slope of the albedo - ice age parameterization (1 per 100 kyrs)
    alpha_quad::Real = 1e-10             # [yr⁻²] Sensitivity of the albedo - ice age quadratic parameterization (1 per 100 kyrs)        
    tau_alpha::Real = 1e3                # [yr] Characteristic time of albedo evolution w.r.t reference value

    ci::Real = 0.1                       # [K/Wm²] climate sensitivity to insolation
    cc::Real = 0.65                      # [K/Wm²] climate sensitivity to pCO2
    cz::Real = 0.00685                   # [K/m³] climate sensitivity to ice sheet elevation

    Tref0::Real = 0.0 + degK            # [K] Reference climatic temperature 

    tau_T::Real = 900.0                 # [yr] Characteristic time for temperature evolution w.r.t radiative forcing for northern hemisphere

    # -- Geophysical parameters (northern hemisphere -> n) (southern hemisphere -> s)
    Beq::Real = 500.0                    # [m] equilibrium altitude of the bedrock
    tau_bed::Real = 5000.0               # [years] relaxation time for the astenosphere
    Tmantle::Real = 0.0 + degK           # [K] Mantle temperature 
    Hmantle::Real = 1000.0               # [m] Mantle thickness
    Qgeo::Real = 50.0                    # [mW/m²] Geothermal heat flux

    ## Ice Parameters   (northern hemisphere -> n) (southern hemisphere -> s)
    # -- Geometry
    L::Real = 1e6                       # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
    Eref::Real = 2.5e7                  # [km²] Maximum ice surface ~ [1-2]e7 km² (Northern Hemisphere)

    # -- Dynamics
    A_flow::Real = 1e-16                  # [Pa-3 a−1] Flow parameter of the Glen's flow law
    glen_n::Real = 3.0                    # [--] Glen's flow law exponent
    v_kin::Real = 1000.0                  # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
    tau_kin::Real = L / v_kin             # kinematic wave typical time (time in which the streams are propagated towards the interior of the ice sheet)
    f1::Real = 5.0e-8                     # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr::Real ==> f1 ~ 1e-6)
    f2::Real = 1.0e-6                     # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
    fstream_min::Real = 0.4               # [--] Mininimal Fraction of the ice sheet considered streaming
    fstream_max::Real = 0.4               # [--] Maximal Fraction of the ice sheet considered streaming
    Cs::Real = 1e-7                       # [m yr⁻¹ Pa⁻²] Raw sliding parameter (between 10^(-10) and 10^(-5) in Pollard and deConto (2012) after Schoof (2007)) 1e-7 by default?? 

    # -- Thermodynamics
    Tsb::Real = 5.0 + degK                # [K] Represents the thermal state of the boundary between deformational and streaming part -- jas 
    kt::Real = 2.1                        # [J s⁻¹ m⁻¹ K⁻¹] Ice thermal conductivity, (2.1, EISMINT value from Huybrecths et al. (1996))
    pr_ref::Real = 1.0                    # [m/yr] Reference value for precipitation
    A_pr::Real = 0.5                      # Amplitude of the cosinus for M (surface mass balance) # 0.1 default -- jas
    Tsnow::Real = -11.6 + degK            # [K] Air temperature below which we consider full snowfall (Bales et al. 2009 take -11.4ºC, Robinson et al. 2010 take -7ºC)
    Train::Real = 7.4 + degK              # [K] Air temperature above which we consider full rain (Bales et al. 2009 take 7.4ºC, Robinson et al. 2010 take 7ºC)
    lambda::Real = 0.05                   # [m yr⁻¹ K⁻¹] Proportionality between positive temperatures and surface melt
    Tthreshold::Real = -5.0 + degK        # [K] Temperature threshold that allows melting (default::Real = -5.0ºC)
    c::Real = 2009.0                      # [J Kg⁻¹K⁻¹] Ice specific heat capacity, EISMINT value 

    km::Real = 0.0                        # [m/yr] offset melting in ITM-like calculation
    ki::Real = 0.025                      # [m/yr/Wm²] sensitivity parameter of insolation melting ! 0.006 the default?
    ka::Real = 0.008                      # [m/yr/K] sensitivity parameter of accumulation to temperature (Clasuius clapeyron like) ! 0.004 the default?

    Aref::Real = 0.4                      # [m/yr] Reference Accumulation for northern hemisphere (Aref = Amean + ka * ΔTmean)

    Ate::Real = 20.0                      # [K] Thermal amplitude due to ice extent (Northern Hemisphere)
end

## Insolation (orbital) parameters
global OrbData = Insolation.OrbitalData(Insolation.datadir())