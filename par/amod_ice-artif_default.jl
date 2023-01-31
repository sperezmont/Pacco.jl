### Run Settings (CTL)
time_init = -2e6                # [yr] Starting time (model years)
time_end = 0                  # [yr] Ending time (model years)  
dt = 10.0                       # [yr] Loop timestep 
dt_out = 1000.0                 # [yr] Frequency of writing
hemisphere = ["n"]              # list with hemispheres to use: ["n"] for northern, ["s"] for southern or ["n", "s"] for both hemispheres

### Initial Conditions (INCOND) (northern hemisphere -> n) (southern hemisphere -> s)
H_init_n = 0.0                  # [m] Initial condition for ice thickness 
H_init_s = 0.0                  # 
Hsed_init_n = 1.0               # [--] Initial condition for sediments thickness
Hsed_init_s = 1.0               # 
B_init_n = 0.0                  # [m] Initial condition for bedrock elevation
B_init_s = 0.0                  # 
t_init_n = 0.0                  # [ºC] Initial condition for regional temperature 
t_init_s = -5.0                 # 
t_ice_init_n = -20.0            # [ºC] Initial condition for ice temperature 
t_ice_init_s = -20.0            # 
A_init = 1e-16                  # [yr⁻¹Pa⁻³??] Initial condition for the flow parameter of the Glen's flow law

### Run parameters (PAR)
# -- Switches
active_outout = false           # Switch: generate out.out?
active_iso = true               # Switch: include active isostatic rebound?
active_sed = true              # Switch: include interactive sediments?
active_climate = false           # Switch: include climate routines?
active_ice = true              # Switch: include ice sheet dynamics?

# -- Cases
ins_case = "artificial"         # Insolation case: "artificial", "laskar" 
ud_case = "sia"                 # Ice flow approximation: "sia"
ub_case = "weertmanq"           # parameterization of basal velocity: "weertman"
tsurf_case = "linear"           # Surface temperature method: "linear"
ac_case = "ins"                 # Clausius-Clapeyron approximation: "ins", "linear"
sm_case = "PDD"                 # Surface melting case: "PDD", "ITM"

# -- Orbital forcing parameters 
P_obl = 0.9                     # Power of obliquity (normalised to At)
tau_obl = 41e3                  # [yr] Obliquity period
P_pre = 0.0                     # Power of precession (normalised to At)
tau_pre = 23e3                  # [yr] Precession period
P_exc = 0.1                     # Power of excentricity (normalised to At)
tau_exc = 100e3                 # [yr] Excentricity period

ins_day = 170.0                 # Day in which calculate insolation
ins_lat_n = 65.0                # [ºN] Latitude in which calculate the insolation in the northern hemisphere 
ins_lat_s = -60.0               # [ºN] Latitude in which calculate the insolation in the southern hemisphere

# -- Radiative forcing parameters   (northern hemisphere -> n) (southern hemisphere -> s)
ins_min = 425.0                 # [W/m²] Insolation minimum value 
ins_max = 565.0                 # [W/m²] Insolation maximum value

ins_ref_n = 480.0                 # [W/m²] Present reference value for insolation (Northern Hemisphere)
ins_ref_s = 526.0                 # [W/m²] Present reference value for insolation (Southern Hemisphere)
ins_prei = 480.0                 # [W/m²] Default preindustrial 
co2_prei = 280.0                 # [ppm] Default preindustrial 
co2_ref = 280.0                  # [ppm] Reference value for co2dot
tau_co2 = 10.0                  # [yr] Characteristic time for co2 evolution
ktco2 = 10.0                    # [ppm/K] Proportionality between temperature and co2 forcing

A_t = 20.0                      # [ºC or K] Amplitude of temperature forcing (surface temperatures)

time_anth = 2000.0              # [yr] Year in which we take into account the anthropogenic forcing
co2_anth = 100.0                # [ppm] Nthropogenic amoun of co2 produced
At_anth = 55.0                  # [ºC or K] Amplitude of anthropogenic temperature forcing (sea level temperatures)
Ac_anth = 20.0                   # [ppm] Amplitude of anthropogenic radiative forcing
tau_anth = 200e3                # [yr] Relaxation time for anthropogenic forcing

albedo_land = 0.2               # [--] ground albedo
albedo_newice = 0.9             # [--] new ice albedo 
albedo_slope = 5e-6             # [yr⁻¹] Slope of the albedo - ice age parameterization (1 per 100 kyrs)
albedo_quad = 1e-10             # [yr⁻²] Sensitivity of the albedo - ice age quadratic parameterization (1 per 100 kyrs)        
tau_albedo = 1e3                # [yr] Characteristic time of albedo evolution w.r.t reference value

cs = 0.15                       # [K/Wm²] climate sensitivity
csz = 0.0065                    # [K/m³] climate sensitivity to ice sheet elevation

t_ref_n = 0.0                    # [ºC] Reference climatic temperature for northern hemisphere
t_ref_s = -5.0                   # [ºC] Reference climatic temperature for southern hemisphere

tau_rf_n = 900.0                 # [yr] Characteristic time for temperature evolution w.r.t radiative forcing for northern hemisphere
tau_rf_s = 1000.0                # [yr] Characteristic time for temperature evolution w.r.t radiative forcing for southern hemisphere

# -- Geophysical parameters (northern hemisphere -> n) (southern hemisphere -> s)
B_eq_n = 500.0                    # [m] equilibrium altitude of the bedrock
B_eq_s = 500.0                    # [m] equilibrium altitude of the bedrock
tau_bed_n = 5000.0                # [years] relaxation time for the astenosphere
tau_bed_s = 5000.0                # [years] relaxation time for the astenosphere
t_mantle_n = 0.0                  # [ºC] Mantle temperature 
t_mantle_s = 0.0                  # [ºC] Mantle temperature 
H_mantle_n = 1000.0               # [m] Mantle thickness
H_mantle_s = 1000.0               # [m] Mantle thickness
Q_geo = 50.0                    # [mW/m²] Geothermal heat flux

## Ice Parameters   (northern hemisphere -> n) (southern hemisphere -> s)
# -- Geometry
L = 1e6                         # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
E_ref_n = 2.5e7                  # [km²] Maximum ice surface ~ [1-2]e7 km² (Northern Hemisphere)
E_ref_s = 14e6                   # [km²] Maximum ice surface ~ [1-2]e7 km² (Southern Hemisphere)

# -- Dynamics
v_kin = 1000.0                  # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
f_1 = 1.5e-8                    # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr ==> f1 ~ 1e-6)
f_2 = 1e-5                      # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
fstream_min_n = 0.4               # [--] Mininimal Fraction of the ice sheet considered streaming
fstream_min_s = 0.4               # [--] Mininimal Fraction of the ice sheet considered streaming
fstream_max_n = 0.4               # [--] Maximal Fraction of the ice sheet considered streaming
fstream_max_s = 0.4               # [--] Maximal Fraction of the ice sheet considered streaming
C_s = 1e-7                      # [m yr⁻¹ Pa⁻²] Raw sliding parameter (between 10^(-10) and 10^(-5) in Pollard and deConto (2012) after Schoof (2007)) 1e-7 by default?? 
glen_n = 3.0                    # [--] Glen's flow law exponent

# -- Thermodynamics
t_sb = 5.0                      # [ºC] Represents the thermal state of the boundary between deformational and streaming part -- jas 
kt = 2.1                        # [J s⁻¹ m⁻¹ K⁻¹] Ice thermal conductivity, (2.1, EISMINT value from Huybrecths et al. (1996))
pr_ref = 0.7                      # [m/yr] Reference value for precipitation
A_pr = 0.5                      # Amplitude of the cosinus for M (surface mass balance) # 0.1 default -- jas
t_snow = -11.6                  # [ºC] Air temperature below which we consider full snowfall (Bales et al. 2009 take -11.4ºC, Robinson et al. 2010 take -7ºC)
t_rain = 7.4                    # [ºC] Air temperature above which we consider full rain (Bales et al. 2009 take 7.4ºC, Robinson et al. 2010 take 7ºC)
lambda = 0.07                   # [m yr⁻¹ K⁻¹] Proportionality between positive temperatures and surface melt 
melt_offset = -5.0              # [ºC or K] Temperature threshold that allows melting (default = -5.0ºC)
c = 2009.0                      # [J Kg⁻¹K⁻¹] Ice specific heat capacity, EISMINT value 

km = 0.0                        # [m/yr] offset melting in ITM-like calculation
ki = 0.009                      # [m/yr/Wm²] sensitivity parameter of insolation melting ! 0.006 the default?
ka = 0.008                      # [m/yr/K] sensitivity parameter of accumulation to temperature (Clasuius clapeyron like) ! 0.004 the default?

Acc_ref_n = 0.1                 # [m/yr] Reference Accumulation for northern hemisphere
Acc_ref_s = 0.1                 # [m/yr] Reference Accumulation for southern hemisphere

A_te_n = 20.0                   # [K] Thermal amplitude due to ice extent (Northern Hemisphere)
A_te_s = 20.0                   # [K] Thermal amplitude due to ice extent (Southern Hemisphere)

## Proxy files
# -- variable names must be "T", "T_lo", "T_up"
T_proxy = ["T_barker-etal_2011.nc"]
           #"T_snyder_2016.nc"]

# -- variable names must be "co2", "co2_lo", "co2_up"
co2_proxy = ["co2_luthi-etal_2008.nc"]

# -- variable names must be "V", "V_lo", "V_up"
V_proxy = [#"V_waelbroeck-etal_2002.nc",
           "V_spratt-lisiecki_2016.nc"]