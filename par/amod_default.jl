### Run Settings (CTL)
time_init = -2e6                # [yr] Starting time (model years)
time_end = 2e6                  # [yr] Ending time (model years)  
dt = 10.0                       # [yr] Loop timestep 
dt_out = 1000.0                 # [yr] Frequency of writing

### Initial Conditions (INCOND)
H_init = 0.0                    # [m] Initial condition for ice thickness
Hsed_init = 1.0                 # [--] Initial condition for sediments thickness
B_init = 0.0                    # [m] Initial condition for bedrock elevation
t_init = -20.0                  # [ºC] Initial condition for ice temperature 
A_init = 1e-16                  # [yr⁻¹Pa⁻³??] Initial condition for the flow parameter of the Glen's flow law

### Run parameters (PAR)
# -- Switches
active_outout = false           # Switch: generate out.out?
active_radco2 = false           # Switch: include co2 radiative forcing?
active_iso = true               # Switch: include active isostatic rebound?
active_sed = true               # Swithc: include interactive sediments?

# -- Cases
ins_case = "artificial"         # Insolation case: "artificial", "laskar" 
ud_case = "sia"                 # Ice flow approximation: "sia"
ub_case = "weertmanq"           # parameterization of basal velocity: "weertman"
tsurf_case = "linear"           # Surface temperature method: "linear"
cc_case = "ins"                 # Clausius-Clapeyron approximation: "ins", "ARM", "cc"
sm_case = "PDD"                 # Surface melting case: "PDD", "ITM"

# -- Orbital forcing parameters
P_obl = 0.9                     # Power of obliquity (normalised to At)
tau_obl = 41e3                  # [yr] Obliquity period
P_pre = 0.0                     # Power of precession (normalised to At)
tau_pre = 23e3                  # [yr] Precession period
P_exc = 0.1                     # Power of excentricity (normalised to At)
tau_exc = 100e3                 # [yr] Excentricity period

# -- Radiative forcing parameters
ins_day = 170.0                 # [day] Day in which calculate insolation (170 = June 21)
ins_lat = 65.0                  # [ºN] Latitude in which calculate insolation [-90, 90]ºN

ins_min = 425.0                 # [W/m²] Insolation minimum value 
ins_max = 565.0                 # [W/m²] Insolation maximum value

ins_prei = 480.0                # [W/m²] Default preindustrial 
co2_prei = 280.0                # [ppm] Default preindustrial 

A_t = 35.0                      # [ºC or K] Amplitude of temperature forcing (surface temperatures)

time_ant = 2000.0               # [yr] Year in which we take into account the anthropogenic forcing
A_ant = 55.0                    # [ºC or K] Amplitude of anthropogenic temperature forcing (sea level temperatures)
tau_ant = 200e3                 # [yr] Relaxation time for anthropogenic forcing

# -- Geophysical parameters
B_eq = 500.0                    # [m] equilibrium altitude of the bedrock
tau_bed = 5000.0                # [years] relaxation time for the astenosphere
t_mantle = 0.0                  # [ºC] Mantle temperature 
H_mantle = 1000.0               # [m] Mantle thickness
Q_geo = 50.0                    # [mW/m²] Geothermal heat flux

## Ice Parameters
# -- Dynamics
L = 1e6                         # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
v_kin = 1000.0                  # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
f_1 = 0.7e-8                    # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr ==> f1 ~ 1e-6)
f_2 = 1e-5                      # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
fstream_min = 0.4               # [--] Mininimal Fraction of the ice sheet considered streaming
fstream_max = 0.4               # [--] Maximal Fraction of the ice sheet considered streaming
C_s = 1e-7                      # [m yr⁻¹ Pa⁻²] Raw sliding parameter (between 10^(-10) and 10^(-5) in Pollard and deConto (2012) after Schoof (2007)) 1e-7 by default?? 
glen_n = 3.0                    # [--] Glen's flow law exponent

# -- Thermodynamics
t_sb = 5.0                      # [ºC] Represents the thermal state of the boundary between deformational and streaming part -- jas 
kt = 2.1                        # [J s⁻¹ m⁻¹ K⁻¹] Ice thermal conductivity, (2.1, EISMINT value from Huybrecths et al. (1996))
pr_ref = 1                      # [m/yr] Reference value for precipitation
A_pr = 0.5                      # Amplitude of the cosinus for M (surface mass balance) # 0.1 default -- jas
e_0 = 332.41                    # [Pa] constant parameter for ARM CC approximation
e_1 = 17.65                     # [Pa] constant parameter 
RH = 0.8                        # [0-1] Relative Humidity -- spm, to calibrate or discuss
k_pr = 50.0                     # [] Precipitation parameter, (Robinson et al. 2010 take 50.0)
tau_w = 5 / 365                 # [yr] Water turnover time in the atmosphere (Robinson et al. 2010 take 5 days)
t_snow = -11.6                  # [ºC] Air temperature below which we consider full snowfall (Bales et al. 2009 take -11.4ºC, Robinson et al. 2010 take -7ºC)
t_rain = 7.4                    # [ºC] Air temperature above which we consider full rain (Bales et al. 2009 take 7.4ºC, Robinson et al. 2010 take 7ºC)
lambda = 0.1                    # [m yr⁻¹ K⁻¹] Proportionality between positive temperatures and surface melt 
melt_offset = -5.0             # [ºC or K] Temperature threshold that allows melting (default = -5.0ºC)
c = 2009.0                      # [J Kg⁻¹K⁻¹] Ice specific heat capacity, EISMINT value 

