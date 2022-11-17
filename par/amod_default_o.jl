### Run Settings (CTL)
time_init = -2e6          # [yr] Starting time (model years)
time_end = 2e6            # [yr] Ending time (model years)  
dt = 10.0                       # [yr] Loop timestep 
dt_out = 1000.0                 # [yr] Frequency of writing

### Initial Conditions (INCOND)
H_init = 0.0                    # [m] Initial condition for ice thickness
Hsed_init = 1.0                 # [--] Initial condition for sediments thickness
t_init = -20.0                  # [ºC] Initial condition for ice temperature 
A_init = 1e-16                  # [] Initial condition for the flow parameter of the Glen's flow law

### Run parameters (PAR)
## Radiative parameters
ins_case = "artificial"         # Insolation case: "artificial", "laskar" 
ins_day = 170.0                 # [day] Day in which calculate insolation (170 = June 21)
ins_lat = 65.0                  # [ºN] Latitude in which calculate insolation [-90, 90]ºN
ins_ref = 480.0                   # [W/m²] Reference insolation value for "artificial" case
ins_min = 425.0                 # [W/m²] Insolation minimum value
ins_max = 565.0                 # [W/m²] Insolation maximum value
ins_prei = 480.0                # [W/m²] Default preindustrial 
co2_prei = 280.0                # [ppm] Default preindustrial 
active_radco2 = false           # Switch: include co2 radiative forcing?
A_ins = 80.0                    # [W/m²] Amplitude of artificial insolation forcing
A_t = 25.0                      # [ºC or K] Amplitude of temperature forcing (surface temperatures)

## Orbital parameters
orb_case = "o"                # o, op, oe, pe, ope # Tsl formula 
P_obl = 1.0                     # Power of obliquity (normalised to At)
tau_obl = 41e3                  # [yr] Obliquity period
P_pre = 1.0                     # Power of precession (normalised to At)
tau_pre = 23e3                  # [yr] Precession period
P_exc = 1.0                     # Power of excentricity (normalised to At)
tau_exc = 100e3                 # [yr] Excentricity period

## Geophysical parameters
active_iso = true               # Switch: include active isostatic rebound?
B_eq = 500.0                    # [m] equilibrium altitude of the bedrock
tau_bed = 5000.0                # [years] relaxation time for the astenosphere
t_mantle = 0.0                  # [ºC] Mantle temperature 
H_mantle = 1000.0               # [m] Mantle thickness
Q_geo = 50.0                    # [mW/m²] Geothermal heat flux

## Ice Parameters
# -- dynamics
v_kin = 1000.0                  # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
k_1 = 1e-6                      # proportionality between Hsed and U
f_1 = 1e-8                      # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr ==> f1 ~ 1e-6)
f_2 = 1e-5                      # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
fstream_min = 0.4               # [--] Mininimal Fraction of the ice sheet considered streaming
fstream_max = 0.4               # [--] Maximal Fraction of the ice sheet considered streaming
L = 1e6                         # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
C_s = 1e-7                      # [m yr⁻¹ Pa⁻²] Raw sliding parameter (between 10^(-10) and 10^(-5) in Pollard and deConto (2012) after Schoof (2007)) 1e-7 by default?? 
ud_case = "sia"                 # Ice flow approximation: "sia"
ub_case = "weertmanq"           # parameterization of basal velocity: "weertman"
glen_n = 3.0                    # [--] Glen's flow law exponent

# -- thermodynamics
t_sb = 5.0                      # [ºC] Represents the thermal state of the boundary between deformational and streaming part -- jas 
tsurf_case = "linear"           # Surface temperature method: "linear"
A_m = 0.5                       # Amplitude of the cosinus for M (surface mass balance) # 0.1 default -- jas
k = 2.1                         # [J s⁻¹ m⁻¹ K⁻¹] Ice thermal conductivity, (2.1, EISMINT value from Huybrecths et al. (1996))
cc_case = "cc"                  # Clausius-Clapeyron approximation: "cc", "AERKi" works well between -80 and 0ºC Alduchov and Eskridge (1996)
RH = 0.8                        # [0-1] Relative Humidity -- spm, to calibrate or discuss
k_pr = 50.0                     # [] Precipitation parameter, (Robinson et al. 2010 take 50.0)
tau_w = 5 / 365                 # [yr] Water turnover time in the atmosphere (Robinson et al. 2010 take 5 days)
t_snow = -11.6                  # [ºC] Air temperature below which we consider full snowfall (Bales et al. 2009 take -11.4ºC, Robinson et al. 2010 take -7ºC)
t_rain = 7.4                    # [ºC] Air temperature above which we consider full rain (Bales et al. 2009 take 7.4ºC, Robinson et al. 2010 take 7ºC)
sm_case = "PDD"                 # Surface melting case: "PDD", "ITM"
lambda = 0.1                    # [m yr⁻¹ K⁻¹] Proportionality between positive temperatures and surface melt 
melt_offset = 10.0              # [ºC] Temperature threshold that allows melting  
c = 2009.0                      # [J Kg⁻¹K⁻¹] Ice specific heat capacity, EISMINT value 

