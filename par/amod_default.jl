# Run Settings
time_init = -2000000.0         # [yr] Starting time (model years)
time_end  =  2000000.0         # [yr] Ending time (model years)  
dt        =       10.0         # [yr] Loop timestep 
dt_out    =     1000.0         # [yr] Frequency of writing

# Initial Conditions
H_init        = 0.0            # [m] Initial condition for ice thickness
Hsed_init     = 1.0            # [--] Initial condition for sediments thickness
T_init        = -20.0          # [degC] Initial condition for ice temperature 
A_init        = 1e-16          # [] Initial condition for Glenns law coefficient

# Radiative parameters
ins_min       = 425.0          # [] Insolation minimum value (June 21 65N)
ins_max       = 565.0          # [] Insolation maximum value (June 21 65N)
ins_prei      = 480.0          # [] Default preindustrial 
co2_prei      = 280.0          # [ppm] Default preindustrial 
T_0           = 0.0            # [degC] Reference freezing temperature 
T_ref         = 0.0            # [degC] Reference climatic temperature (-18 works fine with defaults values) 

# Orbital parameters
orbital_mode  = "ope"          # o, op, oe, pe, ope # Tsl formula 
P_obl         = 1.0            # Power of obliquity (normalised to At)
tau_obl       = 41000.0        # [yr] Obliquity period
P_pre         = 0.0            # Power of precession (normalised to At)
tau_pre       = 23000.0        # [yr] Precession period
P_exc         = 0.1            # Power of excentricity (normalised to At)
tau_exc       = 100000.0       # [yr] Excentricity period

# Geophysical parameters
active_iso    = true         # Switch: include active isostatic rebound?
B_eq          = 500.0          # [m] equilibrium altitude of the bedrock
tau_bed       = 5000.0         # [years] relaxation time for the astenosphere
T_mantle      =  0.0           # [degC] Mantle temperature 
Q_geo         = 50.0           # [mW/m²] Geothermal heat flux

## Ice Parameters
# -- dynamics
v_kin         = 1000.0         # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
k_1           = 1e-6           # proportionality between Hsed and U
f_1           = 1e-8           # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr ==> f1 ~ 1e-6)
f_2           = 1e-5           # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
fstream_min  = 0.4             # [--] Mininimal Fraction of the ice sheet considered streaming
fstream_max  = 0.4             # [--] Maximal Fraction of the ice sheet considered streaming
L             = 1e6            # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
C_s           = 1e-7           # 1e-7 by default?? #raw sliding parameter (m yr^(-1) Pa^(-2), between 10^(-10) and 10^(-5) in Pollard and deConto 
k             = 2.1            # [J s^(-1) m^(-1) K^(-1)] Ice thermal conductivity, EISMINT value
c             = 2009.0         # [J Kg^(-1) K^(-1)] Ice specific heat capacity, EISMINT value 

# -- thermodynamics
A_m           = 0.5            # Amplitude of the cosinus for M (surface mass balance) # 0.1 dfault
A_t           = 25.0           # Amplitude of the cosinus for T (surface temperatures)
lambda        = 0.1            # proportionality between positive temperatures and surface melt (m / yr /K )

