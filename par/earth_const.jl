# =============================
#     Info: 
#          This parameter file contains the terrestrial parameters that are assumed to be constant 
#          The variables are defined as global so that they are accessible by the whole model
# =============================

global year_len = 365.2422          # year length in days
global sec_year = 31556926.0        # [s] Seconds in a year, EISMINT value

global P_sl = 101325                # [Pa] Pressure at sea level
global S_0 = 1365.2                 # [Wm⁻²] Solar constant
global g = 9.81                     # [m/s²] Gravitational acceleration
global A_oc = 3.618*10^8            # [km²] Oceanic surface

global rhoi = 910.0                 # [kg/m³] Ice density 
global rhow = 1000.0                # [kg/m³] Water density
global rhom = 2700.0                # [kg/m³] Lithosphere density 

global grad = 0.0065                # [K/m] Atmospheric lapse rate (Γ)
global degK = 273.15                # [K] kelvin to celsius change 
global t_0 = 0.0                    # [ºC] Reference freezing temperature
global Lv = 2.5e6                   # [J kg⁻¹] Latent heat of vaporization http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
global Rd = 287.0                   # [J K⁻¹kg⁻¹] Dry-air gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
global Rv = 461.0                   # [J K⁻¹kg⁻¹] Water-vapor gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/

# -- changing units to IS
global T_0 = t_0 + degK
