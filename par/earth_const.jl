# =============================
#     Info: 
#          This parameter file contains the terrestrial parameters that are assumed to be constant 
#          The variables are defined as global so that they are accessible by the whole model
# =============================

global epoch = 2000.0               # year defined as present day (default -> January 1, 2000)
global year_len = 360.0             # year length in days (360 as suggested by Alexander Robinson and Mahé Perrette, insolation.f90, 2014)
global sec_year = 31556926.0        # [s] Seconds in a year, EISMINT value

global P_sl = 101325                # [Pa] Pressure at sea level
global S_0 = 1365.0                 # [Wm⁻²] Solar constant
global g = 9.81                     # [m/s²] Gravitational acceleration

global rhoi = 910.0                 # [kg/m³] Ice density 
global rhow = 1000.0                # [kg/m³] Water density
global rhom = 2700.0                # [kg/m³] Lithosphere density 

global grad = 0.0065                # [K/m] Atmospheric lapse rate (Γ)
global degK = 273.15                # [K] kelvin to celsius change 
global t_0 = 0.0                    # [ºC] Reference freezing temperature
global t_ref = -18.0                # [ºC] Reference climatic temperature (-18 works fine with defaults values) 
global Lv = 2.5e6                   # [J kg⁻¹] Latent heat of vaporization http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
global Rd = 287.0                   # [J K⁻¹kg⁻¹] Dry-air gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/
global Rv = 461.0                   # [J K⁻¹kg⁻¹] Water-vapor gas constant http://pressbooks-dev.oer.hawaii.edu/atmo/chapter/chapter-4-water-vapor/

# -- changing units to IS
global T_0 = t_0 + degK
global T_ref = t_ref + degK