# =============================
#     Program: amod_defs.jl
#     Aim: definition of parameters as mutable struct
# =============================

# Help:
#   -- "mutable struct" works similar to fortran classes
#   -- Defining as {T<:Real} make sure that the values are correct
#   -- In this script can be found the struct of parameters and NOW

mutable struct RunSettings{T<:Real}
    time_init :: T       # [yr] Starting time (model years)
    time_end  :: T       # [yr] Ending time (model years)  
    dt        :: T       # [yr] Loop timestep 
    dt_out    :: T       # [yr] Frequency of writing
end

mutable struct InitialConditions{T<:Real}
    H_init      :: T     # [m] Initial condition for ice thickness
    Hsed_init   :: T     # [--] Initial condition for sediments thickness
    T_init      :: T     # [degC] Initial condition for ice temperature 
    A_init      :: T     # [] Initial condition for Glenns law coefficient
end

mutable struct RadiativeParameters{T<:Real}
    ins_min       :: T          # [] Insolation minimum value (June 21 65N)
    ins_max       :: T          # [] Insolation maximum value (June 21 65N)
    ins_prei      :: T          # [] Default preindustrial 
    co2_prei      :: T          # [ppm] Default preindustrial
    T_0           :: T          # [degC] Reference freezing temperature
    T_ref         :: T          # [degC] Reference climatic temperature (-18 works fine with defaults values)  
end

mutable struct OrbitalParameters{T<:Real}
    orbital_mode  :: String          # o, op, oe, pe, ope # Tsl formula 
    P_obl         :: T            # Power of obliquity (normalised to At)
    tau_obl       :: T            # [yr] Obliquity period
    P_pre         :: T            # Power of precession (normalised to At)
    tau_pre       :: T            # [yr] Precession period
    P_exc         :: T            # Power of excentricity (normalised to At)
    tau_exc       :: T            # [yr] Excentricity period
end

mutable struct GeoParameters{T<:Real}
    active_iso    :: Bool         # Switch: include active isostatic rebound?
    B_eq          :: T            # [m] equilibrium altitude of the bedrock
    tau_bed       :: T            # [years] relaxation time for the astenosphere
    T_mantle      :: T            # [degC] Mantle temperature 
    Q_geo         :: T            # [mW/m²] Geothermal heat flux
end

mutable struct IceParameters{T<:Real}
    v_kin         :: T         # [m/yr] kinematic wave velocity (measures the speed of inland propagation of the ice streams for a given stress imbalance) (Payne 2004)
    k_1           :: T         # proportionality between Hsed and U
    f_1           :: T         # [--] fraction of the sediments that is removed beacause of U # Golledge 2013 indicates a typical bed erosion of 10-3 mm/yr (for a speed of ~1 km/yr ==> f1 ~ 1e-6)
    f_2           :: T         # [--] fraction of the surface mass blance (M) that increases the presence of sediments because of weathering (typical denudation rate ~ 10^(-5) m/yr
    fstream_min   :: T         # [--] Mininmal Fraction of the ice sheet considered streaming
    fstream_max   :: T         # [--] Maximal Fraction of the ice sheet considered streaming
    L             :: T         # [m] Aspect ratio of the model representing the necessary scaling for converting dS/dx to H/L (typically 10³ km)
    C_s           :: T         # 1e-7 by default?? #raw sliding parameter (m yr^(-1) Pa^(-2), between 10^(-10) and 10^(-5) in Pollard and deConto 
    k             :: T         # [J s^(-1) m^(-1) K^(-1)] Ice thermal conductivity, EISMINT value
    c             :: T         # [J Kg^(-1) K^(-1)] Ice specific heat capacity, EISMINT value 

    A_m           :: T            # Amplitude of the cosinus for M (surface mass balance) # 0.1 dfault
    A_t           :: T            # Amplitude of the cosinus for T (surface temperatures)
    lambda        :: T            # proportionality between positive temperatures and surface melt (m / yr /K )
end

mutable struct NOW{T<:Real}     # calculated variables for t = time (it is updated with each time step)
    time        :: T       # simulation time
    H           :: T       
    Hsed        :: T       
    T           :: T       
    A           :: T       
    T_sl        :: T      
    TMB         :: T       
    S           :: T       
    B           :: T              
    M           :: T       
    Acc         :: T        
    SMB         :: T       
    U_d         :: T       
    U_b         :: T    
    U           :: T       
    T_surf      :: T 
    tau_b       :: T    
    tau_d       :: T     
    Q_dif       :: T   
    Q_difup     :: T 
    Q_difdown   :: T        
    Q_drag      :: T   
    alpha       :: T   
    Q_adv       :: T               
    fstream     :: T
    fstream_ref :: T
    Hdot        :: T   
    Hseddot     :: T  
    Bdot        :: T   
    Tdot        :: T   
    fstreamdot  :: T
    ins         :: T
    co2         :: T 
end





