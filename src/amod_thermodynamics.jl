# =============================
#     Program: amod_thermodynamics.jl
#     Aim: This program contains functions to calculate thermodynamics
# =============================

@doc """
    calc_T_surf: calculates surface temperature
"""
function calc_T_surf(now_t, par_t)
    if par_t["tsurf_case"] == "linear"
        return now_t["T_sl"] - grad * now_t["S"]
    else
        write(f, "ERROR, T_surf option not recognized")
        return now_t["T_surf"]
    end
end

@doc """
    calc_snowfall: calculates snowfall rate
"""
function calc_snowfall(now_t, par_t)
    t = now_t["T_surf"] - degK    # conversion to ºC

    # First, calculate the saturation vapor pressure e_s
    if par_t["cc_case"] == "AERKi"
        e_s = 6.1121 * exp((22.587 * t) / (t + 273.86))     # Alduchov and Eskridge (1996)
    else
        write(f, "ERROR, e_s option not recognized")
    end

    # Now, compute the amount of water


    # Then, calculate the percentage of snow in precipitation

end

@doc """
    calc_SMB: calculates surface mass balance
"""
function calc_SMB(now_t, par_t)

    # First, calculates Accumulation
    now["Acc"] = calc_snowfall(now_t, par_t) 


    now%M = par%Am*cos(2.0*PI*now%time/40000.0)
    now%Acc = par%Am*(1.1+cos(2.0*PI*now%time/40000.0))+par%Am*0.1*cos(2.0*PI*now%time/100000.0)

    # DEFINIR DE TAL FORMA QUE ACC DEPENDA DE LA T_SL
    # 10 es un offset para dejar fundir entre -10 y 0
    if (now%Tsurf.ge.-10.0) then          ! assuming there is some melt when mean annual T > -5º C 
        now%M = par%lambda*(now%Tsurf+10.0)     ! The 5.0 offset calibrates the Melt to be 0 for a mean annual T of -5ºC  
    else
        now%M = 0.0 
    endif
    
            now%SMB = now%Acc - now%M
end