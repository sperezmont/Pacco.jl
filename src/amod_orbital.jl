# =============================
#     Program: amod_orbital.jl
#     Aim: This program contains functions to calculate orbital boundaries 
# =============================


@doc """
    calc_orb_par: 
"""
function calc_orb_par()
    # Check if is necessary to initialize orbital parameters
    
end

@doc """
    calc_ins: given the day of the year, latitude and time before present (present=1950), 
        calculates the daily insolation values. 
        Adapted from insolation.f90 "calc_insol_day_pt()" by Alexander Robinson and Mahé Perrette 2014
"""
function calc_ins(ctl_o, now_o, par_o)
    if par_o["orb_case"] == "calc"
        calc_orb_par()
    else
        error("Input insolation not implemented yet")
    end



    return new_ins
end

function calc_insol_day_pt(day,lat,time_bp,S0,day_year,fldr) result(insol)
        ! Calculate orbital parameters for current year 
        call calc_orbital_par(time_bp,PER,ECC,XOBCH,TPERI,ZAVEXPE,fldr)

        ! Get hourly zenith angle values for input into daily insol function
        PYTIME = dble(day)/dble(day_max)*2.0_dp*pi
        do h = 1, nh
            PCLOCK = dble(h)/dble(nh)*2.0_dp*pi
            call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                       PCLOCK,PYTIME,PDISSE(h),PZEN1(h),PZEN2(h),PZEN3(h),PRAE)
        end do 

        ! Calculate daily insolation at latitude of interest
        insol = calc_insol_day_internal(lat,PDISSE,PZEN1,PZEN2,S0_value)

        return 

    end function calc_insol_day_pt