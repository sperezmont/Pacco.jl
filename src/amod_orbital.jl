# =============================
#     Program: amod_radiative.jl
#     Aim: This program contains functions to calculate orbital parameters
# =============================
@doc """
    calc_artificial_insolation: Compute daily average insolation through different parameterizations
"""
function calc_artificial_insolation(now_o, par_o)
    # Calculate reference and amplitude
    ins_ref = (par_o["ins_max"] + par_o["ins_min"]) / 2
    A_ins = (par_o["ins_max"] - par_o["ins_min"]) / 2

    # Return artificial insolation -- I have to discuss this with jas
    ins = ins_ref + A_ins * (
        par_o["P_obl"] * cos(2.0 * pi * now_o["time"] / par_o["tau_obl"]) +
        par_o["P_pre"] * cos(2.0 * pi * now_o["time"] / par_o["tau_pre"]) +
        par_o["P_exc"] * cos(2.0 * pi * now_o["time"] / par_o["tau_exc"]))
    for hm in hemisphere
        now_o["ins"*hm] = ins
    end
    return now_o
end

@doc """
    calc_solar_longitude: Estimates solar longitude from calendar day
        Function adapted to Julia from The Climate Laboratory (Python Notebooks) by Brian E. J. Rose, University at Albany
                                    https://brian-rose.github.io/ClimateLaboratoryBook/home.html
        modifications by Marisa Montoya and Sergio Pérez-Montero 
    
        Notes: 
            * Method is using an approximation from Berger et al., 1978 section 3 (omega = 0 at spring equinox)
            * The calendar is referenced to the vernal equinox which always occurs at day 80
            * Solar longitude is the angle of the Earth's orbit measured from spring equinox (21 March)
            * Note that calendar days and solar longitude are not linearly related because, by Kepler's Second Law, Earth's
              angular velocity varies according to its distance from the sun.

"""
function calc_solar_longitude(now_o, par_o)
    # Convert to radians
    long_peri_rad = deg2rad(now_o["long_peri"])

    # Calculate the distance to vernal equinox
    delta_omega = (day - 80.0) * 2 * pi / year_len

    # Calculate beta
    beta = sqrt(1 - now_o["exc"]^2)

    # Compute solar longitude
    omega_long_m = -2 * ((now_o["exc"] / 2 + (now_o["exc"]^3) / 8) * (1 + beta) * sin(-long_peri_rad) -
                         (now_o["exc"]^2) / 4 * (1 / 2 + beta) * sin(-2 * long_peri_rad) + (now_o["exc"]^3) / 8 *
                                                                                           (1 / 3 + beta) * sin(-3 * long_peri_rad)) + delta_omega
    omega_long = (omega_long_m + (2 * now_o["exc"] - (now_o["exc"]^3) / 4) * sin(omega_long_m - long_peri_rad) +
                  (5 / 4) * (now_o["exc"]^2) * sin(2 * (omega_long_m - long_peri_rad)) + (13 / 12) * (now_o["exc"]^3)
                                                                                         * sin(3 * (omega_long_m - long_peri_rad)))
    return omega_long
end

@doc """
    calc_laskar_insolation: Compute daily average insolation given latitude and time of year
"""
function calc_laskar_insolation(now_o, par_o)
    return now_o
end

@doc """
    calc_insol_day: Compute daily average insolation
"""
function calc_insol_day(now_o, par_o)
    if par_o["ins_case"] == "artificial"
        now_o = calc_artificial_insolation(now_o, par_o)
    elseif par_o["ins_case"] == "laskar"
        error("ERROR, laskar option not implemented yet")
        now_o = calc_laskar_insolation(now_o, par_o)        
    else
        error("ERROR, insolation option not recognized")
    end
    return now_o
end

