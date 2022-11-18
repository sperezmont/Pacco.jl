# =============================
#     Program: amod_radiative.jl
#     Aim: This program contains functions to calculate orbital parameters
# =============================
using Insolation

@doc """
    calc_artificial_insolation: Compute daily average insolation through different parameterizations
"""
function calc_artificial_insolation(now_o, par_o)
    # define switches (0 or 1 to select contribution)
    if par_o["orb_case"] == "ope"
        wo, wp, we = 1, 1, 1
    elseif par_o["orb_case"] == "op"
        wo, wp, we = 1, 1, 0
    elseif par_o["orb_case"] == "oe"
        wo, wp, we = 1, 0, 1
    elseif par_o["orb_case"] == "pe"
        wo, wp, we = 0, 1, 1
    elseif par_o["orb_case"] == "o"
        wo, wp, we = 1, 0, 0
    elseif par_o["orb_case"] == "p"
        wo, wp, we = 0, 1, 0
    elseif par_o["orb_case"] == "e"
        wo, wp, we = 0, 0, 1
    else
        error("ERROR, orbital parameterization not recognized")
    end

    # Return artificial insolation -- I have to discuss this with jas
    return par_o["ins_ref"] - par_o["A_ins"] * (
        wo * par_o["P_obl"] * cos(2.0 * pi * now_o["time"] / par_o["tau_obl"]) +
        wp * par_o["P_pre"] * cos(2.0 * pi * now_o["time"] / par_o["tau_pre"]) +
        we * par_o["P_exc"] * cos(2.0 * pi * now_o["time"] / par_o["tau_exc"]))

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
    calc_laskar_insolation: Compute daily average insolation given latitude, time of year and orbital parameters
    Function adapted to Julia from The Climate Laboratory (Python Notebooks) by Brian E. J. Rose, University at Albany
                                https://brian-rose.github.io/ClimateLaboratoryBook/home.html
    modifications by Marisa Montoya and Sergio Pérez-Montero 
"""
function calc_laskar_insolation(now_o, par_o)
    # now_o["long_peri"], now_o["obliquity"], now_o["exc"] = orbital_params(now_r["time"]) # ϖ = long_peri, γ = obliquity, e = excentricity -- from Laskar 2004 
    #     # Convert precession angle and latitude to radians
    #     phi = deg2rad(lat)

    #     # omega_long (solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
    #     omega_long = calc_solar_longitude(now_o, par_o)

    #     # Compute declination angle of the sun
    #     delta = arcsin(sin(deg2rad(now_o["obl"])) * sin(omega_long))

    #     # Compute Ho, the hour angle at sunrise / sunset
    #     # -- check for no sunrise or no sunset: Berger 1978 eqn (8),(9)
    #     Ho = xr.where( abs(delta)-pi/2+abs(phi) < 0., # there is sunset/sunrise
    #             arccos(-tan(phi)*tan(delta)),
    #             # otherwise figure out if it's all night or all day
    #             xr.where(phi*delta>0., pi, 0.) )
    #     # this is not really the daily average cosine of the zenith angle...
    #     #  it's the integral from sunrise to sunset of that quantity...
    #     coszen = Ho*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(Ho)
    #     # Compute insolation: Berger 1978 eq (10)
    #     Fsw = S0/pi*( (1+ecc*cos(omega_long -deg2rad(long_peri)))**2 / (1-ecc**2)**2 * coszen)
    #     if not (lat_is_xarray or day_is_xarray):
    #         # Dimensional ordering consistent with previous numpy code
    #         return Fsw.transpose().values
    #     else:
    #         return Fsw
end

@doc """
    calc_insol_day: Compute daily average insolation
"""
function calc_insol_day(now_o, par_o)
    if par_o["ins_case"] == "artificial"
        now_o["long_peri"], now_o["obl"], now_o["exc"] = 0, 0, 0
        return calc_artificial_insolation(now_o, par_o)
    elseif par_o["ins_case"] == "laskar"
        error("ERROR, laskar option not implemented yet")
        return calc_laskar_insolation(now_o, par_o)
    else
        error("ERROR, insolation option not recognized")
    end

end

