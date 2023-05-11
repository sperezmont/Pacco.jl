# =============================
#     Program: orbital.jl
#     Aim: This program contains functions to calculate orbital parameters
# =============================

"""
    calc_artificial_insolation(p, t)
Compute daily average insolation through different parameterizations
"""
function calc_artificial_insolation(p::Params, t::Real)
    # Calculate reference and amplitude
    I_ref = (p.I_max + p.I_min) / 2
    A_ins = (p.I_max - p.I_min) / 2

    # Return artificial insolation -- I have to discuss this with jas
    return I_ref + A_ins * (
        p.P_obl * cos(2.0 * pi * t / p.tau_obl) +
        p.P_pre * cos(2.0 * pi * t / p.tau_pre) +
        p.P_exc * cos(2.0 * pi * t / p.tau_exc))
end

"""
    calc_solar_longitude(dat, long_peri, ecc)
Compute solar longitude given day of the year and orbital parameters
Adapted from The Climate Laboratory of Brian E. J. Rose (https://brian-rose.github.io/ClimateLaboratoryBook/home.html)
"""
function calc_solar_longitude(day::Real, long_peri::Real, ecc::Real; days_per_year::Real=365.2422)
    long_peri_rad = deg2rad(long_peri)
    delta_omega = (day - 80.0) * 2 * pi / days_per_year
    beta = sqrt(1 - ecc^2)

    omega_long_m = -2 * ((ecc / 2 + (ecc^3) / 8) * (1 + beta) * sin(-long_peri_rad) -
                         (ecc^2) / 4 * (1 / 2 + beta) * sin(-2 * long_peri_rad) +
                         (ecc^3) / 8 * (1 / 3 + beta) * sin(-3 * long_peri_rad)) + delta_omega

    omega_long = omega_long_m + (2 * ecc - (ecc^3) / 4) * sin(omega_long_m - long_peri_rad) +
                 (5 / 4) * (ecc^2) * sin(2 * (omega_long_m - long_peri_rad)) +
                 (13 / 12) * (ecc^3) * sin(3 * (omega_long_m - long_peri_rad))
    return omega_long
end

"""
    calc_laskar_insolation(t)
Compute daily average insolation given latitude and day of the year
Adapted from The Climate Laboratory of Brian E. J. Rose (https://brian-rose.github.io/ClimateLaboratoryBook/home.html)
"""
function calc_laskar_insolation(t::Real; lat::Real=65.0, day::Real=170.0, S0::Real=1365.2, day_type::Real=1, days_per_year::Real=365.2422)
    # First, calculate the orbital parameters at t (years) since J2000 epoch
    long_peri, obliquity, ecc = orbital_params(t) # -- using Insolation.jl (rad, rad, --)
    long_peri, obliquity = rad2deg(long_peri), rad2deg(obliquity)
    phi = deg2rad(lat)

    # -- calculate solar longitude
    if day_type == 1
        omega_long = calc_solar_longitude(day, long_peri, ecc, days_per_year=days_per_year)
    elseif day_type == 2
        omega_long = deg2rad(day)
    end

    # -- declination angle of the Sun
    delta = asin(sin(deg2rad(obliquity)) * sin(omega_long))

    # -- hour angle at sunrise / sunset
    if abs(delta) - pi / 2 + abs(phi) < 0.0 # there is sunset/sunrise
        Ho = acos(-tan(phi) * tan(delta))
    else # otherwise figure out if it's all night or all day
        if phi * delta > 0.0
            Ho = pi
        else
            Ho = 0.0
        end
    end

    # -- cos zenith angle
    coszen = Ho * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(Ho)

    # Compute insolation as in Berger 1978 (equation 10)
    Fsw = S0 / pi * ((1 + ecc * cos(omega_long - deg2rad(long_peri)))^2 / (1 - ecc^2)^2 * coszen)
    return Fsw
end

"""
    calc_I!(u, p, t)
Compute daily average insolation
"""
function calc_I!(u::Vector, p::Params, t::Real)
    if p.I_case == "constant"
        u[10] = p.I_const
    elseif p.I_case == "artificial"
        u[10] = calc_artificial_insolation(p, t)
    elseif p.I_case == "laskar"
        u[10] = calc_laskar_insolation(t,
            lat=p.I_lat,
            day=p.I_day,
            S0=p.I0,
            day_type=1,
            days_per_year=p.year_len)
    else
        error("ERROR, insolation option not recognized")
    end
    return nothing
end

