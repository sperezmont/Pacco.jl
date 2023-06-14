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
    calc_laskar_insolation(t, lat=65.0, day=170.0, S0=1365.2, day_type=1, days_per_year=365.2422)
Compute daily average insolation given latitude and day of the year
Adapted from The Climate Laboratory of Brian E. J. Rose (https://brian-rose.github.io/ClimateLaboratoryBook/home.html)
"""
function calc_laskar_insolation(t::Real; lat::Real=65.0, day::Real=170.0, S0::Real=1365.2, day_type::Real=1, days_per_year::Real=365.2422)
    # First, calculate the orbital parameters at t (years) since J2000 epoch
    long_peri, obliquity, ecc = Insolation.orbital_params(OrbData, t) # -- using Insolation.jl (rad, rad, --)
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
    calc_summer_insolation(t, τ; lat=65.0, S0=1365.2, days_per_year=365.2422)
calculates the integrated summer insolation (ISI) as defined in Leloup and Paillard (2022) after Huybers (2006)
                    ISI(τ) = sum(βi .* (Wi * 86400))
where τ is an insolation threshold that defines the summer days. It is selected as 275 W/m² in Huybers (2006) in order to produce an analogue of PDD for insolation,
T must be ≥ 0 or it does not account, so in insolation, between 40º and 70ºN T == 0ºC if insolation = [250, 300] W/m². In Leloup and Paillard (2022) 300 and 400 W/m² are analyzed.
"""
function calc_summer_insolation(t::Real, τ::Real; lat::Real=65.0, S0::Real=1365.2, days_per_year::Real=365.2422)
    ISI, days_above_τ = 0.0, 0
    for i in 1:1:365
        Wi = calc_laskar_insolation(t, lat=lat, day=i, S0=S0, day_type=1, days_per_year=days_per_year)
        if Wi > τ
            βi = 1.0
            days_above_τ += 1
        else
            βi = 0.0
        end
        ISI += βi * (Wi * 86400) # J/m²
    end

    return ISI / (days_above_τ * 86400) # redimensionalization to W/m²
end

"""
    calc_caloric_insolation(t, lat=65.0, S0=1365.2, days_per_year=365.2422)
calculates the caloric seasons insolation following  Tzedakis et al. (2017) and Milankovitch (1941)
"""
function calc_caloric_insolation(t::Real; lat::Real=65.0, S0::Real=1365.2, days_per_year::Real=365.2422)
    insolations = Vector{Any}(undef, 365)
    half_year = Int(ceil(days_per_year / 2))
    for i in 1:1:365
        insolations[i] = calc_laskar_insolation(t, lat=lat, day=i, S0=S0, day_type=1, days_per_year=days_per_year)
    end
    sort!(insolations, rev=true)
    return sum(insolations[1:half_year]) / half_year # we only integrate over the half of the days with the higher values
end

"""
    read_insolation_from_file(filename, tspan)
reads a .jld2 file that contains a matrix (time/insolation, length) of precalculated insolation and its time vector

## Attributes
* `filename` is the path to input/ and the file name we want to load
* `tspan` is a tuple with the complete timespan of the simulation
"""
function read_insolation_from_file(filename::String, tspan::Tuple)
    data_matrix = JLD2.load_object(filename)
    index_start = findmin(abs.(data_matrix[1, :] .- tspan[1]))
    index_end = findmin(abs.(data_matrix[1, :] .- tspan[2]))

    if (tspan[1] < data_matrix[1, 1]) || (tspan[2] > data_matrix[1, end]) 
        error("The selected input file does not match the time span required")
    end
    return data_matrix[:, index_start[2]:index_end[2]]
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
    elseif p.I_case == "laskar" # this option should be improved using Insolation.jl alone
        u[10] = calc_laskar_insolation(t,
            lat=p.I_lat,
            day=p.I_day,
            S0=p.I0,
            day_type=1,
            days_per_year=p.year_len)
    elseif p.I_case == "summer"
        u[10] = calc_summer_insolation(t, p.Ithreshold,
        lat=p.I_lat,
        S0=p.I0,
        days_per_year=p.year_len)
    elseif p.I_case == "caloric"
        u[10] = calc_caloric_insolation(t,
        lat=p.I_lat,
        S0=p.I0,
        days_per_year=p.year_len)
    elseif p.I_case == "input"
        u[10] = InsolationData[2, InsolationData[1, :] .== t][1]
    else
        error("ERROR, insolation option not recognized")
    end
    return nothing
end