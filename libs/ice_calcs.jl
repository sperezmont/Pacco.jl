# =============================
#     Program: ice_calcs.jl
#     Aim: This program contains functions to calculate convertions and other stuff
# =============================

# -- IME 
@doc """ 
        IME: Transforms volume data to meters of ice equivalent
            [data] = km^3
    """
function IME(data; rhoi=0.9167, rhow=1.0)
    ime = 

end

# -- SLE
@doc """ 
        SLE: Transforms volume data to m SLE
            [data] = km^3 
            rhow = 1 Gt/m3 -> "disregarding the minor salinity/density effects of mixing fresh meltwater with seawater"
                More about this: https://sealevel.info/conversion_factors.html
    """
function SLE(data; rhoi=0.9167, rhow=1.0, Aoc=3.618*10^8)
    SLE = rhoi/rhow .* 1e3 ./ Aoc .* data
    return SLE
end

#-- SLR
@doc """ 
        SLR: Transforms volume data to m SLE and calculates the SLR for each time step
            [data] = km**3
            rhow = 1 Gt/m3 -> "disregarding the minor salinity/density effects of mixing fresh meltwater with seawater"
                More about: https://sealevel.info/conversion_factors.html 
    """
function SLR(data; rhoi=0.9167, rhow=1, Aoc=3.618*10^8, v2s=false)
    (v2s) ? sle = SLE(data, rhoi=rhoi, rhow=rhow, Aoc=Aoc) : sle = copy(data)
    if ndims(data) == 1             # 1 time series    
        SLR = sle[1] .- sle
    elseif ndims(data) == 2         # >1 time series
        SLR = sle[:, 1] .- sle
    elseif ndims(data) == 3         # 1 2D field
        SLR = sle[1, :, :] .- sle
    elseif ndims(data) == 4         # >1 2D field
        SLR = sle[:, 1, :, :] .- sle
    end
    return SLR
end

