module AMOD
    """
        Module AMOD contains all the dependencies and functions to run AMOD model from REPL
    """
    # Import dependencies
    using Pkg, NCDatasets, DataStructures

    # Determine amod path
    global amod_path = pwd()    # -- this line assumes we are working from the main AMOD directory

    # Include libraries and functions
    # -- libraries
    include("../libs/load_args.jl")
    include("../libs/nc.jl")

    # -- earth constants
    include("../par/earth_const.jl")

    # -- amod functions
    include("./lib_amod.jl")
    include("./amod_defs.jl")
    include("./amod_update.jl")
    include("./amod_radiative.jl")
    include("./amod_orbital.jl")
    include("./amod_dynamics.jl")
    include("./amod_thermodynamics.jl")
end



