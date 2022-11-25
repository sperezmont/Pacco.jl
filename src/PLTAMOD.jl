module PLTAMOD
    """
    Module PLOT_AMOD contains all the dependencies and functions to plot AMOD model from REPL
    """
    # Import dependencies
    using NCDatasets, DSP, Statistics

    # Determine amod path
    global amod_path = pwd()    # -- this line assumes we are working from the main AMOD directory

    # Include libraries and functions
    include("../libs/nc.jl")
    include("../libs/plot_lib.jl")

end