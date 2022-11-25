## AMOD (adimensional ice-sheet-sediment model)
by Jorge Alvarez-Solas (Fortran, 2017) and adapted to Julia by Sergio Pérez-Montero (2022)


# Quick-start guide
Clone `amod`
```bash
git clone https://github.com/sperezmont/amod.git
```
Go to the main `amod` directory and configure the model
```bash
julia config/config.jl
```
This will create the directory `amod_env` which is the virtual environment that includes the required Julia dependencies. Once it finishes you can run the model.

# How to run
First, open `julia` REPL using in your terminal
```bash
julia
```
Now, include AMOD module
```julia
julia> include("run_experiment.jl")
```
Once the module `AMOD` is loaded, you can run experiments using
```julia
julia> AMOD.run_amod(output, input)
```
where `output` is your experiment output name and `input` its input parameter file. Note that if you don't provide the function with any argument, AMOD will use `output/test_default` folder and will run with `par/amod_default.jl`.

# How to quick-plotting results
AMOD store the model variables in a netCDF file called `amod.nc` so they can be visualized through `ncview` but in order to further analyze them `PLTAMOD` module is included. First, open `julia` REPL
```bash
julia 
```
Now, include PLTAMOD module
```julia
julia> include("plot_experiment.jl")
```
Once the module `PLTAMOD` is loaded, you can plot experiments using
```julia
julia> PLTAMOD.plot_amod(experiment, variables)
```
This function takes as default arguments `output=test_default` and `vars=["ins_norm", "SMB", "H", "Hsed"]`.

# How to run ensembles
Currently working on that... :)

# Documentation
I am currently working on the documentation of the model.

# Examples
