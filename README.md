## PACCO (Physical Adimensional Climate-Cryosphere mOdel)
by Jorge Alvarez-Solas (AMOD, Fortran, 2017) and adapted to Julia by Sergio Pérez-Montero (2022)

# Quick-start guide
Clone `pacco`
```bash
git clone https://github.com/sperezmont/pacco.git
```
Select or create your branch (for develomental purpose)
```bash
git checkout branch-name      # for an existing branch
git checkout -b new-branch    # for a new branch
```
Go to the main `pacco` directory and configure the model
```bash
julia config.jl
```
This will configure the virtual environment of the model, It includes the required Julia dependencies. Once it finishes you can run the model. Now download/create the directory for proxy data:

Download from github:
```bash
git clone https://github.com/sperezmont/dataForPACCO.git
mv dataForPACCO data 		# plot functions need data as the directory's name
```

Create a virtual link from your local directory to the folder `data` in `pacco` directory
```bash
ln -s path/to/proxy/directory/ data/
```

# How to run
First, open Julia REPL using in your terminal
```bash
julia
```
Include the main script
```julia
include("Pacco.jl")
```
Now, you can run experiments using
```julia
run_pacco(experiment; p)
```
Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions

# How to quick-plotting results
PACCO store the model variables in a netCDF file called `pacco.nc` so they can be visualized through `ncview` but in order to further analyze them `libs/plot_lib.jl` is included. First, open `julia` REPL and include `Pacco.jl`. Once the script is loaded, you can plot experiments using
```julia
plot_pacco(experiment/experiments, vars2plot)  
```
Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions
# How to run ensembles
To run ensembles we only need to use its own functions in `julia` REPL and wait

If we want to use a Latin Hypercube Sampling method
```julia
run_pacco_lhs(experiment, params2per, nsim)
```
or if we want to use strictly certain values
```julia
run_pacco_ensemble(experiment, params2per) 
```
Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions

# Documentation
I am currently working on the documentation of the model.

# Some additional notes
PACCO parameters are defined as a Julia (non-mutable) `struct` which means it cannot be overdefined after defined. The default structure is stored in `par/` directory under the name `default_params.jl`. Thus, if we want to work with a certain set of parameters we have two options:
* First option: Open `Pacco.jl` and edit the 5th line and change the value of `source_parameter_file` to the filename we want to use (previously we should create the file we want to use copying the default `struct` file).
```bash
vim Pacco.jl    # edit 5th line with for example "example_params.jl"
julia
```
```julia
include("Pacco.jl")   # now we run simulations using as default parameters the ones from "example_params.jl"
```
* Second option: Go to `par/` directory and save the default `struct` file as `default_params_orig.jl`. Then, rename the file we want to use as `default_params.jl`
```bash
cd par/
mv default_params.jl default_params_orig.jl
mv example_params.jl default_params.jl
cd ..
julia
```julia
include("Pacco.jl")   # now it includes the new struct loading the same parameter file
```
Note: However you can always use the original `default_params.jl` and change manually some parameter when running `run_pacco("test1", p = Params(new_values))`. These previous instructions are intended for changes numerous parameters.
