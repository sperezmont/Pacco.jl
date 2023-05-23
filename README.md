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
git clone https://github.com/sperezmont/Proxy.git
mv Proxy data 		# plot functions need data as the directory's name
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
