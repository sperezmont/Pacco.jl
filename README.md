## PACCO (Paleo Adimensional Cryosphere Climate mOdel)
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
julia config/config.jl
```
This will create the directory `env` which is the virtual environment that includes the required Julia dependencies. Once it finishes you can run the model. Now download/create the directory for proxy data:

```bash
git clone https://github.com/sperezmont/Proxy.git
mv Proxy data
```

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
run_pacco(;experiment, par_file, par2change)
```
Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions

# How to quick-plotting results
PACCO store the model variables in a netCDF file called `pacco.nc` so they can be visualized through `ncview` but in order to further analyze them `libs/plot_lib.jl` is included. First, open `julia` REPL and include `Pacco.jl`. Once the script is loaded, you can plot experiments using
```julia
plot_pacco(;experiment/experiments, vars2plot)  
```
Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions
# How to run ensembles
To run ensembles we only need to use its own functions in `julia` REPL and wait

If we want to use a Latin Hypercube Sampling method (preferred methodoly)
```julia
run_pacco_lhs(par2per, nsim; experiment, par_file)
```
or if we want to use strictly certain values
```julia
run_ensemble(par2per; experiment, par_file) 
```
Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions

Example:
```julia
par2per = Dict("ki"=>(0.001, 0.09), "lambda"=>(0.01, 0.1))
nsims = 30
exp2do = "ice-clim_sed0_ki"
parfile = "pacco_ice-climate_def.jl"

run_pacco_lhs(par2per, nsims; experiment=exp2do, par_file=parfile)
```
This example should provide 30 simulations combining random values of parameters `ki` and `lambda` in the selected range. 

```julia
new_parameters = OrderedDict("active_iso"=>[true, false],
                             "H_init"=>[0.0, 1000.0, 2000.0],
                             "active_sed"=>[true, false],
                             "Hsed_init"=>[0.0, 1.0])
run_pacco_ensemble(new_parameters, out_name="test1", par_file="pacco_default.jl")
```
This ensemble will run 24 simulations swapping between the values of `new_parameters` and store their results in subdirectories located in `output/test1/`.

Note: remember to use `?function` in Julia REPL if you need to know how to use model's functions

# Documentation
I am currently working on the documentation of the model.
