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
This will create the directory `amod_env` which is the virtual environment that includes the required Julia dependencies. Once it finishes you can run the model. Now create the directory for proxy data using a symbolic link:
```bash
ln -s path/to/proxy/directory/ data/
```
At the moment this directory must have the following directories and files inside (I am working on something generalistic)
```bash
Snyder_2016/snyder_2016.nc
Luthi-etal_2008/luthi-etal_2008.nc
Waelbroeck-etal_2002/waelbroeck-etal_2002.nc
Spratt-Lisiecki_2016/spratt-lisiecki_2016.nc
```

# How to run
First, open `julia` REPL using in your terminal
```bash
julia
```
Include the main script
```julia
include("amod.jl")
```
Now, you can run experiments using
```julia
run_amod(out_name="test_default", par_file="amod_default.jl", par2change=[])
```
where `out_name` is your experiment output name, `par_file` its input parameter file and `par2change` a dictionary with the parameters we want to change (one value for each parameter, for ensembles see below). Note that if you don't provide the function with any argument, AMOD will use `output/test_default` folder and will run with `par/amod_default.jl`.

# How to quick-plotting results
AMOD store the model variables in a netCDF file called `amod.nc` so they can be visualized through `ncview` but in order to further analyze them `libs/plot_lib.jl` is included. First, open `julia` REPL and include `amod.jl`. Once the script is loaded, you can plot experiments using
```julia
plot_amod(experiment="test_default", vars=["ins_norm", "SMB", "H", "Hsed"])
```
Note: the values of the arguments are the ones set as default in the function.

# How to run ensembles
To run ensembles we only need to use its own function in `julia` REPL and wait
```julia
run_amod_ensemble(par2per; out_name="test_default_ens", par_file="amod_default.jl")
```
where `par2per` is an `OrderedDict` (ordered dictionary) with the keys and values of the parameters we want to be exchanged. `output` is the name of the desired output directory and `namelist` the original namelist we want to employ. 

Example:
```julia
new_parameters = OrderedDict("active_iso"=>[true, false],
                             "H_init"=>[0, 1000, 2000],
                             "active_sed"=>[true, false],
                             "Hsed_init"=>[0, 1])
run_amod_ensemble(new_parameters, out_name="test1", par_file="amod_default.jl")
```
This ensemble will run 24 simulations swapping between the values of `new_parameters` and store their results in subdirectories located in `output/test1/`.

# Documentation
I am currently working on the documentation of the model.
