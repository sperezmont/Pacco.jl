## PACCO (Paleo Adimensional Climate-Cryosphere mOdel)
by Jorge Alvarez-Solas (AMOD, Fortran, 2017) and adapted to Julia by Sergio Pérez-Montero (2022)

![Model scheme](config/pacco_scheme.png?raw=true)

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

If we want to use a Latin Hypercube Sampling method
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

# How to analyze ensembles
In order to analyze large ensembles `analyze_pacco()` is provided:
```julia
analyze_runs(;experiment="test_ensemble_default", experiments=[], rule::Integer=1, reanalyze=[])
```
This function takes `experiment` or a list `experiments` and applies a rule `rule` (number) to them. In order to see what rules are defined, use:
```julia
see_rules()
```
Applying this function to an ensemble will generate a file called `good_runs_rx.txt` which contains the names of the runs that fullfill the conditions applied. `x` is the number of the rule applied. If we use the optional keyword `reanalyze` with the name of a previous txt file it will analyze the previously analyzed runs. Thus, the working flow should be:
```julia
new_parameters = OrderedDict("active_iso"=>[true, false],
                             "H_init"=>[0.0, 1000.0, 2000.0],
                             "active_sed"=>[true, false],
                             "Hsed_init"=>[0.0, 1.0])
run_pacco_ensemble(new_parameters, out_name="test1", par_file="pacco_default.jl")

analyze_pacco(experiment="test1", rule=1) # this rule only considers runs with at least 1 complete deglaciation and generates "good_runs_r1.txt"
analyze_pacco(experiment="test1", rule=3, reanalyze="good_runs_r1.txt") # this line will apply rule 3 to the previously analyzed set stored in "good_runs_r1.txt"
```
Once we have analyzed the ensemble, we can plot results with `fast_plot()`, `fast_histogram()` and `plot_std()`
```julia
fast_plot(experiment, filename; var2plot="H_n", cmap=:darkrainbow, all_runs=true, plot_PSD=false)
fast_histogram(experiment, filename; all_runs=true)
```
`fast_plot()` and `fast_histogram()` require the name of the experiment and the name of the text file where the "good runs" are stored (these arguments do not have default values). The first function plots the variable `var2plot` and can be used to just plot the Power Spectrum Density diagram setting `plot_PSD=true`, and the second function plots a histogram with the amount of simulations as a function of the parameter value employed (It also plots points of the simulation number and its value, with a color code given by `cmap`). In both functions we can plot in two modes: `all_runs=true` and `all_runs=false`. If this argument is set to `true` the function will take into account the total numnber of simulations when applying the colors (this will help to see the clustering of valid runs), and if it is set to `false` the function just uses the number of "good runs". 
Note: I use `fast_histogram()` in order to see the parameters values that have a major amount of valid runs.

```julia
plot_std(experiment; var2plot="H_n", cmap=:darktest)
```
This function plots the standard deviation of each run in the **entire ensemble**. Note: Useful to see how the clouster of permutations increases or decreases the variability in the run.  

Note: `plot_pacco()` can also be employed, but **this function does not apply any filter**, so it will plot every run in the ensemble.

# Documentation
I am currently working on the documentation of the model.
