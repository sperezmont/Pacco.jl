## AMOD
# Quick-start guide
Clone `amod_tools`
```bash
git clone https://github.com/sperezmont/amod.git
```
Go to the main `amod` directory and configure the model
```bash
julia config/config.jl
```
This will install the necessary Julia dependencies that needs the model. Once it finishes you can run `amod`

# How to run
```bash
julia amod.jl [OUT] [PARF]
```
where `[OUT]` is the desired output directory name located in `output/` and `[PARF]` is the desired parameters file name located in `par/` directory.

Note: If you just run
  ```bash
  julia amod.jl
  ```
  or
  ```bash
  julia amod.jl [OUT]
  ```
  the model will generate a directory called `test_default` (or `[OUT]`) and will use the parameter file called `amod_default.jl`.

# How to quick-plotting results
`amod` store the model variables in a netCDF file called `amod.nc` so they can be visualized through `ncview` but in order to further analyze them `plot_amod.jl` is included. This script will plot the variables passed as arguments and their frequencies spectrum. `plot_amod.jl` is executed as follows:
```bash
julia plot_amod.jl [OUT] [VARS]
```
where `[OUT]` is the name of the experiment and `[VARS]` are the variables we want to plot separated by spaces. If none of the arguments are provided, `plot_amod.jl` will plot the results of `test_default` and some default variables. 

# Examples
