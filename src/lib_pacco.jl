# =============================
#     Program: lib_pacco.jl
#     Aim: main functions of PACCO
# =============================
"""
    calc_diagnostic_variables!(u, dudt, p, t)
computes diagnostic variables in `u` for time `t` using parameters in `p`
"""
function calc_diagnostic_variables!(u::Vector, p::Params, t::Real)
    # Compute Orbital Forcing
    calc_insolation!(u, p, t)

    # Translate Orbital Forcing to temperature
    if p.active_climate == false
        calc_sealevel_temperature!(u, p, t)
    end
    calc_reference_temperature!(u, p, t)

    # Compute Ice-Sheet Geometry
    calc_icesheet_elevation!(u, p)
    calc_icesheet_area!(u, p)
    calc_icesheet_volume!(u, p)

    if p.active_climate
        # Compute Climate variables
        calc_reference_albedo!(u, p)
    end
    calc_surface_temperature!(u, p)

    # Compute Ice-Sheet Mass Balance
    calc_snowfall!(u, p)
    calc_ablation!(u, p)

    if p.active_ice
        # Compute Ice-Sheet dynamics
        calc_driving_stress!(u, p)
        calc_basal_stress!(u, p)
        calc_deformational_velocity!(u, p)
        calc_basal_velocity!(u, p)

        (p.active_thermo) && (calc_reference_streaming_fraction!(u, p)) # update reference streaming fraction

        u[25] = u[22] + u[23]    # calculates total velocity in the ice sheet

        if p.active_thermo
            # Compute Ice-Sheet thermodynamics
            calc_conductive_heat!(u, p)
            calc_dragging_heat!(u, p)
            calc_geothermal_heat!(u, p)
            calc_vertical_heat_advection!(u, p)
        end
    end

    return nothing
end

"""
    dudt!(dudt, u, p, t)
computes derivatives of prognostic variables in `u` using parameters in `p` at time `t`
"""
function dudt!(dudt::Vector, u::Vector, p::Params, t::Real)
    ###########################
    # Update diagnostics
    ###########################
    calc_diagnostic_variables!(u, p, t)

    ###########################
    # Update prognostics
    ###########################
    if p.active_climate
        # -- regional air temperature
        dudt[1] = calcdot_regional_temperature(u, p, t)

        # -- C
        dudt[2] = calcdot_carbon_dioxide(u, p, t)

        # -- ice age
        dudt[3] = calcdot_iceage()

        # -- albedo
        dudt[4] = calcdot_albedo(u, p)
    else
        dudt[1:4] .= 0.0
    end

    # -- ice thickness
    dudt[5] = calcdot_icethickness(u, p)

    # -- sediments layer thickness
    if t < p.time_init
        dudt[6] = 0.0
    else
        (p.active_sed) ? (dudt[6] = calcdot_sediment_thickness(u, p)) : (dudt[6] = 0.0)
    end

    # -- bedrock elevation
    (p.active_iso) ? (dudt[7] = calcdot_bedrock_elevation(u, p)) : (dudt[7] = 0.0)

    if p.active_thermo
        # -- ice temperature
        dudt[8] = calcdot_ice_temperature(u, p)

        # -- streaming fraction
        dudt[9] = calcdot_streaming_fraction(u, p)
    else
        dudt[8:9] .= 0.0
    end

    # -- diagnostic derivatives
    view(dudt, lprog+1:lsu) .= 0.0

    # Modify states to ensure physical meaning
    view(u, 1:lprog) .= max.(view(u, 1:lprog), [0.0, 1.0, 0.0, p.albedo_land, 0.0, 0.0, -Inf, 0.0, 0.0])
    view(u, 1:lprog) .= min.(view(u, 1:lprog), [Inf, Inf, Inf, Inf, Inf, Inf, Inf, p.Tmp, Inf])

    if p.regtemp_case == "comp" # two components in the signal
        u[1] += p.kT * (t - p.time_init)
    end

    if p.carbon_case == "comp" # two components in the signal
        u[2] += p.kC * (t - p.time_init)
    end

    if u[5] == p.ice_exists_thr              # (m) no ice
        u[3] = 0.0              # ice age is set to 0
        u[4] = p.albedo_land    # albedo is land albedo
        u[9] = p.fstrmin        # streaming is set to minimum
        u[16] = p.albedo_land   # reference albedo is land albedo

    elseif  u[5] < p.ice_is_big_thr         # (m) the ice sheet is small
        u[3] = 0.0              # ice age is set to 0
        u[4] = p.albedo_newice  # albedo is new ice        
        u[9] = p.fstrmin        # streaming is set to minimum
        u[16] = p.albedo_newice # to stop evolution

    elseif (u[3] < p.ice_is_old_thr)        # (yr) the ice sheet is young 
        u[4] = p.albedo_newice  # albedo is new ice
    end

    return nothing
end

"""
    pacco(u0, p, tspan)
takes initial conditions for `u0`, parameters in `p` and solves PACCO for timespan `tspan`. 
This function uses `BS3()` algorithm to solve the problem and saves results each `p.dt_out` time step. 
This option was selected because of the model's simplicity and the computational speed of the method as recomended by DifferentialEquations.jl documentation:
"For fast solving at higher tolerances, we recommend BS3". This option uses the Bogacki-Shampine 3/2 method:
* https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/ 
* https://en.wikipedia.org/wiki/Bogacki%E2%80%93Shampine_method
"""
function pacco(u0::Vector, p::Params, tspan::Tuple)
    # Define and solve the problem
    prob = ODEProblem(dudt!, u0, tspan, p)

    if (p.dt_case == "fixed") || (p.insol_case == "input")
        return solve(prob, Euler(), saveat=p.dt_out, dt=p.dt)
    elseif p.dt_case == "adaptive"
        return solve(prob, BS3(), saveat=p.dt_out)
    else
        error("`dt_case = $(p.dt_case)` option not recognized")
    end
end

"""
    run_pacco(experiment; p)
runs experiment `experiment` using as source parameter file `nml_ice-clim.jl` and modify the parameters selected in `pars`

## Attributes
* `experiment` String with the path to experiment from `output/` folder (default is `test_default`)
* `p` `Params()` struct filled with the new values of the parameters we want to change
## Examples:
        run_pacco("test1")
Runs experiment `test1` using the default values of `p` located in `par/default_params.jl`

        run_pacco("test1", p = Params(time_init = -2e6, sref = 0.1, lambda = 0.1))
Runs experiment `test1` using the default values of `p` located in `par/default_params.jl` and change those provided by `p = Params(time_init = -2e6, sref = 0.1, lambda = 0.1)`

        run_pacco("test1", p = JLD2.load_object("path/to/julia/object/params.jld2"))
Runs experiment `test1` and uses as parameters the ones stored in `"path/to/julia/object/params.jld2"`
"""
function run_pacco(experiment::String; p::Params=Params(), returnsol=false)
    ## Now, load arguments
    output_path = pwd() * "/output/" * experiment * "/"
    isdir(output_path) || mkdir(output_path)
    write_run_info(output_path, experiment, p)
    JLD2.save_object(output_path * "/params.jld2", p)
    timespan = (p.time_init - p.time_spinup, p.time_end)

    ## Load input if desired
    if p.insol_case == "input"
        global InsolationData = read_insolation_from_file(p.insol_input, timespan)    # defined as global in order to be accesible from the entire model
    end

    ## Run pacco()
    u0, out_attr = load_defs(p)
    out = pacco(u0, p, timespan)

    ## Create output nc file
    # if outfile exists remove it
    isfile(output_path * "/pacco.nc") && rm(output_path * "/pacco.nc")
    genout_nc(output_path, "/pacco.nc", out, p, out_attr)

    ## Done!
    if returnsol
        return out
    end
end

"""
    run_pacco_ensemble(experiment, params2per)
run the ensemble `experiment` given by the permutations of values in dictionary `params2per`
"""
function run_pacco_ensemble(experiment::String, params2per::Dict)
    # First, obtain simulations
    perm = calc_permutations(params2per)

    # Second, get number of simulations
    nperms = length(perm)

    # Now, create ensemble directory
    if isdir(pwd() * "/output/" * experiment)
        rm(pwd() * "/output/" * experiment, recursive=true)
    end
    mkdir(pwd() * "/output/" * experiment)
    mkdir(pwd() * "/output/" * experiment * "/runs/")
    mkdir(pwd() * "/output/" * experiment * "/results")

    # Save combinations and names in permutations.txt
    file_perm = open(pwd() * "/output/" * experiment * "/results/permutations.txt", "w")

    # Third, run each permutation
    printstyled("Running $(experiment): $(nperms) runs\n", color=:green)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in 1:nperms
        experimenti = "/runs/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        line2print = "run $(i)/$(nperms) s$(repeat("0", ndigits - length(digits(i))))$i $(perm[i])"
        write(file_perm, "$(line2print) \n")

        new_pi = Params(; (Symbol(k) => v for (k, v) in perm[i])...)
        run_pacco(experiment * experimenti, p=new_pi)
        println("run $(i)/$(nperms) s$(repeat("0", ndigits - length(digits(i))))$i")
    end
    close(file_perm)

end

"""
    run_pacco_lhs(experiment, params2per, nsim)
run the ensemble `experiment` given by the permutations of values in dictionary `params2per` using Latin Hypercube Sampling and subdividing the given ranges in `nsim` elements
"""
function run_pacco_lhs(experiment::String, params2per::Dict, nsim::Int)
    parnames = collect(keys(params2per))
    # First, calculate LHS
    permutations, permutations_dict = gen_lhs(params2per, nsim; pars_type=2)

    # Now, create ensemble directory
    if isdir(pwd() * "/output/" * experiment)
        rm(pwd() * "/output/" * experiment, recursive=true)
    end
    mkdir(pwd() * "/output/" * experiment)
    mkdir(pwd() * "/output/" * experiment * "/runs/")

    # Save combinations and names in permutations.txt
    file_perm = open(pwd() * "/output/" * experiment * "/results/permutations.txt", "w")

    # Plot (if possible) the LHS
    mkdir(pwd() * "/output/" * experiment * "/results")
    if length(params2per) == 2
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2])
        scatter!(ax, permutations)
        save(pwd() * "/output/" * experiment * "/results/lhs.png", fig)
    elseif length(params2per) == 3
        fig = Figure()
        ax = Axis3(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2], zlabel=parnames[3])
        scatter!(ax, permutations)
        save(pwd() * "/output/" * experiment * "/results/lhs.png", fig)
    end

    # Run each permutation
    nperms = length(permutations_dict)
    printstyled("Running $(experiment): $(nperms) runs\n", color=:green)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in 1:nperms
        experimenti = "/runs/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        line2print = "s" * repeat("0", ndigits - length(digits(i))) * "$i $(permutations_dict[i])"
        line2print = "run $(i)/$(nperms) s$(repeat("0", ndigits - length(digits(i))))$i $(permutations_dict[i])"
        write(file_perm, "$(line2print) \n")

        new_pi = Params(; (Symbol(k) => v for (k, v) in permutations_dict[i])...)
        run_pacco(experiment * experimenti, p=new_pi)
        println("run $(i)/$(nperms) s$(repeat("0", ndigits - length(digits(i))))$i")
    end
    # Done!
end

"""
    run_fixed_ensemble(experiment, vop)
creates the ensemble `experiment` given a vector of parameter structs (vector of Params())
"""
function run_fixed_ensemble(experiment::String, vop::Vector)
    # First, obtain simulations
    perm = vop
    nperms = length(vop)

    # Now, create ensemble directory
    if isdir(pwd() * "/output/" * experiment)
        rm(pwd() * "/output/" * experiment, recursive=true)
    end
    mkdir(pwd() * "/output/" * experiment)
    mkdir(pwd() * "/output/" * experiment * "/runs/")
    mkdir(pwd() * "/output/" * experiment * "/results")

    # Save combinations and names in permutations.txt
    file_perm = open(pwd() * "/output/" * experiment * "/results/permutations.txt", "w")

    # Third, run each permutation
    printstyled("Running $(experiment): $(nperms) runs\n", color=:green)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in 1:nperms
        experimenti = "/runs/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        line2print = "run $(i)/$(nperms) s$(repeat("0", ndigits - length(digits(i))))$i"
        write(file_perm, "$(line2print) \n")

        run_pacco(experiment * experimenti, p=perm[i])
        println(line2print)
    end
    close(file_perm)
end

############
# Shortcuts
############
"""
    runplot_pacco(experiment)
runs and plot given `experiment` using parameters in `p`. Use `?run_pacco` for help.
"""
function runplot_pacco(experiment::String; p::Params=Params())
    run_pacco(experiment; p=p)
    plot_pacco(experiment)
    plot_pacco_states(experiment)
    plot_pacco_comp_states(experiment)
end

"""
    runplot_pacco_ensemble(experiment, params2per)
runs and plot given an ensemble named `experiment` using parameters in `p`. Use `?run_pacco` for help.
"""
function runplot_pacco_ensemble(experiment::String, params2per::Dict)
    run_pacco_ensemble(experiment, params2per)
    plot_pacco(experiment)
end

"""
    runplot_pacco_lhs(experiment, params2per, nsim)
runs and plot given an ensemble named `experiment` using parameters in `p` (Latin Hypercube Sampling mode). Use `?run_pacco` for help.
"""
function runplot_pacco_lhs(experiment::String, params2per::Dict, nsim::Int)
    run_pacco_lhs(experiment, params2per, nsim)
    plot_pacco(experiment)
end



function display_shortcuts()
    printstyled("Some available shortcuts: \n", color=:blue)
    println("* run_fixed_ensemble(experiment, vop)")
    println("* runplot_pacco(experiment)")
    println("* runplot_pacco_ensemble(experiment, params2per)")
    println("* runplot_pacco_lhs(experiment, params2per, nsim)")
    println("* fastplot(experiment, y; x, use_colormap, plot_function)")
    println("* plot_pacco_states(experiment)")
    println("* plot_pacco_comp_states(experiment)")
    println("")
    println("Note: Use ?function in REPL for more info about the arguments")
end