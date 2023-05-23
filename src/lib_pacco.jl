# =============================
#     Program: lib_pacco.jl
#     Aim: main functions of PACCO
# =============================
"""
    calc_diagnostic_variables!(u, p, t)
computes diagnostic variables in `u` for time `t` using parameters in `p`
"""
function calc_diagnostic_variables!(u::Vector, p::Params, t::Real)
    # Compute Orbital Forcing
    calc_I!(u, p, t)

    # Translate Orbital Forcing to temperature
    if p.active_climate
        calc_R!(u, p)
    else
        calc_Tsl!(u, p, t)
    end
    calc_Tref!(u, p, t)

    # Compute Ice-Sheet Geometry
    calc_Z!(u, p)
    calc_E!(u, p)
    calc_V!(u, p)

    if p.active_climate
        # Compute Climate variables
        calc_alpha_ref!(u, p)
    end
    calc_Tsurf!(u, p)

    # Compute Ice-Sheet Mass Balance
    calc_A!(u, p)
    calc_M!(u, p)

    if p.active_ice
        # Compute Ice-Sheet dynamics
        calc_taud!(u, p)
        calc_taub!(u, p)
        calc_Ud!(u, p)
        calc_Ub!(u, p)
        calc_fstream_ref!(u, p)
        u[26] = u[23] + u[9] * u[24]

        # Compute Ice-Sheet thermodynamics
        calc_Qdif!(u, p)
        calc_Qdrag!(u, p)
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
        dudt[1] = calc_Tdot(u, p)

        # -- co2
        dudt[2] = calc_co2dot(u, p, t)

        # -- ice age
        dudt[3] = calc_iceagedot()

        # -- albedo
        dudt[4] = calc_alphadot(u, p)
    end

    # -- ice thickness
    dudt[5] = calc_Hdot(u, p)

    # -- sediments layer thickness
    (p.active_sed) && (dudt[6] = calc_Hseddot(u, p))

    # -- bed elevation
    (p.active_iso) && (dudt[7] = calc_Bdot(u, p))

    if p.active_ice
        # -- ice temperature
        (p.active_ice) && (dudt[8] = calc_Ticedot(u))

        # -- streaming fraction
        (p.active_ice) && (dudt[9] = calc_fstreamdot(u, p))
    end

    # -- diagnostic derivatives
    view(dudt, lprog+1:lsu) .= 0.0

    # Modify states to ensure physical meaning
    view(u, 1:lprog) .= max.(view(u, 1:lprog), [0.0, 1.0, 0.0, p.alpha_land, 0.0, 0.0, -Inf, 0.0, 0.0])
    view(u, 1:lprog) .= min.(view(u, 1:lprog), [Inf, Inf, Inf, Inf, Inf, 1.0, Inf, p.degK, Inf])

    if u[5] == 0.0  # no ice
        u[3] = 0.0
        u[4] = p.alpha_land
    elseif u[3] < 10.0 # first ice
        u[4] = p.alpha_newice
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
    solution = solve(prob, BS3(), saveat=p.dt_out)
    return solution
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

        run_pacco("test1", p = Params(time_init = -2e6, Aref = 0.1, lambda = 0.1))
Runs experiment `test1` using the default values of `p` located in `par/default_params.jl` and change those provided by `p = Params(time_init = -2e6, Aref = 0.1, lambda = 0.1)`

        run_pacco("test1", p = JLD2.load_object("path/to/julia/object/params.jld2"))
Runs experiment `test1` and uses as parameters the ones stored in `"path/to/julia/object/params.jld2"`
"""
function run_pacco(experiment::String; p::Params=Params())
    ## Now, load arguments
    output_path = pwd() * "/output/" * experiment * "/"
    isdir(output_path) || mkdir(output_path)
    write_run_info(output_path, experiment, p)
    JLD2.save_object(output_path * "/params.jld2", p)

    ## Run pacco()
    u0, out_attr = load_defs(p)
    out = pacco(u0, p, (p.time_init, p.time_end))
    JLD2.save_object(output_path * "/pacco.jld2", out)

    ## Create output nc file
    # if outfile exists remove it
    isfile(output_path * "/pacco.nc") && rm(output_path * "/pacco.nc")
    genout_nc(output_path, "/pacco.nc", out, p, out_attr)

    ## Done!
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
        println(line2print)
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

        new_pi = Params(; (Symbol(k) => v for (k, v) in permutations_dict[i])...)
        run_pacco(experiment * experimenti, p=new_pi)
        println(line2print)
    end
    # Done!
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
    runplot_pacco_ensemble(experiment)
runs and plot given an ensemble named `experiment` using parameters in `p`. Use `?run_pacco` for help.
"""
function runplot_pacco_ensemble(experiment::String, params2per::Dict)
    run_pacco_ensemble(experiment, params2per)
    plot_pacco(experiment)
end

"""
    runplot_pacco_lhs(experiment)
runs and plot given an ensemble named `experiment` using parameters in `p` (Latin Hypercube Sampling mode). Use `?run_pacco` for help.
"""
function runplot_pacco_lhs(experiment::String, params2per::Dict, nsim::Int)
    run_pacco_lhs(experiment, params2per, nsim)
    plot_pacco(experiment)
end