# =============================
#     Program: lib_pacco.jl
#     Aim: main functions of PACCO
# =============================

"""
    update_forward(now, par, ctl, vars2update) 

takes vars2update and calculates their time evolution vars2update is a vector with the names of the variables to update vars2update is computed in run_pacco() function 

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters
* `ctl` Run control parameters
* `vars2update` Vector of variables to update in time as var = var + vardot * dt
## Return
updated `now` dictionary
"""
function update_forward(now::OrderedDict, par::OrderedDict, ctl::OrderedDict, vars2update::Vector)
    for hm in par["hemisphere"], v in vars2update
        # -- calculate time evolution
        variab, vardot = v * "_" * hm, v * "dot_" * hm
        now[variab] += now[vardot] * ctl["dt"]    # now = now + dnow/dt * dt

        # -- modify if desired
        if variab in ["H_n", "H_s"]
            now[variab] = max(now[variab], 0.0)
        elseif variab in ["Hsed_n", "Hsed_s"]
            now[variab] = min(max(now[variab], 0.0), 1.0) # sediments go from 0 to 1
        elseif variab in ["B_n", "B_s"]
            (~par["active_iso"]) && (now[variab] = par["B_eq_"*hm]) # reupdate to equilibrium value
        elseif variab in ["T_ice_n", "T_ice_s"]
            now[variab] = min(now[variab], degK)
        elseif variab in ["co2_n", "co2_s"]
            now[variab] = max(now[variab], 1.0)
        end
    end
    return now
end

"""
     pacco(now, par, ctl, vars2update)

PACCO model 

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `par` Dictionary with run parameters
* `ctl` Run control parameters
* `vars2update` Vector of variables to update in time as var = var + vardot * dt
## Return
updated `now` dictionary
"""
function pacco(now::OrderedDict, par::OrderedDict, ctl::OrderedDict, vars2update::Vector)
    # Update variables
    now = update_forward(now, par, ctl, vars2update)
    now = update_Z(now, par)

    now = calc_E(now, par)
    now = calc_V(now, par)

    # -- calculate reference value for albedo
    (par["active_climate"]) && (now = calc_albedo_ref(now, par))

    # -- thermomechanical coupling (For now, the rate factor A is constant -- jas)

    ## Ice Dynamics
    if par["active_ice"]
        now = calc_tau_d(now, par)           # -- driving stress
        now = calc_tau_b(now, par)           # -- basal shear stress

        now = calc_U_d(now, par)             # -- mean deformational velocity
        now = calc_U_b(now, par)             # -- mean basal velocity

        now = calc_fstream(now, par, ctl)   # -- stream fraction

        now = calc_U(now, par)              # -- total velocity
    end

    # -- surface temperature
    now = calc_T_surf(now, par)

    # -- bedrock temperature (currently prescribed -- jas)

    # -- mass balance
    now = calc_SMB(now, par)
    now = calc_TMB(now, par)

    ## Ice thermodynamics
    if par["active_ice"]
        now = calc_Qdif(now, par)       # -- diffusion
        now = calc_Qdrag(now, par)      # -- drag
    end

    # Calculate time evolution
    # -- ice thickness
    now = calc_Hdot(now, par)

    # -- sediment thickness
    if par["active_ice"]
        if par["active_sed"]
            now = calc_Hseddot(now, par)
        end
        # -- ice temperature
        now = calc_T_icedot(now, par)
    end

    # -- bedrock elevation
    now = calc_Bdot(now, par)

    # -- regional temperature
    (par["active_climate"]) && (now = calc_Tdot(now, par))

    # -- albedo
    now = calc_albedodot(now, par)

    # -- co2
    now = calc_co2dot(now, par)

    return now
end

"""
    pacco_loop(now, out, par, ctl, vars2update, file)
calculates the time loop of the model

## Arguments
* `now` Dictionary with values of the model variables at current timestep
* `out` Dictionary with stored values of the model variables (output)
* `par` Dictionary with run parameters
* `ctl` Run control parameters
* `vars2update` Vector of variables to update in time as var = var + vardot * dt
* `file` out.out file  
## Return
updated `out` dictionary
"""
function pacco_loop(now::OrderedDict, out::OrderedDict, par::OrderedDict, ctl::OrderedDict, vars2update::Vector, file)
    time_length = ceil((ctl["time_end"] - ctl["time_init"]) / ctl["dt"])
    for n in 1:time_length
        # update simulation time
        now["time"] = ctl["time_init"] + n * ctl["dt"]

        # update contour variables
        if par["active_climate"]
            now = calc_T_rf(now, par)     # -- rf due to coupled climate    
        else
            now = calc_T_sl(now, par)    # -- T_sl from insolation
        end

        # run PACCO
        now = pacco(now, par, ctl, vars2update)

        # only update output variable at desired frequency
        if mod(now["time"], ctl["dt_out"]) == 0
            out = update_pacco_out(now, out)
            if par["active_outout"]
                write(file, "time = " * string(now["time"]) * "\n")
            end
        end

        # Check if run presents instabilities
        break_iter = check_run(now)
        (break_iter) && (break)
    end
    return out
end


function write_run_info(outpath, expname, parname, changed_pars)
    isfile(outpath * "/run-info.txt") && rm(outpath * "/run-info.txt")
    f = open(outpath * "/run-info.txt", "w")
    write(f, "$(expname)\n")
    write(f, "Source parameter file: $(parname)\n")
    write(f, "Parameters used: \n")
    write(f, "$(changed_pars)\n")
    close(f)
end


"""
    run_pacco(; experiment, par_file, par2change)
        
main function of PACCO (all arguments are optional)

## Arguments 
* `experiment`  Name of the experiment to perform
* `par_file`    Name of the parameter file to use
* `pars2change` Dictionary with the (parameter => new value) to change in the original `par_file` parameter file
## Return
nothing
"""
function run_pacco(; experiment::String="test_default", par_file::String="pacco_default.jl", par2change=[])
    ## Now, load arguments
    output_path = load_out(pacco_path, experiment)
    load_parf(pacco_path, output_path, par_file)
    write_run_info(output_path, experiment, par_file, par2change)

    # -- check if some parameters need to be changed
    if par2change != []
        change_namelist(pacco_path * "/output/" * experiment, "namelist.jl", par2change)
    end

    # Assign parameters
    # -- assign parameters
    ctl, now, par, out, out_prec, out_attr = load_defs(output_path * "namelist.jl")

    ## Open out.out
    # if out.out exists remove it
    isfile(output_path * "/out.out") && rm(output_path * "/out.out")
    f = open(output_path * "/out.out", "w")
    if par["active_outout"]
        write(f, "**** Starting PACCO " * experiment * " ... ****\n")
    end

    ## Initialize
    out = update_pacco_out(now, out) # update output

    if par["active_outout"]
        write(f, "time = " * string(now["time"]) * "\n")
    end

    ## Let's run!
    # -- define updatable (in time) variables
    vars2update = ["H", "B"]
    (par["active_ice"]) && (vars2update = vcat(vars2update, ["Hsed", "T_ice"]))
    (par["active_climate"]) && (vars2update = vcat(vars2update, ["T", "co2", "albedo", "ice_time"]))

    # -- run PACCO
    out = pacco_loop(now, out, par, ctl, vars2update, f)

    ## Create output nc file
    # if outfile exists remove it
    isfile(output_path * "/pacco.nc") && rm(output_path * "/pacco.nc")
    genout_nc(output_path, "pacco.nc", out, out_prec, out_attr)

    if par["active_outout"]
        write(f, "**** PACCO " * experiment * " done! ****" * "\n")
    end
    close(f)
    (par["active_outout"] == false) && rm(output_path * "/out.out")

    ## Done!
end

"""
    run_pacco_ensemble(par2per; experiment, par_file)
        
Takes an (ordered) dictionary of parameters to be exchanged in order to run an ensemble and runs PACCO for each combination
Deprecated, use `run_pacco_lhs()` instead

## Arguments 
* `par2per`     Dictionary with parameters to be exchanged
## Optional Arguments
* `experiment`  Name of the experiment to perform
* `par_file`    Name of the parameter file to use
"""
function run_pacco_ensemble(par2per::OrderedDict; experiment::String="test_default_ens", par_file::String="pacco_default.jl")
    # First, obtain simulations
    perm = calc_permutations(par2per)

    # Second, get number of simulations
    nperms = length(perm)

    # Now, create ensemble directory
    if isdir(pacco_path * "/output/" * experiment)
        rm(pacco_path * "/output/" * experiment, recursive=true)
    end
    mkdir(pacco_path * "/output/" * experiment)
    mkdir(pacco_path * "/output/" * experiment * "/runs/")
    mkdir(pacco_path * "/output/" * experiment * "/results")

    # Save combinations and names in permutations.txt
    file_perm = open(pacco_path * "/output/" * experiment * "/results/permutations.txt", "w")

    # Third, run each permutation
    printstyled("Running $(experiment): $(nperms) runs\n", color=:green)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in 1:nperms
        experimenti = "/runs/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        line2print = "s" * repeat("0", ndigits - length(digits(i))) * "$i $(perm[i])"
        write(file_perm, "$(line2print) \n")

        if i == 1
            time1run = @elapsed run_pacco(experiment=experiment * experimenti, par_file=par_file, par2change=perm[i])
            time_needed, t_units = nperms * time1run / 60 / 60, "h"
            if time_needed < 1
                time_needed, t_units = nperms * time1run / 60, "m"
                if time_needed < 1
                    time_needed, t_units = nperms * time1run, "s"
                end
            end

            printstyled("It is expected to take $(time_needed) $(t_units) ($(time1run)s/it)\n", color=:blue)
        else
            run_pacco(experiment=experiment * experimenti, par_file=par_file, par2change=perm[i])
        end
        println(line2print)
    end
    close(file_perm)
    # Done!
end

"""
    run_pacco_lhs(par2per, nsim; experiment, par_file)
        
Takes a dictionary of parameters and create LHS, key => (min, max) to run an ensemble and runs PACCO for each combination

## Arguments 
* `par2per`     Dictionary with parameters to be exchanged, key => (min, max)
* `nsim`        Number of simulations to perform

## Optional Arguments
* `experiment`  Name of the experiment to perform
* `par_file`    Name of the parameter file to use
## Return
nothing
"""
function run_pacco_lhs(par2per::Dict, nsim::Int; experiment::String="test_default_ens", par_file::String="pacco_default.jl")
    parnames = collect(keys(par2per))
    # First, calculate LHS
    permutations, permutations_dict = gen_lhs(par2per, nsim; pars_type=2)

    # Now, create ensemble directory
    if isdir(pacco_path * "/output/" * experiment)
        rm(pacco_path * "/output/" * experiment, recursive=true)
    end
    mkdir(pacco_path * "/output/" * experiment)
    mkdir(pacco_path * "/output/" * experiment * "/runs/")

    # Plot (if possible) the LHS
    mkdir(pacco_path * "/output/" * experiment * "/results")
    if length(par2per) == 2
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2])
        scatter!(ax, permutations)
        save(pacco_path * "/output/" * experiment * "/results/lhs.png", fig)
    elseif length(par2per) == 3
        fig = Figure()
        ax = Axis3(fig[1, 1], xlabel=parnames[1], ylabel=parnames[2], zlabel=parnames[3])
        scatter!(ax, permutations)
        save(pacco_path * "/output/" * experiment * "/results/lhs.png", fig)
    end

    # Run each permutation
    printstyled("Running $(experiment): $(nperms) runs\n", color=:green)
    nperms = length(permutations_dict)
    ndigits = Int(round(log10(nperms)) + 1)
    for i in 1:nperms
        experimenti = "/runs/s" * repeat("0", ndigits - length(digits(i))) * "$i/"
        line2print = "s" * repeat("0", ndigits - length(digits(i))) * "$i $(perm[i])"

        if i == 1
            time1run = @elapsed run_pacco(experiment=experiment * experimenti, par_file=par_file, par2change=permutations_dict[i])
            time_needed, t_units = nperms * time1run / 60 / 60, "h"
            if time_needed < 1
                time_needed, t_units = time_needed * 60, "m"
                if time_needed < 1
                    time_needed, t_units = time_needed * 60 * 60, "s"
                end
            end
            printstyled("It is expected to take $(time_needed) $(t_units) ($(time1run)s/it)\n", color=:blue)
        else
            run_pacco(experiment=experiment * experimenti, par_file=par_file, par2change=permutations_dict[i])
        end
        println(line2print)
    end
    # Done!
end


# Some shortcuts
