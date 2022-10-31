function run_amod(ctl, par, now)
    # =============================
    #     Module: run_amod
    #     Model: AMOD (adimensional ice-sheet-sediment model)
    #            by Jorge Alvarez-Solas (Fortran, 2017)
    # =============================
    #     Julia adapted by Sergio Pérez-Montero

    # Load src modules
    #include("./amod_orbital.jl")
    #include("./amod_dynamics.jl")
    #include("./amod_thermodynamics.jl")

    if ctl["dt"] <= 0.0 
        write(f, "Error with timestepping...")
        write(f, "dt   =  "*string(ctl["dt"])) 
        write(f,"time = "*string(now["time"]))
    end

end