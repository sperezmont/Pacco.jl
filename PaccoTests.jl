"""
    run_tests(; testDyn, testClim, testDynClim)
Shortcut to run some general tests (dev), all attributes are optional
"""
function run_tests(; testDyn=false, testClim=false, testDynClim=false)
    isdir(pwd() * "/output/tests/") || mkdir(pwd() * "/output/tests/")

    # Test 1. Just ice dynamics and artificial insolation
    if testDyn
        expname = "tests/test1_JustDynamics"
        PaccoDynParams = Params(time_init=-2.0e6,
            Hsed0=1.0,
            active_sed=true, active_iso=true, active_climate=false, active_ice=true,
            I_case="artificial", M_case="PDD", At=25, f1=5e-8, f2=1e-6,
            lambda=0.05, ka=0.008, Aref=0.4)
        run_pacco(expname, p=PaccoDynParams)
        plot_pacco(expname, vars2plot=["I", "H", "B", "Hsed", "U"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 1 $(expname) completed\n", color=:green)
    end

    # Test 2. Just climate model, laskar forced
    if testClim
        expname = "tests/test2_JustClimatics"
        PaccoClimParams = Params(time_init=-8.0e5,
            Hsed0=1.0,
            active_sed=false, active_iso=true, active_climate=true, active_ice=false,
            I_case="laskar", M_case="ITM",
            alpha_slope=2.5e-6,
            ci=0.065, cc=0.65, cz=0.0065,
            lambda=0.01, ki=0.009, ka=0.0085, Aref=0.15)

        run_pacco(expname, p=PaccoClimParams)
        plot_pacco(expname, vars2plot=["I", "H", "T", "co2", "V"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 2 $(expname) completed\n", color=:green)

        expname = "tests/test2_JustClimatics_tests"
        PaccoClimParams = Dict("time_init" => [-8.0e5],
            "Hsed0" => [1.0],
            "active_sed" => [false], "active_iso" => [true], "active_climate" => [true], "active_ice" => [false],
            "I_case" => ["laskar"], "M_case" => ["ITM"],
            "alpha_slope" => [2.5e-6],
            "ci" => [0.065], "cc" => [0.65], "cz" => [0.0065],
            "lambda" => [0.01], "ki" => [0.009], "ka" => [0.0085], "Aref" => [0.15])

        run_pacco_ensemble(expname, PaccoClimParams)
        plot_pacco(expname, vars2plot=["I", "H", "T", "co2", "V"], plot_MPT=true, plot_PSD=true)
    end

    # Test 3. Climate-Cryosphere coupled
    if testDynClim
        # expname = "tests/test3_DynClim"
        # PaccoDynClimParams = Params(time_init = -2.0e6,
        #                          Hsed0 = 1.0,
        #                          active_sed = true, active_iso = true, active_climate = true, active_ice = true,
        #                          I_case = "laskar", M_case = "ITM",
        #                          ci = 0.1, cc = 0.65, cz = 0.008,
        #                          lambda = 0.064, ki = 0.025, ka = 0.02, Aref = 0.3)
        # run_pacco(expname, p = PaccoDynClimParams)
        # plot_pacco(expname, vars2plot=["I", "H", "Z", "B", "U", "Hsed", "alpha", "T"], plot_MPT=true, plot_PSD=true)
        # plot_wavelet(expname, MPT=true, fs=1 / 1000, sigma=π)
        # printstyled("test 3 $(expname) completed\n", color=:green)

        # expname = "tests/test3_DynClim"
        # PaccoDynClimParams = Params(time_init=-2.0e6,
        #     Hsed0=1.0,
        #     active_sed=true, active_iso=true, active_climate=true, active_ice=true,
        #     I_case="laskar", M_case="ITM",
        #     ci=0.1, cc=0.65, cz=0.0065,
        #     lambda=0.01, ki=0.05, ka=0.008, Aref=0.4)
        # run_pacco(expname, p=PaccoDynClimParams)
        # plot_pacco(expname, vars2plot=["Ianom", "A", "M", "U", "H", "Hsed", "T"], plot_MPT=true, plot_PSD=true)
        # plot_wavelet(expname, MPT=true, fs=1 / 1000, sigma=π)
        # printstyled("test 3 $(expname) completed\n", color=:green)

        expname = "tests/test3_DynClim_tests"
        PaccoDynClimParams = Dict("time_init" => [-2.0e6],
            "Hsed0" => [1.0],
            "active_sed" => [true], "active_iso" => [true], "active_climate" => [true], "active_ice" => [true],
            "I_case" => ["laskar"], "M_case" => ["ITM"],
            #"alpha_slope" => 2e-6:0.1e-6:5e-6,
            "ci" => [0.1], "cc" => [0.65], "cz" => [0.0065],
            "lambda" => [0.01], "ki" => [0.009, 0.025, 0.05], "ka" => [0.001, 0.005, 0.008], "Aref" => [0.4])
        run_pacco_ensemble(expname, PaccoDynClimParams)
        plot_pacco(expname, vars2plot=["Ianom", "A", "M", "U", "H", "Hsed", "T"], plot_MPT=true, plot_PSD=true)
    end


end