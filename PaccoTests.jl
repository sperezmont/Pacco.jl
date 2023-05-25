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
        plot_wavelet(expname, plot_MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 1 $(expname) completed\n", color=:green)
    end

    # Test 2. Just climate model, laskar forced
    if testClim
        expname = "tests/test2_JustClimate"
        PaccoClimParams = Params(time_init=-1.0e6,
            Hsed0=1.0,
            active_sed=false, active_iso=true, active_climate=true, active_ice=false,
            I_case="laskar", M_case="ITM",
            alpha_slope=5e-6, ktco2=10.0,
            ci=0.08, cc=0.65, cz=0.0065,
            lambda=0.1, ki=0.0095, ka=0.008, Aref=0.15)
        run_pacco(expname, p=PaccoClimParams)
        plot_pacco(expname, vars2plot=["Ianom", "H", "T", "co2", "V"], times=(-8e5, 0), plot_MPT=true, plot_PSD=true)
        plot_wavelet(expname, plot_MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 2 $(expname) completed\n", color=:green)
    end

    # Test 3. Coupled Climate-Cryosphere
    if testDynClim
        expname = "tests/test3_DynClim"
        PaccoDynClimParams = Params(time_init = -2e6,
                                 Hsed0 = 1.0,
                                 active_sed = true, active_iso = true, active_climate = true, active_ice = true,
                                 I_case = "laskar", M_case = "ITM",
                                 alpha_slope=5e-6, ktco2=10.0,
                                 ci = 0.08, cc = 0.65
                                 , cz = 0.005,
                                 lambda = 0.1, ki = 0.025, ka = 0.02, Aref = 0.3)
        run_pacco(expname, p = PaccoDynClimParams)
        plot_pacco(expname, vars2plot=["I", "H", "T", "co2", "V", "Hsed"], times=(-8e5,0), plot_MPT=true, plot_PSD=true)
        plot_wavelet(expname, plot_MPT=true, fs=1 / 1000, sigma=π)
        plot_pacco_comp_states("tests/test3_DynClim")
        printstyled("test 3 $(expname) completed\n", color=:green)
    end
end