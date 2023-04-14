"""
    run_tests(; test1a, test1b, test1c, test2, test3)
        
Shortcut to run some general tests (dev), all attributes are optional

## Arguments (default value is `false`)
* `test1a`  Just ice dynamics and artificial insolation
* `test1b`  Just ice dynamics, artificial insolation and linear Clausius-clapeyron
* `test1c`  Just ice dynamics, laskar insolation and linear Clausius-clapeyron
* `test2`   Just climate
* `test3`   Default mode, ice dynamics + coupled climate
* `test4`   Runs an ensemble of ice+climate with different weights in ice divergence
* `all_tests` Runs all tests if set to `true`
## Return
nothing
"""
function run_tests(; test1=false, test2=false, test3=false, test4=false, testx=false, testEGU23=false, all_tests=false)
    isdir(pwd() * "/output/tests/") || mkdir(pwd() * "/output/tests/")
    if all_tests == true
        test1 = true
        test2 = true
        test3 = true
        test4 = true
    end

    # Test 1. Just ice dynamics and artificial insolation
    if test1
        expname = "tests/test1_ice-artif"
        parfile = "test_par_files/test_ice-artif.jl"
        run_pacco(experiment=expname, par_file=parfile)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "Hsed_n", "B_n"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 1 $(expname) completed, needs further improvement in calibration\n", color=:green)
    end

    # Test 2. Just ice dynamics and real insolation
    if test2
        expname = "tests/test2_ice-real"
        parfile = "test_par_files/test_ice-real.jl"
        run_pacco(experiment=expname, par_file=parfile)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "Hsed_n", "B_n"], plot_MPT=true, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 2 $(expname) completed, needs further improvement in calibration\n", color=:green)
    end

    # Test 3. Just climate
    if test3
        expname = "tests/test3_clim"
        parfile = "test_par_files/test_clim.jl"
        run_pacco(experiment=expname, par_file=parfile)
        plot_pacco(experiment=expname, vars2plot=["ins_anom_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 3 $(expname) completed \n Note: needs further improvement in calibration\n", color=:green)

        expname = "tests/test3_clim_2"  # this version includes some really good combinations with H > 2000 m
        parfile = "test_par_files/test_clim.jl"
        par2change = OrderedDict("Acc_ref_n" => [0.25, 0.3], "csi" => [0.1], "ki" => [0.009, 0.02], "ka" => [0.008, 0.01]) # 0.05 > ki > 0.01
        run_pacco_ensemble(par2change, experiment=expname, par_file=parfile)
        plot_pacco(experiment=expname, vars2plot=["ins_anom_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        printstyled("test 3 $(expname) completed \n Note: needs further improvement in calibration\n", color=:green)

    end

    if test4 # default mode, ice + coupled climate 
        parfile = "test_par_files/test_ice-clim.jl"
        pars2change = OrderedDict(
            "active_sed" => false,
            "height_temp" => "useZ",
            "f_1" => 1e-7,      
            "Acc_ref_n" => 0.3, # 0.3
            "csi" => 0.1,       # 0.1
            "cs" => 0.65,     # 0.65 
            "csz" => 0.0065,    # 0.0065
            "ki" => 0.025,      # 0.025
            "ka" => 0.02,       # 0.02
            "lambda" => 0.05)   # 0.05

        expname = "tests/test4_ice-clim_Hsed0"  # no sediments
        pars2change["time_init"], pars2change["time_end"] = -1e6, 0.0
        pars2change["active_sed"], pars2change["Hsed_init_n"] = false, 0.0
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        expname = "tests/test4_ice-clim_Hsed0_2"  # no sediments
        pars2change["time_init"], pars2change["time_end"] = -2e6, 0.0
        pars2change["active_sed"], pars2change["Hsed_init_n"] = false, 0.0
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        expname = "tests/test4_ice-clim_Hsed1" # with sediments
        pars2change["time_init"], pars2change["time_end"] = -2e6, -1e6
        pars2change["active_sed"], pars2change["Hsed_init_n"] = false, 1.0
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        expname = "tests/test4_ice-clim_Hsed1_2" # with sediments
        pars2change["time_init"], pars2change["time_end"] = -2e6, 0
        pars2change["active_sed"], pars2change["Hsed_init_n"] = false, 1.0
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        expname = "tests/test4_ice-clim" # active sediments
        pars2change["time_init"], pars2change["time_end"] = -2e6, 0.0
        pars2change["active_sed"], pars2change["Hsed_init_n"] = true, 1.0
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        expname = "tests/test4_ice-clim_Hsed0_3"  # no sediments
        pars2change["time_init"], pars2change["time_end"] = -2e6, 0.0
        pars2change["active_sed"], pars2change["Hsed_init_n"] = false, 0.0
        pars2change["lambda"] = 0.05
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        expname = "tests/test4_ice-clim_Hsed1_3"  # no sediments
        pars2change["time_init"], pars2change["time_end"] = -2e6, 0.0
        pars2change["active_sed"], pars2change["Hsed_init_n"] = false, 1.0
        pars2change["lambda"] = 0.05
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true)
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

        plot_pacco(experiments=["tests/test4_ice-clim", "tests/test4_ice-clim_Hsed1", "tests/test4_ice-clim_Hsed0"], vars2plot=["ins_n", "H_n", "T_n"], plot_PSD=true, plot_MPT=true)
        plot_pacco(experiments=["tests/test4_ice-clim_Hsed0", "tests/test4_ice-clim_Hsed0_2", "tests/test4_ice-clim_Hsed0_3"], vars2plot=["ins_n", "H_n", "T_n"], plot_PSD=true, plot_MPT=true)
        plot_pacco(experiments=["tests/test4_ice-clim_Hsed1", "tests/test4_ice-clim_Hsed1_2", "tests/test4_ice-clim_Hsed1_3"], vars2plot=["ins_n", "H_n", "T_n"], plot_PSD=true, plot_MPT=true)

    end

    if testEGU23
        parfile = "test_par_files/test_ice-clim.jl"
        pars2change = OrderedDict(
            "time_init" => -1.5e6,
            "time_end" => 0.0,
            "active_sed" => false,
            "height_temp" => "useZ",
            "Hsed_init_n" => 0.0,
            "f_1" => 1e-7,      
            "Acc_ref_n" => 0.3, # 0.3
            "csi" => 0.1,       # 0.1
            "cs" => 0.65,     # 0.65 
            "csz" => 0.00685,    # 0.0065
            "ki" => 0.025,      # 0.025
            "ka" => 0.02,       # 0.02
            "lambda" => 0.064)   # 0.05

        expname = "tests/test4_ice-clim_EGU23"  # no sediments
        run_pacco(experiment=expname, par_file=parfile, par2change=pars2change)
        plot_pacco(experiment=expname, vars2plot=["ins_n", "H_n", "T_n", "co2_n", "V_n"], plot_MPT=false, plot_PSD=true, plot_MIS=true, times=(-8e5, 0))
        plot_wavelet(experiment=expname, MPT=true, fs=1 / 1000, sigma=π)
        printstyled("test 4 $(expname) completed \n Note: Work in progress\n", color=:red)

    end

    if testx
        expname = "tests/testx" 
        parfile = "test_par_files/test_ice-clim.jl"
        pars2perm = OrderedDict(
            "active_sed" => [false],
            "height_temp" => ["useZ"],     
            "Acc_ref_n" => [0.3], # 0.3
            "csi" => [0.1],       # 0.1
            "cs" => [0.65],     # 0.65 
            "csz" => [0.00685],    # 0.0065
            "ki" => [0.025],      # 0.025
            "ka" => [0.02],       # 0.02
            "lambda" => 0.059:0.001:0.064)   # 0.05
        run_pacco_ensemble(pars2perm, experiment=expname, par_file=parfile)
        plot_pacco(experiment=expname, vars2plot=["ins_anom_n", "H_n", "T_n", "co2_n", "V_n"], plot_PSD=true, plot_MIS=true)
        printstyled("test x $(expname) completed \n Note: Work in progress\n", color=:red)
    end
end