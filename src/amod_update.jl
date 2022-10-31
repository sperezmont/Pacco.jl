# =============================
#     Program: amod_update.jl
#     Aim: functions to update amod variables
# =============================
using DataStructures

function update_amod(old_vals::OrderedDict)
    return OrderedDict("time"      => old_vals["time"]::Real,            # updates model variables at t = time
                "H"         => old_vals["H"]::Real,
                "Hsed"      => old_vals["Hsed"]::Real,
                "T"         => old_vals["T"]::Real,
                "A"         => old_vals["A"]::Real,    
                "T_sl"      => old_vals["T_sl"]::Real, 
                "TMB"       => old_vals["TMB"]::Real,       
                "S"         => old_vals["S"]::Real,       
                "B"         => old_vals["B"]::Real,           
                "M"         => old_vals["M"]::Real,    
                "Acc"       => old_vals["Acc"]::Real,   
                "SMB"       => old_vals["SMB"]::Real,
                "U_d"       => old_vals["U_d"]::Real,       
                "U_b"       => old_vals["U_b"]::Real,  
                "U"         => old_vals["U"]::Real,      
                "T_surf"    => old_vals["T_surf"]::Real,
                "tau_b"     => old_vals["tau_b"]::Real,
                "tau_d"     => old_vals["tau_d"]::Real,     
                "Q_dif"     => old_vals["Q_dif"]::Real,  
                "Q_difup"   => old_vals["Q_difup"]::Real,
                "Q_difdown" => old_vals["Q_difdown"]::Real,      
                "Q_drag"    => old_vals["Q_drag"]::Real,   
                "alpha"     => old_vals["alpha"]::Real,  
                "Q_adv"     => old_vals["Q_adv"]::Real,               
                "fstream"   => old_vals["fstream"]::Real,
                "fstream_ref"=> old_vals["fstream_ref"]::Real,
                "Hdot"      => old_vals["Hdot"]::Real,   
                "Hseddot"   => old_vals["Hseddot"]::Real,  
                "Bdot"      => old_vals["Bdot"]::Real,  
                "Tdot"      => old_vals["Tdot"]::Real,   
                "fstreamdot"=> old_vals["fstreamdot"]::Real,
                "ins"       => old_vals["ins"]::Real,
                "co2"       => old_vals["co2"]::Real
                )
end

function update_amod_out(d::OrderedDict, vals::OrderedDict)
    for (key, value) in vals
        push!(d[key], vals[key])
    end
    return d
end

# OLD CODE

# function amod_init(runset_vars::Dict, incond_vars::Dict, par_vars::Dict)
#     # This function assigns initial conditions to the model variables
#     init_vals = Dict(
#     "time"    => runset_vars["time_init"],        
#     "H"       => incond_vars["H_init"],
#     "Hsed"    => incond_vars["Hsed_init"],
#     "T"       => incond_vars["T_init"],
#     "A"       => incond_vars["A_init"],
#     "T_sl"    => par_vars["T_ref"],
#     "TMB"     => 0.0,
#     "S"       => 0.0,
#     "TMB"     => 0.0,
#     "S"       => 0.0,                
#     "B"       => 0.0,
#     "M"       => 0.0,
#     "Acc"     => 0.0,
#     "SMB"     => 0.0,
#     "U_d"     => 0.0,
#     "U_b"     => 0.0,
#     "U"       => 0.0,
#     "T_surf"  => 0.0,
#     "tau_b"   => 0.0,
#     "tau_d"   => 0.0,
#     "Q_dif"   => 0.0,
#     "Q_difup" => 0.0,
#     "Q_difdown" => 0.0,
#     "Q_drag"  => 0.0,
#     "alpha"   => 0.0,
#     "Q_adv"   => 0.0,
#     "fstream" => 0.0,
#     "fstream_ref" => 0.0,
#     "Hdot"    => 0.0,
#     "Hseddot" => 0.0,
#     "Bdot"    => 0.0,
#     "Tdot"    => 0.0,
#     "fstreamdot" => 0.0,
#     "ins"     => par_vars["ins_prei"],
#     "co2"     => par_vars["co2_prei"]
#     )
#     now = update_amod(init_vals)
#     return now
# end

# function assign_parameters()    # assign the parameters of the namelist
#     ctl_assign = RunSettings(time_init, time_end, dt, dt_out)
#     inicond_assign = InitialConditions(H_init, Hsed_init, T_init, A_init)
#     par_assign = Parameters(ins_min, ins_max, ins_prei, co2_prei, T_0, T_ref,                               # radiative parameters
#                             orbital_mode, P_obl, tau_obl, P_pre, tau_pre, P_exc, tau_exc,                   # orbital parameters
#                             active_iso, B_eq, tau_bed, T_mantle, Q_geo,                                     # geophysical parameters
#                             v_kin, k_1, f_1, f_2, fstream_min, fstream_max, L, C_s, k, c, A_m, A_t, lambda) # ice parameters
#     return ctl_assign, inicond_assign, par_assign
# end

# function amod_init()    # initialize the model
    
#     now_init = NOW(ctl.time_init,        # time
#                 inicond.H_init,     # H
#                 inicond.Hsed_init,  # H_sed
#                 inicond.T_init,     # T
#                 inicond.A_init,     # A
#                 par.T_ref,       # T_sl
#                 0.0,                # TMB                
#                 0.0,                # S
#                 0.0,                # B
#                 0.0,                # M
#                 0.0,                # Acc
#                 0.0,                # SMB
#                 0.0,                # U_d
#                 0.0,                # U_b
#                 0.0,                # U
#                 0.0,                # T_surf
#                 0.0,                # tau_b
#                 0.0,                # tau_d
#                 0.0,                # Q_dif
#                 0.0,                # Q_difup
#                 0.0,                # Q_difdown
#                 0.0,                # Q_drag
#                 0.0,                # alpha
#                 0.0,                # Q_adv
#                 0.0,                # fstream
#                 0.0,                # fstream_ref
#                 0.0,                # Hdot
#                 0.0,                # Hseddot
#                 0.0,                # Bdot
#                 0.0,                # Tdot
#                 0.0,                # fstreamdot
#                 par.ins_prei,    # ins
#                 par.co2_prei)    # co2
#     return now_init
# end
