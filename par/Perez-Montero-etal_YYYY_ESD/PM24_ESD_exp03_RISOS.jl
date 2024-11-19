Base.@kwdef mutable struct Params{T<:AbstractFloat}
  time_init::T = -1.0e6
  time_end::T = 0.0
  dt::T = 10.0
  dt_out::T = 1000.0
  time_spinup::T = 0.0
  degK::T = 273.15
  T0::T = 273.15
  C0::T = 280.0
  iceage0::T = 0.0
  albedo0::T = 0.2
  H0::T = 0.0
  Hsed0::T = 30.0
  B0::T = 0.0
  Tice0::T = 258.15
  fstr0::T = 0.2
  active_iso::Bool = true
  active_sed::Bool = false
  active_climate::Bool = false
  active_ice::Bool = true
  active_thermo::Bool = false
  active_aging::Bool = true
  active_snow_on_ice::Bool = true
  dt_case::String = "fixed" 
  insol_case::String = "laskar" 
  regtemp_case::String = "dynamic" 
  carbon_case::String = "dynamic" 
  deformationalflow_case::String = "profile" 
  plugflow_case::String = "profile" 
  deformational_case::String = "glen58" 
  basal_case::String = "weertmanq" 
  albedo_case::String = "prognostic" 
  snowfall_case::String = "linear" 
  ablation_case::String = "PDD" 
  diffusion_case::String = "2pts" 
  streaming_case::String = "fixed" 
  ref_streaming_case::String = "thermo" 
  ice_hrz_scale_case::String = "dynamic" 
  ice_discharge_case::String = "fixed" 
  temperate_layer_case::String = "fixed" 
  insol_day::T = 170.0
  insol_lat::T = 65.0
  insol_input::String = "data/insolation/solstice_insolation_65N170_10yr_5MyrBP-0.jld2" 
  insol_const::T = 400.0
  insol_min::T = 425.0
  insol_max::T = 565.0
  insol_threshold::T = 300.0
  insol_ref::T = 480.0
  Ppre::T = 0.6
  taupre::T = 23000.0
  Pobl::T = 0.4
  tauobl::T = 41000.0
  Pexc::T = 0.0
  tauexc::T = 100000.0
  At::T = 15.0
  cI::T = 0.1
  cISI::T = 1.0e-7
  cCAL::T = 1.7e-8
  cC::T = 0.65
  cZ::T = 0.007
  tauT::T = 900.0
  kT::T = -1.2e-6
  Tref::T = 273.15
  deltaT::T = 0.0
  Cref::T = 280.0
  tauC::T = 10.0
  kCT::T = 2.0
  kTC::T = 7.0
  kC::T = -1.0e-5
  deltaC::T = 0.0
  time_anth::T = 2000.0
  C_anth::T = 3000.0
  AT_anth::T = 55.0
  AC_anth::T = 20.0
  tau_anth::T = 200000.0
  albedo_land::T = 0.2
  albedo_oldice::T = 0.4
  albedo_newice::T = 0.9
  k_albedo::T = 5.0e-6
  albedo_quad::T = 1.0e-10
  tau_albedo::T = 1000.0
  pr_ref::T = 1.0
  Apr::T = 0.5
  Tsnow::T = 261.54999999999995
  Train::T = 280.54999999999995
  sref::T = 0.2
  ks::T = 0.0025
  lambda::T = 0.7
  Tthreshold::T = 268.15
  km::T = 0.0
  kI::T = 0.027
  Ath::T = 20.0
  L::T = 1.0e6
  hrz_scale_coeff::T = 0.9
  hrz_scale_exp::T = 2.0
  hrz_scale_lb::T = 223606.79774997896
  hrz_scale_ub::T = 4.0e6
  hrz_vel_scale::T = 300.0
  vrt_vel_scale::T = -0.01
  Aflow::T = 1.0e-16
  glen_n::T = 3.0
  Cs::T = 1.0e-10
  beta_min::T = 0.0
  sintheta::T = 0.001
  Hsed_max::T = 30.0
  Hsed_min::T = 0.0
  fv::T = 1.5e-7
  fa::T = 5.0e-6
  Beq::T = 0.0
  taubedrock::T = 5000.0
  hgeo::T = 0.05
  Tstr::T = 278.15
  cice::T = 2009.0
  kthr::T = 2.1
  Tmp::T = 273.15
  basal_scale::T = 10.0
  Pe::T = 10.0
  fstrmin::T = 0.4
  fstrmax::T = 0.4
  vkin::T = 1000.0
  taukin::T = 1000.0
  year_len::T = 365.2422
  sec_year::T = 3.1556926e7
  I0::T = 1365.2
  g::T = 9.81
  Surfoc::T = 3.618e8
  rhoice::T = 910.0
  muice0::T = 1.928e-5
  rhowater::T = 1000.0
  rhobed::T = 2700.0
  Γ::T = 0.0065
  Tfreezing::T = 273.15
  Lv::T = 2.5e6
  Rd::T = 287.0
  Rv::T = 461.0
  ice_exists_thr::T = 0.0
  ice_is_big_thr::T = 100.0
  ice_is_old_thr::T = 10.0
end 
