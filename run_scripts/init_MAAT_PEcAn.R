################################
#
# MAAT Leaf Model - initialisation script
# 
# AWalker (walkerap@ornl.gov) 
# December 2015
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script initialises the leaf photosynthesis version of MAAT
# - by setting the values of these lists:

# Static Variables:

# fnames.static
# pars.static     
# env.static


# Dynamic Variables:

# fnames.var
# pars.var     
# env.var

# This script must set the values of the above 6 lists, even if their value is NA

###################################################################



### Define the functions, parameters, and environment variables that are to remain static through a simulation 
###############################

# define lists
leaf.fnames.static <- list(
  gstar_tcor  = 'f_temp_scalar_quadratic_bf1985',
  Kc_tcor     = 'f_temp_scalar_Arrhenius',
  Ko_tcor     = 'f_temp_scalar_Arrhenius',
  vcmax       = 'f_constant_vcmax',
  jmax        = 'f_jmax_walker2014',
  tpu         = 'f_constant_tpu',
  vcmax_tcor  = 'f_temp_scalar_no_response',
  jmax_tcor   = 'f_temp_scalar_no_response',
  etrans      = 'f_j_harley1992',
  wc          = 'f_wc_farquhar1980',
  wj          = 'f_wj_generic',
  wp          = 'f_wp_vonc2000',            
  gas_diff    = 'f_ficks_ci_bound0',
  respiration = 'f_rd_collatz1991',
  fwdw_ratio  = 'f_none',                   
  ri          = 'f_r_zero',
  rs          = 'f_r_zero',
  rb          = 'f_r_zero',
  solver      = 'f_R_Brent_solver',
  solver_func = 'f_A_r_leaf',
  Alim        = 'f_lim_farquhar1980'
)

leaf.pars.static <- list(
  a             = 0.80,       # fraction of PAR absorbed                               (unitless)  --- this should equal 1 - leaf scattering coefficient, there is potential here for improper combination of models
  f             = 0.23,       # fraction of absorbed PAR not collected by photosystems (unitless)
  theta         = 0.90,       # curvature of J quadratic in Farqhuar & Wong 1984       (unitless)
  theta_collatz = 0.98,       # curvature of 1st limitation quadratic in Collatz 1991  (unitless)
  beta_collatz  = 0.95,       # curvature of 2nd limitation quadratic in Collatz 1991  (unitless)
  avn_25        = 10,         # intercept of linear vcmax25 to leaf N relationship     (umolm-2s-1)
  bvn_25        = 30,         # slope of linear vcmax25 to leaf N relationship         (umolm-2s-1g-1 N)
  ajv_25        = 29,         # intercept of linear jmax25 to vcmax25 relationship     (umolm-2s-1)
  bjv_25        = 1.63,       # slope of linear jmax25 to vcmax25 relationship         (unitless)
  e_ajv_25      = 1.01,       # intercept of log-log jmax25 to vcmax25 relationship    (log(umolm-2s-1))
  e_bjv_25      = 0.89,       # slope of log-log jmax25 to vcmax25 relationship        (unitless)
  flnr          = 0.09,       # fraction of leafN in RuBisCO -- PFT specific           (unitless)
  fnr           = 7.16,       # ratio of RuBisCO molecular mass to N in RuBisCO        (g RuBisCO g-1 N)
  Rsa           = 60,         # specific activity of RuBisCO                           ()
  wp_alpha      = 0.05,       # alpha in tpu limitation eq, often set to zero check Ellesworth PC&E 2014 (unitless)
  # resistance parameters
  g0            = 0.01,       # Medlyn 2011 min gs                                     (molm-2s-1)
  g1_medlyn     = 5,          # Medlyn 2011 gs slope                                   ()
  g1_leuning    = 5,          # Leuning 1995 gs slope                                  ()
  d0            = 2,          # Leuning 1995 D0                                        ()
  g1_ball       = 5,          # Ball 1987 gs slope                                     ()
  gi            = 0.035,      # mesophyll conductance                                                (molm-2s-1Pa-1)
  ri            = 1/0.035,    # mesophyll resistance                                                 (m2sPa mol-1)
  co2_diff      = 1.7e-9,     # CO2 diffusivity in water                      - these three parameters are from Evans etal 2009 and the diffusivities are temp dependent  
  hco_co2_ratio = 0,          # ratio of HCO and CO2 concentration in water, assumed 0 for bog pH i.e. below 4.5   
  hco_co2_diff_ratio = 0.56,  # ratio of HCO and CO2 diffusivity in water  
  fwdw_wl_slope = -0.022,     # delta sphagnum fwdw ratio per mm of decrease in water level    (mm-1) , currently from Adkinson & Humpfries 2010, Rydin 1985 has similar intercept but slope seems closer to -0.6 
  fwdw_wl_sat   = 16,         # sphagnum fwdw ratio at 0 water level, currently from Adkinson & Humpfries 2010     
  fwdw_wl_exp_a = -0.037,     # decrease in sphagnum fwdw ratio as an exponential f of water level (cm), currently from Strack & Price 2009
  fwdw_wl_exp_b = 3.254,      # decrease in sphagnum fwdw ratio as an exponential f of water level (cm) 
  # respiration parameters
  rd            = 1.5,        # rd as a constant,                                      (umolm-2s-1)
  rd_prop_vcmax = 0.015,      # rd as a proportion of Vcmax, Williams & Flannagan 1998 ~ 0.1         (unitless)
  rd_prop_N     = 0.15,       # rd as a proportion of leaf N                           (umols-1g-1)
  # temperature response parameters
  reftemp.rd    = 25,         # reference temperature at which rd scalar = 1            (oC) 
  reftemp.vcmax = 25,         # reference temperature at which Vcmax scalar = 1         (oC) 
  reftemp.jmax  = 25,         # reference temperature at which Jmax scalar = 1          (oC)
  reftemp.Kc    = 25,         # reference temperature at which Kc scalar = 1            (oC)
  reftemp.Ko    = 25,         # reference temperature at which Ko scalar = 1            (oC)
  reftemp.gstar = 25,         # reference temperature at which gamma star scalar = 1    (oC)
  atref.rd      = 0,          # rd at ref temp (usually 25oC)    - used to set rd as a parameter                        (umolm-2s-1) 
  atref.vcmax   = @VCMAX@,    # vcmax at ref temp (usually 25oC) - used to set Vcmax as a parameter instead of an f(N)  (umolm-2s-1) 
  atref.jmax    = 0,          # jmax at ref temp (usually 25oC)  - used to set Jmax as a parameter instead of an f(N)   (umolm-2s-1)
  atref.tpu     = 10,         # tpu at ref temp (usually 25oC)   - used to set TPU as a parameter                       (umolm-2s-1)
  atref.Kc      = 40.49,      # Kc for RuBisCO at ref temp (usually 25oC)               ( Pa)
  atref.Ko      = 27.84,      # Kc for RuBisCO at ref temp (usually 25oC)               (kPa)
  atref.gstar   = 4.325,      # Gamma star at ref temp (usually 25oC), 4.325 is Farquhar & Brooks value converted to Pa (Pa)
  Ha.vcmax      = 69830,      # activation energy of Vcmax                              (J mol-1)
  Ha.jmax       = 100280,     # activation energy of Jmax                               (J mol-1)
  Ha.Kc         = 79430,
  Ha.Ko         = 36380,
  Ha.gstar      = 37830,
  Ha.vomax      = 60110,
  Hd.vcmax      = 200000,     # deactivation energy of Vcmax                            (J mol-1)
  Hd.jmax       = 200000,     # deactivation energy of Jmax                             (J mol-1)
  Topt.vcmax    = 27.56,      # temperature optimum of Vcmax                            (oC)
  Topt.jmax     = 19.89       # temperature optimum of Jmax                             (oC)
)

leaf.env.static  <- list(
  ca_conc   = 400,                 # (umol mol-1)
  o2_conc   = 0.21,                # ( mol mol-1)
  par       = 1000,                # (umol photons m-2 s-1)
  temp      = 25,                  # (oC)
  vpd       = 2,                   # (kPa)
  atm_press = 101325               # ( Pa)
)


### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths

# define lists
# leaf.fnames.var <- list(
#   etrans = c('f_j_farquhar1980','f_j_farquharwong1984','f_j_harley1992','f_j_collatz1991'),
#   ri     = c('f_r_zero','f_ri_constant')
# )
# 
# leaf.pars.var <- list(
#   avn_25 = rnorm(10,10,1)
# )
# 
# leaf.env.var <- list(
#   temp = c(5,20)
# )

leaf.fnames.var <- NA
leaf.pars.var   <- NA
leaf.env.var    <- NA



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_ls <- list(
  lfs = leaf.fnames.static,
  lps = leaf.pars.static,
  les = leaf.env.static,
  lfv = leaf.fnames.var,
  lpv = leaf.pars.var,
  lev = leaf.env.var
)


