################################
#
# Leaf object 
# 
# AWalker Dec 2018
#
################################

source('leaf_system_functions.R')
source('leaf_functions.R')
source('leaf_solver_functions.R')



# LEAF OBJECT
###############################################################################

# use generic template
setwd('..')
source('generic_model_object.R')
leaf_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('leaf')



# assign object functions
###########################################################################
leaf_object$name <- 'leaf'


# function to configure unique elements of the object
# - adds functions to fns that are not in fnames
# - or functions that are derivations of other functions, in this case teh rs derived fuinctions like rs_r0 and rs_fe 
####################################
leaf_object$configure_unique <- function(., init=F, flist=NULL ) {
  if(init) {
    source('../../functions/general_functions.R')
    .$fns$puniroot            <- puniroot
    .$fns$assimilation        <- f_assimilation
    .$fns$assim_no_resistance <- f_solver_analytical_leaf_no_r
    .$fns$transition_cc       <- transition_cc
    .$fns$analytical_simple   <- f_solver_analytical_leaf_simple 
    .$fns$analytical_quad     <- f_solver_analytical_leaf_quad 
    .$fns$iterate_guesses     <- f_iterate_guesses 
    .$fns$residual_iterate    <- f_residual_iterate 
    .$fns$write_residual      <- f_write_residual 
    .$fns$solver_numerical    <- f_solver_brent 
  }

  if(any(names(flist)=='rs')) {
   .$fns$rs_fe <- get(paste0(.$fnames$rs,'_fe'), pos=1 )
   .$fns$rs_r0 <- get(paste0(.$fnames$rs,'_r0'), pos=1 )
  }
}



# assign object variables 
###########################################################################

# function names
####################################
leaf_object$fnames <- list(
  sys            = 'f_sys_enzymek', 
  solver         = 'f_solver_brent',
  residual_func  = 'f_residual_func_leaf_Ar',
  semiana        = 'f_semiana_quad',
  Acg            = 'f_Acg_farquhar1980',
  Ajg            = 'f_Ajg_generic',
  Apg            = 'f_Apg_vonc2000',            
  etrans         = 'f_j_harley1992',
  gas_diff       = 'f_ficks_ci',
  Alim           = 'f_lim_farquhar1980',
  vcmax          = 'f_vcmax_lin',
  jmax           = 'f_jmax_power',
  tcor_jmax      = 'f_scalar_none',
  tpu            = 'f_tpu_lin',
  rd             = 'f_rd_lin_vcmax',
  rl_rd_scalar   = 'f_scalar_none',
  gstar          = 'f_gstar_f1980',
  ri             = 'f_r_zero',
  rs             = 'f_rs_medlyn2011',
  rb             = 'f_rb_leafdim',
  cica_ratio     = 'f_cica_constant',             
  d13c           = 'f_d13c_classical',             
  tcor_asc = list(
    vcmax          = 'f_tcor_asc_Arrhenius',
    jmax           = 'f_tcor_asc_Arrhenius',
    tpu            = 'f_tcor_asc_Arrhenius',
    rd             = 'f_tcor_asc_Arrhenius',
    gstar          = 'f_tcor_asc_quadratic_bf1985',
    tau            = 'f_tcor_asc_Q10',
    Kc             = 'f_tcor_asc_Arrhenius',
    Ko             = 'f_tcor_asc_Arrhenius'
  ),
  tcor_des = list(
    vcmax          = 'f_tcor_des_modArrhenius',
    jmax           = 'f_tcor_des_modArrhenius',
    tpu            = 'f_tcor_des_modArrhenius',
    rd             = 'f_scalar_none'
  ),
  tcor_dep = list(
    tpu            = 'f_tcor_dep_independent',
    rd             = 'f_tcor_dep_independent',
    tau            = 'f_tcor_dep_independent'
  ),
  deltaS   = list(  
    rd             = 'f_deltaS',
    vcmax          = 'f_deltaS',
    jmax           = 'f_deltaS',
    tpu            = 'f_deltaS'
  ),
  q10_func = list(
    rd             = 'f_q10_constant',
    vcmax          = 'f_q10_constant',
    jmax           = 'f_q10_constant',
    tau            = 'f_q10_constant',
    Kc             = 'f_q10_constant',
    Ko             = 'f_q10_constant'
  )      
)


# environment
####################################
leaf_object$env <- list(
  ca_conc   = 400,                 # (umol mol-1)
  o2_conc   = 0.21,                # ( mol mol-1)
  par       = 1000,                # (umol photons m-2 s-1)
  water_l   = numeric(1),          # (mm) water level relative to hollow surfwce
  sphag_l   = 0,                   # (mm) Sphagnum surface relative to hollow surfwce
  temp      = 25,                  # (oC)
  vpd       = 1,                   # (kPa)
  rh        = numeric(1),          # (unitless - proportion)
  atm_press = 101325,              # ( Pa)
  wind      = 1                    # (m s-1)
)


# state
####################################
leaf_object$state <- list(
  #environmental state
  oi = numeric(1),                 # atmospheric & internal O2  (kPa)
  ca = numeric(1),                 # atmospheric CO2            ( Pa)
  cb = numeric(1),                 # boundary layer CO2         ( Pa)
  ci = numeric(1),                 # leaf internal CO2          ( Pa) 
  cc = numeric(1),                 # chloroplast CO2            ( Pa)
  leaf_temp = numeric(1),          # leaf temperature           (oC)
  
  #leaf state - calculated by canopy object so need initialisation
  leafN_area = 2,                  # leaf N per unit area       (g N m-2)
  fwdw_ratio = 5,                  # fresh weight dry weight ratio, used for Sphagnum conductance term 
  
  #calculated state
  J   = numeric(1),                # electron transport rate                            (umol electrons m-2 s-1) 
  Acg = numeric(1),                # Carboxylaton limited rate of net asssimilation     (umol m-2 s-1)
  Ajg = numeric(1),                # light limited rate of carboxylation                (umol m-2 s-1)
  Apg = numeric(1),                # TPU limited rate of carboxylation                  (umol m-2 s-1)
  A   = numeric(1),                # actual rate of carboxylation                       (umol m-2 s-1)
  rd  = numeric(1),                # actual rate of respiration                         (umol m-2 s-1)
  lim = numeric(1),                # flag indicationg limitation state of assimilation, wc = wc limited, wj = wj limited, wp = wp limited

  # diagnostic state
  A_ana_rbzero   = numeric(1),     # rate of carboxylation assuming zero boundary layer resistance to CO2 diffusion (umol m-2 s-1)
  A_ana_rbg0zero = numeric(1),     # rate of carboxylation assuming zero boundary layer resistance & zero minimum stomatal conductance to CO2 diffusion (umol m-2 s-1)
  aguess         = numeric(4),     # value of three guesses and solution from semi-analytical solver  (umol m-2 s-1)
  faguess        = numeric(4),     # value of solver function on three guesses and solution from semi-analytical solver  (umol m-2 s-1)
  aguess_flag    = numeric(1),     # flag for first guess from semi-analytical solver
  iter           = numeric(1),     # number of iterations to solve solver function in Brent uniroot
  estimprec      = numeric(1),     # estimated precision of solve from uniroot solver function 
  assim          = numeric(2),     # roots of fitted quadratic in semi-analytical solver
  fA_ana_final   = numeric(2),     # value of solver function for final gusee from semi-analytical solver  (umol m-2 s-1)
  A_noR          = numeric(1),     # rate of carboxylation assuming zero resistance to CO2 diffusion (umol m-2 s-1)
  transition     = numeric(1),     # cc at the transition point where wc = wj                        (Pa)
  d13c           = numeric(1)      # delta 13 C concentration of assimilated C 
)


# results from solver
####################################
leaf_object$solver_out = NULL


# state parameters (i.e. calculated parameters)
####################################
leaf_object$state_pars <- list(
  vcmax    = numeric(1),   # umol m-2 s-1
  vcmaxlt  = numeric(1),   # umol m-2 s-1
  jmax     = numeric(1),   # umol m-2 s-1
  jmaxlt   = numeric(1),   # umol m-2 s-1
  tpu      = numeric(1),   # umol m-2 s-1
  tpult    = numeric(1),   # umol m-2 s-1
  rd       = numeric(1),   # umol m-2 s-1; respiration at reftemp
  Kc       = numeric(1),   #  Pa
  Ko       = numeric(1),   # kPa
  Km       = numeric(1),   #  Pa
  gstar    = numeric(1),   #  Pa
  tau      = numeric(1),   #  -
  gamma    = numeric(1),   #  Pa
  rb       = numeric(1),   # m2s mol-1 H2O
  rs       = numeric(1),   # m2s mol-1 H2O
  ri       = numeric(1),   # m2s mol-1 CO2     
  alpha    = numeric(1),   # mol electrons mol-1 absorbed photosynthetically active photons
  cica_chi = numeric(1)    # Ci:Ca ratio 
)


# parameters
####################################
leaf_object$pars   <- list(
  diag          = F,          # calculate diagnostic output during runtime and add to output, such as cc transition point and non-stomatal limited assimilation rate 
  d13c          = F,          # calculate d13c and add to output
  deltaA_prop   = 0.15,       # proportion of first guess in A to use as delta in semi-analytical solver (unitless)  
  solver_min    = -0.0029834, # lower bracket for numerical solver (arbitrary decimal places to avoid numerical errors) 
  solver_max    = 51.8364435, # upper bracket for numerical solver  

  # photosynthetic parameters
  a             = 0.80,       # fraction of PAR absorbed by leaf                       (unitless)  --- this should equal 1 - leaf scattering coefficient, there is potential here for improper combination of models
  f             = 0.23,       # fraction of absorbed PAR not collected by photosystems (unitless)
  ko_kc_ratio   = 0.21,       # ratio of RuBisCO turnover numbers for oxgenation and carboxylation (unitless)
  theta_j       = 0.90,       # curvature of J quadratic in Farqhuar & Wong 1984       (unitless)
  theta_col_cj  = 0.95,       # curvature of 1st limitation quadratic, theta, in Collatz 1991  (unitless)
  theta_col_cjp = 0.98,       # curvature of 2nd limitation quadratic, beta, in Collatz 1991   (unitless)
  avn_25        = 10,         # intercept of linear vcmax25 to leaf N relationship     (umolm-2s-1)
  bvn_25        = 30,         # slope of linear vcmax25 to leaf N relationship         (umolm-2s-1g-1 N)
  ajv_25        = 29,         # intercept of linear jmax25 to vcmax25 relationship     (umolm-2s-1)
  bjv_25        = 1.63,       # slope of linear jmax25 to vcmax25 relationship         (unitless)
  a_jvt_25      = 2.59,       # intercept of linear jmax25:vcmax25 relationship to temperature   (e C-1)
  b_jvt_25      = -0.035,     # slope of linear jmax25:vcmax25 relationship to temperature       (e C-1 oC-1)
  e_ajv_25      = 1.01,       # intercept of log-log jmax25 to vcmax25 relationship    (log(umolm-2s-1))
  e_bjv_25      = 0.89,       # slope of log-log jmax25 to vcmax25 relationship        (unitless)
  atv_25        = 0,          # intercept of linear tpu25 to vcmax25 relationship      (umolm-2s-1)
  btv_25        = 1/6,        # slope of linear tpu25 to vcmax25 relationship          (unitless)
  flnr          = 0.09,       # fraction of leafN in RuBisCO -- PFT specific           (unitless)
  fnr           = 7.16,       # ratio of RuBisCO molecular mass to N in RuBisCO        (g RuBisCO g-1 N)
  Rsa           = 60,         # specific activity of RuBisCO                           ()
  Apg_alpha     = 0,          # alpha in tpu limitation eq, often set to zero check Ellesworth PC&E 2014 (unitless)

  # resistance parameters
  g0            = 0.01,       # Medlyn 2011 min gs                                     (molm-2s-1)
  g1_medlyn     = 6,          # Medlyn 2011 gs slope                                   (kPa^0.5)
  g1_leuning    = 10,         # Leuning 1995 gs slope                                  (unitless - likely higher than medlyn and ball g1)
  d0            = 1,          # Leuning 1995 D0                                        (kPa)
  g1_ball       = 6,          # Ball 1987 gs slope                                     (unitless - multiplier on RH as a proportion)
  g_a1_yin      = 0.85,       # Yin and Struik 2009 VPD response intercept             (unitless)
  g_b1_yin      = 0.14,       # Yin and Struik 2009 VPD response slope                 (kPa-1)
  rs            = 1/0.15,     # stomatal resistance                                    (m2s mol-1 h2o)
  cica_chi      = 0.7,        # constant Ci:Ca ratio                                   (unitless)
  rb            = 1/10,       # leaf boundary layer resistance                         (m2s mol-1 h2o)
  can_ttc       = 0.01,       # turbulent transfer coefficient between canopy surface and canopy air (m s-0.5)
  leaf_width    = 0.1,        # leaf dimension perpendicular to wind direction         (m)
  ri            = 1/0.15,     # mesophyll resistance                                   (m2s mol-1 - expressed in these units for consistency with other resistance terms, often expressed in the literature multiplied by Pa)
  co2_diff      = 1.7e-9,     # CO2 diffusivity in water                      - these three parameters are from Evans etal 2009 and the diffusivities are temp dependent  
  hco_co2_ratio = 0,          # ratio of HCO and CO2 concentration in water, assumed 0 for bog pH i.e. below 4.5   
  hco_co2_diff_ratio = 0.56,  # ratio of HCO and CO2 diffusivity in water  
  d13c_a        = 4.4,        # 13 C discrimination during CO2 diffusion through stomata 
  d13c_b        = 27,         # 13 C discrimination during carboxylation in the simple model 
  d13c_b_prime  = 30,         # 13 C discrimination during carboxylation 
  d13c_am       = 1.8,        # 13 C discrimination during diffusion to site of carboxylation 
  d13c_f        = 12,         # 13 C discrimination during photorespiration 

  # respiration parameters
  a_rdv_25      = 0,          # intercept of linear rd25 to vcmax25 relationship        (umolm-2s-1)
  b_rdv_25      = 0.015,      # slope of linear rd25 to vcmax25 relationship            (unitless)
  a_rdn_25      = 0.5,        # intercept of linear rd25 to leaf N area relationship    (umolm-2s-1)
  b_rdn_25      = 0.15,       # slope of linear rd25 to leaf N area relationship        (unitless)
  rl_rd_ratio   = 1,          # ratio of non-photorespiratory respiration in the light to respiration in the dark  (unitless)
  rl_rd_lloyd_a = 0.5,        # intercept of rl to rd scalar relationship to ln(PAR) Lloyd 1995 taken from mercado 2007 (unitless) 
  rl_rd_lloyd_b = 0.05,       # slope of rl to rd scalar relationship to ln(PAR) Lloyd 1995 taken from mercado 2007 (?)
  a_rdv_25_t    = 0.015,      # intercept of b_rdv_25 relationship to temperature       (umolm-2s-1)
  b_rdv_25_t    = -0.0005,    # slope of b_rdv_25 relationship to temperature           (unitless)
  
  # temperature response parameters
  reftemp = list(
    rd    = 25,               # reference temperature at which rd scalar = 1            (oC) 
    vcmax = 25,               # reference temperature at which Vcmax scalar = 1         (oC) 
    jmax  = 25,               # reference temperature at which Jmax scalar = 1          (oC)
    tpu   = 25,               # reference temperature at which TPU scalar = 1           (oC)
    Kc    = 25,               # reference temperature at which Kc scalar = 1            (oC)
    Ko    = 25,               # reference temperature at which Ko scalar = 1            (oC)
    gstar = 25,               # reference temperature at which gamma star scalar = 1    (oC)
    tau   = 25                # reference temperature at which tau scalar = 1           (oC)
  ),
  atref = list(
    rd      = 2,              # rd at ref temp (usually 25oC)    - used to set rd as a parameter                        (umolm-2s-1) 
    vcmax   = 50,             # vcmax at ref temp (usually 25oC) - used to set Vcmax as a parameter instead of an f(N)  (umolm-2s-1) 
    jmax    = 100,            # jmax at ref temp (usually 25oC)  - used to set Jmax as a parameter instead of an f(N)   (umolm-2s-1)
    tpu     = 5,              # tpu at ref temp (usually 25oC)   - used to set TPU as a parameter                       (umolm-2s-1)
    Kc      = 40.49,          # Kc for RuBisCO at ref temp (usually 25oC)               ( Pa)
    Ko      = 27.84,          # Kc for RuBisCO at ref temp (usually 25oC)               (kPa)
    gstar   = 4.325,          # Gamma star at ref temp (usually 25oC), 4.325 is Farquhar & Brooks value converted to Pa (Pa)
    tau     = 2600,           # CO2/O2 specificity ratio at ref temp (usually 25oC), Collatz 1991 (-)
    vomax   = numeric(1) 
  ),
  Ha = list(
    rd         = 69830,       # activation energy of respiration                        (J mol-1)
    vcmax      = 69830,       # activation energy of Vcmax                              (J mol-1)
    jmax       = 100280,      # activation energy of Jmax                               (J mol-1)
    tpu        = 69830,       # activation energy of TPU                                (J mol-1)
    Kc         = 79430,       # activation energy of Kc                                 (J mol-1)
    Ko         = 36380,       # activation energy of Ko                                 (J mol-1)
    gstar      = 37830,       # activation energy of gamma star                         (J mol-1)
    tau        = -41572,      # activation energy of tau                                (J mol-1)
    vomax      = 60110        # activation energy of Vomax                              (J mol-1)i
  ),
  Hd = list(
    rd         = 200000,      # deactivation energy of rd                               (J mol-1)
    vcmax      = 200000,      # deactivation energy of Vcmax                            (J mol-1)
    jmax       = 200000,      # deactivation energy of Jmax                             (J mol-1)
    tpu        = 200000       # deactivation energy of TPU                              (J mol-1)i
  ),
  Topt = list(
    rd       = 27.56,         #  temperature optimum of rd                               (oC)
    vcmax    = 27.56,         #  temperature optimum of Vcmax                            (oC)
    jmax     = 19.89,         #  temperature optimum of Jmax                             (oC)
    tpu      = 27.56          #  temperature optimum of TPU                              (oC)
  ),
  deltaS = list(
    rd     = numeric(1),      # 
    vcmax  = numeric(1),      # 
    jmax   = numeric(1),      #
    tpu    = numeric(1)       #
  ),
  a_deltaS_t = list(
    rd     = 490,             # linear temperature response of rd deltaS   
    vcmax  = 668,             # linear temperature response of vcmax deltaS (Kattge & Knorr)  
    jmax   = 660,             # linear temperature response of jmax  deltaS (Kattge & Knorr)
    tpu    = 485              # linear temperature response of tpu   deltaS (Kattge & Knorr)
  ),
  b_deltaS_t = list(
    rd     = 0,               # linear temperature response of rd deltaS
    vcmax  = -1.07,           # linear temperature response of vcmax deltaS (Kattge & Knorr)
    jmax   = -0.75,           # linear temperature response of jmax deltaS (Kattge & Knorr)
    tpu    = 0                # linear temperature response of tpu deltaS (Kattge & Knorr)
  ),
  q10 = list(
    rd        = 2,            # Q10 of Rd                                               (-)
    vcmax     = 2,            # Q10 of Vcmax                                            (-)
    jmax      = 2,            # Q10 of Jmax                                             (-)
    tpu       = 2,            # Q10 of TPU                                              (-)
    Kc        = 2,            # Q10 of Kc                                               (-)
    Ko        = 2,            # Q10 of Ko                                               (-)
    tau       = 0.57          # Q10 of tau                                              (-)
  ),
  a_q10_t = list(
    rd    = 3.22              # linear temperature response of rd Q10 (Tjoelker etal 2001)
  ),
  b_q10_t = list(
    rd    = -0.046            # linear temperature response of rd Q10 (Tjoelker etal 2001)
  ),
  tupp_cox = list(
    vcmax = 36,               # upper leaf T for Vcmax temp scaling from Cox 2001       (oC)
    rd    = 45                # upper leaf T for rd temp scaling from Cox 2001 (LM3)    (oC)
  ),
  tlow_cox = list(
    vcmax= 0,                 # lower leaf T for Vcmax temp scaling from Cox 2001       (oC)
    rd   = 5                  # lower leaf T for rd temp scaling from Cox 2001          (oC)
  ),
  exp_cox = list(
    vcmax = 0.3,              # exponent for Vcmax temp scaling from Cox 2001           (-)
    rd    = 0.4               # exponent for rd temp scaling from Cox 2001              (-)
  ),
  gstar_bf_a    = 0.012,      # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
  gstar_bf_b    = 1.68,       # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
  gstar_bf_c    = 42.7,       # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
  
  #physical constants
  R   = 8.31446               # molar gas constant                                      (m2 kg s-2 K-1 mol-1  ==  Pa m3 mol-1K-1)
)


# run control parameters
####################################
leaf_object$cpars <- list(
  verbose       = F,          # write diagnostic output during runtime 
  verbose_loop  = F,          # write diagnostic output on the solver during runtime 
  cverbose      = F,          # write configuration output during runtime 
  output        = 'slim'      # type of output from run function
)


# output functions
#######################################################################        
leaf_object$output <- function(.){

  lout <- 

    if(.$cpars$output=='run') {
      
      c(.$state_retrive(snames=c('A','cc','ci','rd','lim')),
                .$state_retrive(snames=c('ri','rs','rb'), state='state_pars') )
      
    } else if(.$cpars$output=='slim') {
      
      .$state_retrive(snames=c('A','ci','lim')) 
      
    } else if(.$cpars$output=='main'|.$cpars$output=='mcmc') {
      
      c(A=.$state$A)
      
    } else if(.$cpars$output=='state') {
      
      unlist(.$state)
      
    } else if(.$cpars$output=='full') {
      
      unlist(c(.$state, .$state_pars ))
      
    } else if(.$cpars$output=='all_lim') {
      
      c(.$state_retrive(snames=c('A','Acg','Ajg','Apg','cc','ci','ca','rd','lim')), 
                .$state_retrive(snames=c('ri','rs','rb'), state='state_pars' ) )
      
    } else if(.$cpars$output=='canopy') {
      
      c(.$state_retrive(snames=c('A','Acg','Ajg','Apg','cc','ci','ca','rd','lim')), 
                gi=1/.$state_pars$ri, gs=1/.$state_pars$rs, gb=1/.$state_pars$rb, 
                g =1/ (.$state_pars$rs + .$state_pars$rb) )
      
    } else if(.$cpars$output=='sphagnum') {
      
      c(.$state_retrive(snames=c('A','Acg','Ajg','Apg','cc','ci','ca','rd','lim','fwdw_ratio')), 
              .$state_retrive(snames=c('ri','rs','rb'), state='state_pars' ) )

    } else {

      stop(paste('No', .$name, 'model output option:', .$cpars$output, ', defined in output function.' ))
    }
  
  if(.$pars$diag&.$cpars$output!='full') c(lout, A_noR=.$state$A_noR, transition=.$state$transition ) 
  else                                   lout
  
}    


# test functions
#######################################################################        

leaf_object$.test <- function(., verbose=T, verbose_loop=T, leaf.par=1000, leaf.ca_conc=300, rs='f_rs_medlyn2011' ) {
  
  if(verbose) {
    str(.)
    print(.$env)
  }
  .$cpars$verbose       <- verbose
  .$cpars$verbose_loop  <- verbose_loop
  .$cpars$output        <-'full'
  
  .$fnames$rb            <- 'f_r_zero'
  .$fnames$ri            <- 'f_r_zero'
  .$fnames$rs            <- rs
  .$fnames$residual_func <- 'f_residual_func_leaf_Ar'
  
  .$env$par     <- leaf.par
  .$env$ca_conc <- leaf.ca_conc

  .$configure_test()
  .$run()
}


leaf_object$.test_residual_func <- function(., verbose=T, verbose_loop=T,
                                            centrala=5, range=3, inc=range/20,
                                            leaf.par=200, leaf.ca_conc=300, rs='f_rs_medlyn2011' ) {
  
  if(verbose) str(.)
  
  .$cpars$verbose       <- verbose
  .$cpars$verbose_loop  <- verbose_loop

  .$fnames$ri            <- 'f_r_zero'
  .$fnames$rs            <- rs
  .$fnames$rb            <- 'f_r_zero'
  .$fnames$residual_func <- 'f_residual_func_leaf_Ar'
  .$pars$g0              <- 0.01 
  
  # configure methods
  .$configure_test()

  # initialise the model without running the solution by setting PAR to zero  
  .$env$ca_conc        <- leaf.ca_conc
  .$env$par            <- 0 
  .$run()

  # calculate electron transport rate
  .$env$par <- leaf.par
  .$state$J <- .$fns$etrans()

  if(verbose) {
    print(.$fnames)
    print(.$state_pars)
    print(.$env)
  }

  # run the residual function within a loop
  out <- .$fns$residual_iterate(centrala, range, inc)

  # output
  print(xyplot(resid~a, as.data.frame(out), type='b', abline=0 ))
  out
}

    
leaf_object$.test_tscalar <- function(., leaf.temp=0:50, leaf.par=c(1000), leaf.ca_conc=400, rs='f_rs_medlyn2011',
                          tcor_asc='f_tcor_asc_Arrhenius', tcor_des='f_scalar_none', Ha=70000, 
                          verbose=F,verbose_loop=F) {
  
  .$cpars$verbose       <- verbose
  .$cpars$verbose_loop  <- verbose_loop
  .$cpars$output        <- 'full'
  
  if(verbose) str(.)
  
  .$fnames$tcor_asc[['vcmax']]  <- tcor_asc
  .$fnames$tcor_des[['vcmax']]  <- tcor_des
  .$pars$Ha$vcmax        <- Ha
  .$fnames$ri            <- 'f_r_zero'
  .$fnames$rs            <- rs
  .$fnames$residual_func <- 'f_residual_func_leaf_Ar'
  .$fnames$solver        <- 'f_solver_brent'
  
  # configure methods
  .$configure_test()

  # configure met data and run
  .$dataf     <- list()
  .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par','leaf.temp')))      
  .$dataf$out <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
  
  p1 <- xyplot(I(vcmaxlt/vcmax) ~ leaf.temp | as.factor(paste(tcor_asc,tcor_des)), .$dataf$out, abline=list(h=c(0,1), v=.$pars$reftemp[['vcmax']]),
               ylab=expression('scalar'),xlab=expression(T*' ['^o*C*']'))
  print(p1)
}


leaf_object$.test_aci <- function(., leaf.par=c(100,1000), leaf.ca_conc=seq(0.1,1500,50), rs='f_rs_medlyn2011', rb='f_r_zero', 
                      verbose=F, verbose_loop=F, diag=F, output='all_lim' ) {
  
  .$cpars$verbose       <- verbose
  .$cpars$verbose_loop  <- verbose_loop
  .$pars$diag           <- diag
  .$cpars$output        <- output
  
  if(verbose) str(.)
  
  .$fnames$residual_func <- 'f_residual_func_leaf_Ar'
  .$fnames$solver        <- 'f_solver_brent'
  .$fnames$ri            <- 'f_r_zero'
  .$fnames$rb            <- rb
  .$fnames$rs            <- rs
  
  # configure methods
  .$configure_test()

  # configure met data and run
  .$dataf     <- list()
  .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))      
  .$dataf$out <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
  
  print(cbind(.$dataf$met,.$dataf$out))
  p1 <- xyplot(A~cc|as.factor(.$dataf$met$leaf.par),.$dataf$out,groups=unlist(lim),abline=0,
               ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),xlab=expression(C[c]*' [Pa]'),
               panel=function(subscripts=subscripts,...) {
                 if(diag) {
                   panel.abline(v=.$dataf$out$transition[subscripts][1])
                   panel.points(y=.$dataf$out$A_noR[subscripts],x=.$dataf$out$cc[subscripts],col='black')                       
                 }
                 panel.xyplot(subscripts=subscripts,...)
               })
  print(p1)
}


leaf_object$.test_aci_light <- function(.,leaf.par=seq(10,2000,50),leaf.ca_conc=seq(1,1200,50),rs='f_rs_medlyn2011',
                                        verbose=F,verbose_loop=F,output=F,diag=F) {
  
  .$cpars$verbose       <- verbose
  .$cpars$verbose_loop  <- verbose_loop
  .$pars$diag           <- diag
  .$cpars$output        <- 'all_lim'
  
  if(verbose) str(.)
  
  .$fnames$ri            <- 'f_r_zero'
  .$fnames$rs            <- rs
  .$fnames$residual_func <- 'f_residual_func_leaf_Ar'
  .$fnames$solver        <- 'f_solver_brent'
  
  # configure methods
  .$configure_test()

  # configure met data and run
  .$dataf          <- list()
  .$dataf$met      <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))
  .$dataf$out      <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
  .$dataf$out_full <- cbind(.$dataf$met,.$dataf$out)
  
  p1 <- xyplot(A~cc,.$dataf$out_full,subset=leaf.par==1010,abline=0,groups=unlist(lim),
               ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),xlab=expression(C[c]*' [Pa]'),
               panel=function(subscripts=subscripts,...) {
                 if(diag) {
                   panel.abline(v=.$dataf$out_full$transition[subscripts][1])
                   panel.points(y=.$dataf$out_full$A_noR[subscripts],x=.$dataf$out_full$cc[subscripts],col='black')                       
                 }
                 panel.xyplot(subscripts=subscripts,...)
               })
  
  p2 <- xyplot(A~leaf.par,.$dataf$out_full,subset=leaf.ca_conc==401,abline=0,groups=unlist(lim),
               ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),xlab=expression('PAR ['*mu*mol*' '*m^-2*s-1*']'),
               panel=function(subscripts=subscripts,...) {
                 if(diag) {
                   #                        panel.abline(v=.$dataf$out_full$transition[subscripts][1])
                   panel.points(y=.$dataf$out_full$A_noR[subscripts],x=.$dataf$out_full$leaf.par[subscripts],col='black')                       
                 }
                 panel.xyplot(subscripts=subscripts,...)
               })
  
  print(p1,split=c(1,1,2,1),more=T)
  print(p2,split=c(2,1,2,1),more=F)
  if(output) .$dataf$out_full
}


leaf_object$.test_aci_analytical <- function(., rs='f_rs_medlyn2011', 
                                             leaf.par=c(100,1000), leaf.ca_conc=seq(50,1200,50), 
                                             leaf.rb=0, leaf.g0=0.01, 
                                             ana_only=F, verbose=F, verbose_loop=F, diag=F ) {
  
  if(verbose) str(.)
  .$cpars$verbose      <- verbose
  .$cpars$verbose_loop <- verbose_loop
  .$pars$diag          <- diag
  .$cpars$output       <- 'all_lim'
  
  .$fnames$rs <- rs
  .$fnames$ri <- 'f_r_zero'
  .$fnames$rb <- 'f_rb_constant'
  #.$fnames$rb <- 'f_rb_leafdim'
  #.$fnames$rb <- 'f_r_zero'
  .$pars$rb   <- leaf.rb
  .$pars$g0   <- leaf.g0
  
  # configure met data 
  .$dataf     <- list()
  .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))      
  
  if(!ana_only) {
    .$fnames$solver <- 'f_solver_brent'
    .$configure_test()
    .$dataf$out     <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)) )
    num_soln        <- cbind(.$dataf$met, .$dataf$out, sol=.$fnames$solver )
  }

  .$fnames$solver <- 'f_solver_analytical_leaf_simple'
  .$configure_test()
  .$dataf$out     <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)) )
  ana_soln        <- cbind(.$dataf$met, .$dataf$out, sol=.$fnames$solver)

  .$fnames$solver <- 'f_solver_analytical_leaf_quad'
  .$configure_test()
  .$dataf$out     <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)) )
  ana_soln        <- rbind(ana_soln,  cbind(.$dataf$met, .$dataf$out, sol=.$fnames$solver) )


  odf <- if(ana_only) ana_soln else rbind(num_soln,ana_soln)
  
  if(!ana_only) {
    p1 <- xyplot(A~cc|as.factor(odf$leaf.par),odf,groups=sol,abline=0,
                 ylab=expression('A ['*mu*mol*' '*m^-2*s^-1*']'),xlab=expression(C[c]*' [Pa]'),
                 panel=function(subscripts=subscripts,...) {
                   if(diag) {
                     panel.abline(v=odf$transition[subscripts][1])
                     panel.points(y=odf$A_noR[subscripts],x=odf$cc[subscripts],col='black')                       
                   }
                   panel.xyplot(subscripts=subscripts,...)
                 })
    
    p2 <- xyplot(A~leaf.ca_conc|as.factor(odf$leaf.par),odf,groups=sol,abline=0,
                 main=rs,auto.key=T,
                 ylab=expression('A ['*mu*mol*' '*m^-2*s^-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'),
                 panel=function(subscripts=subscripts,...) {
                   if(diag) {
                     panel.abline(v=odf$transition[subscripts][1])
                     panel.points(y=odf$A_noR[subscripts],x=odf$cc[subscripts],col='black')                       
                   }
                   panel.xyplot(subscripts=subscripts,...)
                 })
    
    ol <- list(odf,p2,p1)        

  } else {
    p1 <- NULL
    p2 <- xyplot(A~leaf.ca_conc,odf,groups=as.factor(odf$leaf.par),abline=0,
                 main=rs,
                 ylab=expression('A ['*mu*mol*' '*m^-2*s^-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'),
                 panel=function(subscripts=subscripts,...) {
                   if(diag) {
                     panel.abline(v=odf$transition[subscripts][1])
                     panel.points(y=odf$A_noR[subscripts],x=odf$cc[subscripts],col='black')                       
                   }
                   panel.xyplot(subscripts=subscripts,...)
                 })
    
    ol <- list(odf,p2,p1)        
  }

  print(p2)
  ol
}


leaf_object$.test_aci_lim <- function(.,rs='f_rs_medlyn2011',et='f_j_farquharwong1984',leaf.par=c(100,1000),leaf.ca_conc=seq(100,1500,50), 
                                      ana_only=F,verbose=F,verbose_loop=F,diag=F) {
  
  .$cpars$verbose       <- verbose
  .$cpars$verbose_loop  <- verbose_loop
  .$pars$diag           <- diag
  .$cpars$output        <- 'all_lim'
  
  if(verbose) str(.)
  
  .$fnames$rs           <- rs
  .$fnames$ri           <- 'f_r_zero'
  .$fnames$gas_diff     <- 'f_ficks_ci'
  .$fnames$etrans       <- et
  
  .$dataf     <- list()
  .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))      
  
  .$fnames$Alim <- 'f_lim_farquhar1980'
  .$configure_test()
  .$dataf$out   <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
  Fout          <- cbind(.$dataf$met,.$dataf$out,Alim='F1980')
  
  .$fnames$Alim <- 'f_lim_collatz1991'
  .$configure_test()
  .$dataf$out   <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
  Cout          <- cbind(.$dataf$met,.$dataf$out,Alim='C1991')
  
  odf <- rbind(Fout,Cout)
  
  p1 <- xyplot(A~cc|as.factor(et),odf,groups=unlist(Alim),abline=0,
               ylab=expression('A ['*mu*mol*' '*m^-2*s^-1*']'),xlab=expression(C[c]*' [Pa]'),
               type='l',auto.key=list(x=0.9,y=0.1,corner=c(1,0),points=F,lines=T),
               strip=strip.custom(bg='grey90'),
               panel=function(subscripts=subscripts,...) {
                 if(diag) {
                   panel.abline(v=odf$transition[subscripts][1])
                   panel.points(y=odf$A_noR[subscripts],x=odf$cc[subscripts],col='black')                       
                 }
                 panel.xyplot(subscripts=subscripts,...)
               })
  
  p2 <- xyplot(A~leaf.ca_conc|as.factor(odf$leaf.par),odf,groups=unlist(Alim),abline=0,
               auto.key=T,
               ylab=expression('A ['*mu*mol*' '*m^-2*s^-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'),
               panel=function(subscripts=subscripts,...) {
                 if(diag) {
                   panel.abline(v=odf$transition[subscripts][1])
                   panel.points(y=odf$A_noR[subscripts],x=odf$cc[subscripts],col='black')                       
                 }
                 panel.xyplot(subscripts=subscripts,...)
               })
  
  print(p1)
  # print(p2)
  odf
}



### END ###
