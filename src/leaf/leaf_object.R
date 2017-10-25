################################
#
# Leaf object 
# 
# AWalker March 2014
#
################################

library(proto)
library(stringr)
#library(plyr)

source('leaf_functions.R')



# LEAF OBJECT
###############################################################################

leaf_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'leaf'
    
    # no expected child objects
    
    # build function - no expected child objects so returns null
    build <- function(.) {
      # as.proto(.$as.list())
      NULL
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.) {
      
#       if(.$cpars$cverbose) {
#         print('',quote=F)
#         print('Leaf configuration:',quote=F)
#         print(unlist(.$fnames),quote=F)
#         print(unlist(.$pars),quote=F)
#         print(unlist(.$env),quote=F)
#       }
      
      # calculate environmental state
      # CO2 partial pressure (Pa)
      .$state$ca <- .$env$ca_conc * .$env$atm_press * 1e-6
      # O2 partial pressure  (kPa)
      .$state$oi <- .$env$o2_conc * .$env$atm_press * 1e-3
      # leaf temperature
      .$state$leaf_temp <- .$env$temp               
      # RH from VPD
      .$env$rh   <- f_rh_from_vpd(.)

      # calculate state parameters
      # photosynthetic parameters
      .$state_pars$vcmax   <- get(.$fnames$vcmax)(.)
      .$state_pars$jmax    <- get(.$fnames$jmax)(.)
      .$state_pars$tpu     <- get(.$fnames$tpu)(.)
      .$state_pars$rd      <- get(.$fnames$respiration)(.)
      .$state_pars$alpha   <- 0.5 * (1-.$pars$f)
      
      # kinetic pars & temperature dependence
      # - if the solver involves the energy balance all of this crap needs to go in the solver! 
      .$state_pars$Kc      <- .$pars$atref.Kc * get(.$fnames$Kc_tcor)(.,parlist=list(Tr=.$pars$reftemp.Kc,Ha=.$pars$Ha.Kc,q10=.$pars$q10.Kc,q10_func='f_q10_constant'))
      .$state_pars$Ko      <- .$pars$atref.Ko * get(.$fnames$Ko_tcor)(.,parlist=list(Tr=.$pars$reftemp.Ko,Ha=.$pars$Ha.Ko,q10=.$pars$q10.Ko,q10_func='f_q10_constant')) 
      .$state_pars$Km      <- .$state_pars$Kc * (1+(.$state$oi/.$state_pars$Ko)) 
      .$state_pars$gstar   <- get(.$fnames$gstar)(.) 
      .$state_pars$vcmaxlt <- .$state_pars$vcmax * 
        get(.$fnames$vcmax_tcor_asc)(.,parlist=list(Tr=.$pars$reftemp.vcmax,Ha=.$pars$Ha.vcmax,q10=.$pars$q10.vcmax,q10_func=.$fnames$q10_func.vcmax)) *
        get(.$fnames$vcmax_tcor_des)(.,parlist=list(Tr=.$pars$reftemp.vcmax,Ha=.$pars$Ha.vcmax,Hd=.$pars$Hd.vcmax,Topt=.$pars$Topt.vcmax,
                                                    tupp=.$pars$tupp_cox.vcmax,tlow=.$pars$tlow_cox.vcmax,exp=.$pars$exp_cox.vcmax,
                                                    a_deltaS_t=.$pars$a_deltaS_t.vcmax,b_deltaS_t=.$pars$b_deltaS_t.vcmax,deltaS=.$pars$deltaS.vcmax))
      .$state_pars$jmaxlt  <- .$state_pars$jmax  * 
        get(.$fnames$jmax_tcor_asc)(.,parlist=list(Tr=.$pars$reftemp.jmax,Ha=.$pars$Ha.jmax,q10=.$pars$q10.jmax,q10_func=.$fnames$q10_func.jmax)) *
        get(.$fnames$jmax_tcor_des)(.,parlist=list(Tr=.$pars$reftemp.jmax,Ha=.$pars$Ha.jmax,Hd=.$pars$Hd.jmax,Topt=.$pars$Topt.jmax,q10=.$pars$q10.jmax,
                                                   a_deltaS_t=.$pars$a_deltaS_t.jmax,b_deltaS_t=.$pars$b_deltaS_t.jmax,deltaS=.$pars$deltaS.jmax))
      .$state_pars$tpult   <- .$state_pars$tpu   * get(.$fnames$tpu_tcor_dependence)(.) 

      # conductance/resistance terms
      # - if either of these functions become a function of co2 or assimilation they can be easily moved into the solver
      .$state_pars$rb      <- get(.$fnames$rb)(.)
      .$state_pars$ri      <- get(.$fnames$ri)(.)  
      
      # print state parameters to screen
      if(.$cpars$verbose) {
        print(.$state_pars)
      }  
      
      # calculate physiological state
      # electron transport rate
      .$state$J <- get(.$fnames$etrans)(.)
      # respiration
      .$state$respiration  <- .$state_pars$rd * get(.$fnames$rd_tcor_dependence)(.)
      
      # if PAR > 0
      if(.$env$par > 0) {
        # diagnostic calculations
        if(.$pars$diag) {
          .$state$A_noR      <- f_A_r_leaf_noR(.)
          .$state$transition <- transition_cc(.)
        }
        # run photosynthesis
        # account for decreased respiration in the light
        .$state$respiration  <- get(.$fnames$rl_rd_scalar)(.) * .$state$respiration
        .$state_pars$gamma   <- (-.$state_pars$vcmaxlt * .$state_pars$gstar - .$state$respiration * .$state_pars$Km) / (.$state$respiration - .$state_pars$vcmaxlt)
        # determine rate limiting step - this is done based on carboxylation, not net assimilation (Gu etal 2010).
        .$state$A       <- get(.$fnames$solver)(.)      
        # assign the limitation state a numerical code - assumes the minimum is the dominant limiting rate
        .$state$lim     <- c(2,3,7)[which(c(.$state$Acg,.$state$Ajg,.$state$Apg)==min(c(.$state$Acg,.$state$Ajg,.$state$Apg),na.rm=T))]       
        # after the fact calculations
        .$state$cb      <- f_ficks_ci(.,A=.$state$A,r=1.4*.$state_pars$rb,c=.$state$ca)
        if(!grepl('analytical',.$fnames$solver)) .$state_pars$rs <- get(.$fnames$rs)(.) 
        .$state$ci      <- f_ficks_ci(.,A=.$state$A,r=1.6*.$state_pars$rs,c=.$state$cb)
        .$state$cc      <- f_ficks_ci(.,A=.$state$A,r=.$state_pars$ri,c=.$state$ci)        
      }
      # if PAR < 0
      else {
        # assume infinite conductances when concentration gradient is small
        # - this ignores the build up of CO2 within the leaf due to respiration and high rs 
        .$state$cc <- .$state$ci <- .$state$cb <- .$state$ca
        .$state$A_noR      <- NA
        .$state$transition <- NA
        .$state$A          <- -.$state$respiration      
        .$state$lim        <- 0       
        .$state_pars$rs    <- get(.$fnames$rs)(.)
      }
      
      # print to screen
      if(.$cpars$verbose) {
        print(.$state)
      }
      
      # output
      .$output()
    } 
    
    
    
    ###########################################################################
    # Output functions

    # -- returns a vector of outputs
    state_retrive <- function(.,snames) {
      lsubs <- match(snames,names(.$state))
      unlist(.$state[lsubs])
    }

    output <- function(.){
      if(.$cpars$output=='slim') {
        
        lout <- .$state_retrive(snames=c('A','ci','lim')) 
        
      } else if(.$cpars$output=='run') {
        
        lout <- .$state_retrive(snames=c('A','cc','ci','ri','rs','rb','respiration','lim')) 
        
      } else if(.$cpars$output=='all_lim') {
        
        lout <- .$state_retrive(snames=c('A','Acg','Ajg','Apg','cc','ci','ca','ri','rs','rb','respiration','lim')) 
        
      } else if(.$cpars$output=='full') {
        
        lout <- unlist(c(.$state, .$state_pars ))
        
      } else if(.$cpars$output=='sphagnum') {
        
        lout <- .$state_retrive(snames=c('A','cc','ci','ri','rs','rb','respiration','lim','fwdw_ratio')) 
        
      }
      
      if(.$pars$diag&.$cpars$output!='full') c( lout, A_noR=.$state$A_noR, transition=.$state$transition ) 
      else                                    lout
      
    }    

    
    ###########################################################################
    # Variables etc
    
    # function names
    fnames <- list(
      gstar               = 'f_gstar_constref',
      gstar_tcor          = 'f_temp_scalar_quadratic_bf1985',
      tau_tcor            = 'f_temp_scalar_Q10',
      Kc_tcor             = 'f_temp_scalar_Arrhenius',
      Ko_tcor             = 'f_temp_scalar_Arrhenius',
      vcmax               = 'f_vcmax_lin',
      jmax                = 'f_jmax_power',
      tpu                 = 'f_tpu_lin',
      vcmax_tcor_asc      = 'f_temp_scalar_Arrhenius',
      jmax_tcor_asc       = 'f_temp_scalar_Arrhenius',
      vcmax_tcor_des      = 'f_temp_scalar_modArrhenius_des',
      jmax_tcor_des       = 'f_temp_scalar_modArrhenius_des',
      tpu_tcor_dependence = 'f_tpu_tcor_dependent',
      tpu_tcor_asc        = 'f_temp_scalar_Arrhenius',
      tpu_tcor_des        = 'f_temp_scalar_modArrhenius_des',
      deltaS              = 'f_deltaS',
      etrans              = 'f_j_harley1992',
      Acg                 = 'f_Acg_farquhar1980',
      Ajg                 = 'f_Ajg_generic',
      Apg                 = 'f_Apg_vonc2000',            
      gas_diff            = 'f_ficks_ci',
      respiration         = 'f_rd_lin_vcmax',
      rl_rd_scalar        = 'f_scalar_none',
      rd_tcor_dependence  = 'f_rd_tcor_dependent',
      rd_tcor_asc         = 'f_temp_scalar_Q10',
      rd_tcor_des         = 'f_temp_scalar_cox2001_des',
      q10_func.rd         = 'f_q10_constant',
      q10_func.vcmax      = 'f_q10_constant',
      q10_func.jmax       = 'f_q10_constant',
      fwdw_ratio          = 'f_none',                   
      cica_ratio          = 'f_cica_constant',             
      ri                  = 'f_r_zero',
      rs                  = 'f_r_zero',
      rb                  = 'f_r_zero',
      solver              = 'f_R_Brent_solver',
      solver_func         = 'f_A_r_leaf',
      Alim                = 'f_lim_farquhar1980'
    )
    
    # leaf environment
    env <- list(
      ca_conc   = numeric(1),          # (umol mol-1)
      o2_conc   = 0.21,                # ( mol mol-1)
      par       = numeric(1),          # (umol photons m-2 s-1)
      water_l   = numeric(1),          # (mm) water level relative to hollow surfwce
      sphag_l   = 0,                   # (mm) Sphagnum surface relative to hollow surfwce
      temp      = 25,                  # (oC)
      vpd       = 2,                   # (kPa)
      rh        = numeric(1),          # (unitless - proportion)
      atm_press = 101325               # ( Pa)
    )

    # leaf state
    state <- list(
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
      A            = numeric(1),       # actual rate of carboxylation                       (umol m-2 s-1)
      respiration  = numeric(1),       # actual rate of respiration                         (umol m-2 s-1)
      lim          = numeric(1),       # flag indicationg limitation state of assimilation, wc = wc limited, wj = wj limited, wp = wp limited

      # diagnostic state
      A_noR        = numeric(1),       # rate of carboxylation assuming zero resistance to CO2 diffusion (umol m-2 s-1)
      transition   = numeric(1)        # cc at the transition point where wc = wj                        (Pa)
    )
    
    # results from solver
    solver_out = NULL
    
    #leaf state parameters (i.e. calculated parameters)
    state_pars <- list(
      vcmax    = numeric(1),   # umol m-2 s-1
      vcmaxlt  = numeric(1),   # umol m-2 s-1
      jmax     = numeric(1),   # umol m-2 s-1
      jmaxlt   = numeric(1),   # umol m-2 s-1
      tpu      = numeric(1),   # umol m-2 s-1
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
    
    #leaf parameters
    pars   <- list(
      diag          = F,          # calculate diagnostic output during runtime and add to output, such as cc transition point and non-stomatal limited assimilation rate 
      # photosynthetic parameters
      # deprecated    alpha    = 0.24,         # harley 1992 alpha - Williams & Flannagan 1998 use 0.21 but calculate 0.25 
      a             = 0.80,       # fraction of PAR absorbed                               (unitless)  --- this should equal 1 - leaf scattering coefficient, there is potential here for improper combination of models
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
      btv_25        = 1/8.2,      # slope of linear tpu25 to vcmax25 relationship          (unitless)
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
      rs            = 1/0.15,     # stomatal resistance                                    (m2s mol-1 h2o)
      cica_chi      = 0.7,        # constant Ci:Ca ratio                                   (unitless)
      rb            = 1/10,       # leaf boundary layer resistance                         (m2s mol-1 h2o)
      ri            = 1/0.15,     # mesophyll resistance                                   (m2s mol-1 - expressed in these units for consistency with other resistance terms, often expressed in the literature multiplied by Pa)
      co2_diff      = 1.7e-9,     # CO2 diffusivity in water                      - these three parameters are from Evans etal 2009 and the diffusivities are temp dependent  
      hco_co2_ratio = 0,          # ratio of HCO and CO2 concentration in water, assumed 0 for bog pH i.e. below 4.5   
      hco_co2_diff_ratio = 0.56,  # ratio of HCO and CO2 diffusivity in water  
      fwdw_wl_slope = -0.022,     # delta sphagnum fwdw ratio per mm of decrease in water level      (mm-1), currently from Adkinson & Humpfries 2010, Rydin 1985 has similar intercept but slope seems closer to -0.6 
      fwdw_wl_sat   = 16,         # sphagnum fwdw ratio at 0 water level, currently from Adkinson & Humpfries 2010     
      fwdw_wl_exp_a = -0.037,     # decrease in sphagnum fwdw ratio as an exponential f of water level (cm), currently from Strack & Price 2009
      fwdw_wl_exp_b = 3.254,      # decrease in sphagnum fwdw ratio as an exponential f of water level (cm) 
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
      reftemp.rd    = 25,         # reference temperature at which rd scalar = 1            (oC) 
      reftemp.vcmax = 25,         # reference temperature at which Vcmax scalar = 1         (oC) 
      reftemp.jmax  = 25,         # reference temperature at which Jmax scalar = 1          (oC)
      reftemp.tpu   = 25,         # reference temperature at which TPU scalar = 1           (oC)
      reftemp.Kc    = 25,         # reference temperature at which Kc scalar = 1            (oC)
      reftemp.Ko    = 25,         # reference temperature at which Ko scalar = 1            (oC)
      reftemp.gstar = 25,         # reference temperature at which gamma star scalar = 1    (oC)
      reftemp.tau   = 25,         # reference temperature at which tau scalar = 1           (oC)
      atref.rd      = 2,          # rd at ref temp (usually 25oC)    - used to set rd as a parameter                        (umolm-2s-1) 
      atref.vcmax   = 50,         # vcmax at ref temp (usually 25oC) - used to set Vcmax as a parameter instead of an f(N)  (umolm-2s-1) 
      atref.jmax    = 100,        # jmax at ref temp (usually 25oC)  - used to set Jmax as a parameter instead of an f(N)   (umolm-2s-1)
      atref.tpu     = 5,          # tpu at ref temp (usually 25oC)   - used to set TPU as a parameter                       (umolm-2s-1)
      atref.Kc      = 40.49,      # Kc for RuBisCO at ref temp (usually 25oC)               ( Pa)
      atref.Ko      = 27.84,      # Kc for RuBisCO at ref temp (usually 25oC)               (kPa)
      atref.gstar   = 4.325,      # Gamma star at ref temp (usually 25oC), 4.325 is Farquhar & Brooks value converted to Pa (Pa)
      atref.tau     = 2600,       # CO2/O2 specificity ratio at ref temp (usually 25oC), Collatz 1991 (-)
      atref.vomax   = numeric(1),
      Ha.rd         = 69830,      # activation energy of respiration                        (J mol-1)
      Ha.vcmax      = 69830,      # activation energy of Vcmax                              (J mol-1)
      Ha.jmax       = 100280,     # activation energy of Jmax                               (J mol-1)
      Ha.tpu        = 69830,      # activation energy of TPU                                (J mol-1)
      Ha.Kc         = 79430,      # activation energy of Kc                                 (J mol-1)
      Ha.Ko         = 36380,      # activation energy of Ko                                 (J mol-1)
      Ha.gstar      = 37830,      # activation energy of gamma star                         (J mol-1)
      Ha.tau        = -41572,     # activation energy of tau                                (J mol-1)
      Ha.vomax      = 60110,      # activation energy of Vomax                              (J mol-1)
      Hd.rd         = 200000,     # deactivation energy of rd                               (J mol-1)
      Hd.vcmax      = 200000,     # deactivation energy of Vcmax                            (J mol-1)
      Hd.jmax       = 200000,     # deactivation energy of Jmax                             (J mol-1)
      Hd.tpu        = 200000,     # deactivation energy of TPU                              (J mol-1)
      Topt.rd       = 27.56,      # temperature optimum of rd                               (oC)
      Topt.vcmax    = 27.56,      # temperature optimum of Vcmax                            (oC)
      Topt.jmax     = 19.89,      # temperature optimum of Jmax                             (oC)
      Topt.tpu      = 27.56,      # temperature optimum of TPU                              (oC)
      deltaS.rd     = numeric(1), # 
      deltaS.vcmax  = numeric(1), # 
      deltaS.jmax   = numeric(1), #
      deltaS.tpu    = numeric(1), #
      a_deltaS_t.rd     = 490,    # linear temperature response of rd deltaS   
      a_deltaS_t.vcmax  = 668,    # linear temperature response of vcmax deltaS (Kattge & Knorr)  
      a_deltaS_t.jmax   = 660,    # linear temperature response of jmax  deltaS (Kattge & Knorr)
      a_deltaS_t.tpu    = 485,    # linear temperature response of tpu   deltaS (Kattge & Knorr)
      b_deltaS_t.rd     = 0,      # linear temperature response of rd deltaS
      b_deltaS_t.vcmax  = -1.07,  # linear temperature response of vcmax deltaS (Kattge & Knorr)
      b_deltaS_t.jmax   = -0.75,  # linear temperature response of jmax deltaS (Kattge & Knorr)
      b_deltaS_t.tpu    = 0,      # linear temperature response of tpu deltaS (Kattge & Knorr)
      q10.rd        = 2,          # Q10 of Rd                                               (-)
      q10.vcmax     = 2,          # Q10 of Vcmax                                            (-)
      q10.jmax      = 2,          # Q10 of Jmax                                             (-)
      q10.tpu       = 2,          # Q10 of TPU                                              (-)
      q10.Kc        = 2,          # Q10 of Kc                                               (-)
      q10.Ko        = 2,          # Q10 of Ko                                               (-)
      q10.tau       = 0.57,       # Q10 of tau                                              (-)
      a_q10_t.rd    = 3.22,       # linear temperature response of rd Q10 (Tjoelker etal 2001)
      b_q10_t.rd    = -0.046,     # linear temperature response of rd Q10 (Tjoelker etal 2001)
      tupp_cox.vcmax= 36,         # upper leaf T for Vcmax temp scaling from Cox 2001       (oC)
      tupp_cox.rd   = 45,         # upper leaf T for rd temp scaling from Cox 2001 (LM3)    (oC)
      tlow_cox.vcmax= 0,          # lower leaf T for Vcmax temp scaling from Cox 2001       (oC)
      tlow_cox.rd   = 5,          # lower leaf T for rd temp scaling from Cox 2001          (oC)
      exp_cox.vcmax = 0.3,        # exponent for Vcmax temp scaling from Cox 2001           (-)
      exp_cox.rd    = 0.4,        # exponent for rd temp scaling from Cox 2001              (-)
      gstar_bf_a    = 0.012,      # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
      gstar_bf_b    = 1.68,       # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
      gstar_bf_c    = 42.7,       # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
      
      #physical constants
      R   = 8.31446               # molar gas constant                                      (m2 kg s-2 K-1 mol-1  ==  Pa m3 mol-1K-1)
    )
    
    # run control parameters
    cpars <- list(
      verbose       = F,          # write diagnostic output during runtime 
      verbose_loop  = F,          # write diagnostic output on the solver during runtime 
      cverbose      = F,          # write configuration output during runtime 
      output        = 'slim'      # type of output from run function
    )
    
    
    
    ###########################################################################
    # configure & run_met functions

    configure <- function(.,vlist,df,o=T){
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
      
      # name and assign the UQ variables
      uqvars <- names(df)
      prefix <- substr(uqvars,1,str_locate(uqvars,'\\.')[,2]-1)
      modobj <- .$name
      dfss   <- which(prefix==modobj)
      vlss   <- match(uqvars[dfss], paste0(modobj,'.',names(.[[vlist]])) )
      # could write a line to catch NAs in vlss
      .[[vlist]][vlss] <- df[dfss]
      
      if(.$cpars$cverbose&o) {
        print('',quote=F)
        print('Leaf configure:',quote=F)
        print(prefix,quote=F)
        print(df,quote=F)
        print(.$fnames,quote=F)
      }
    }
    
    run_met <- function(.,l){
      # This wrapper function is called from an lapply function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are sequential
      # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
      # typically used to run the model using data collected at a specific site and to compare against observations
      
      # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
      # any "env" variables specified in the "dataf$env" dataframe but also specified in .$dataf$met will be overwritten by the .$dataf$met values 
      
      # met data assignment
      .$configure(vlist='env',df=.$dataf$met[l,],F)

      # run model
      .$run()              
    }
    
    
    
    #######################################################################        
    # Test functions
    # - not copied when the object is cloned

    .test_leaf <- function(.,verbose=T,verbose_loop=T,leaf.par=1000,leaf.ca_conc=300,rs='f_rs_medlyn2011',gd='f_ficks_ci') {
      
      if(verbose) {
        str.proto(.)
        print(.$env)
      }
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- rs
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$fnames$gas_diff    <- gd
      
      .$env$par     <- leaf.par
      .$env$ca_conc <- leaf.ca_conc
      
      .$run()
    }


    .test_solverFunc <- function(.,verbose=T,verbose_loop=T,leaf.par=200,leaf.ca_conc=300,rs='f_rs_medlyn2011') {
      
      if(verbose) {
        str.proto(.)
        print(.$env)
      }
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop

      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- rs
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$fnames$gas_diff    <- 'f_ficks_ci'
      
      .$env$ca_conc        <- leaf.ca_conc
      .$env$par            <- 0
      .$run()

      # proper calc of electron transport rate
      .$env$par            <- leaf.par
      .$state$J <- get(.$fnames$etrans)(.)
      f_A_r_leaf(.,A=-10:100)
    }
    
        
    .test_tscalar <- function(.,leaf.temp=0:50,leaf.par=c(1000),leaf.ca_conc=400,rs='f_rs_medlyn2011',
                              tcor_asc='f_temp_scalar_Arrhenius',tcor_des='f_scalar_none', 
                              verbose=F,verbose_loop=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <- 'full'
      
      if(verbose) str.proto(.)
      
      .$fnames$vcmax_tcor_asc  <- tcor_asc
      .$fnames$vcmax_tcor_des  <- tcor_des
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- rs
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$fnames$solver      <- 'f_R_Brent_solver'
      
      .$dataf     <- list()
      .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par','leaf.temp')))      
      .$dataf$out <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      
      print(cbind(.$dataf$met,.$dataf$out))
      p1 <- xyplot(I(unlist(vcmaxlt)/unlist(vcmax))~leaf.temp|as.factor(paste(tcor_asc,tcor_des)),.$dataf$out,abline=list(h=c(0,1),v=.$pars$reftemp.vcmax),
                   ylab=expression('scalar'),xlab=expression(T*' ['^o*C*']'))
      print(p1)
    }
    
    
    .test_aci <- function(.,leaf.par=c(100,1000),leaf.ca_conc=seq(0.1,1500,50),rs='f_rs_medlyn2011', 
                          verbose=F,verbose_loop=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$pars$diag           <- diag
      # .$cpars$output        <- 'all_lim'
      .$cpars$output        <- 'full'
      
      if(verbose) str.proto(.)
      
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- rs
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$fnames$solver      <- 'f_R_Brent_solver'
      
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
    

    .test_aci_light <- function(.,leaf.par=seq(10,2000,50),leaf.ca_conc=seq(1,1200,50),rs='f_rs_medlyn2011',
                                verbose=F,verbose_loop=F,output=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$pars$diag           <- diag
      .$cpars$output        <- 'all_lim'
      
      if(verbose) str.proto(.)
      
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- rs
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$fnames$solver      <- 'f_R_Brent_solver'
      
      .$dataf     <- list()
      .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))
      
      .$dataf$out <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
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

    .test_aci_analytical <- function(.,rs='f_rs_medlyn2011',leaf.par=c(100,1000),leaf.ca_conc=seq(100,1200,50),leaf.rb=0, 
                                     ana_only=F,verbose=F,verbose_loop=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$pars$diag           <- diag
      .$cpars$output        <- 'all_lim'
      
      if(verbose) str.proto(.)
      
      .$fnames$rs           <- rs
      .$fnames$ri           <- 'f_r_zero'
      .$fnames$rb           <- 'f_rb_constant'
      .$pars$rb             <- leaf.rb
      .$fnames$gas_diff     <- 'f_ficks_ci'
      
      .$dataf     <- list()
      .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))      
      
      if(!ana_only) {
        .$fnames$solver <- 'f_R_Brent_solver'
        .$dataf$out     <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
        num_soln        <- cbind(.$dataf$met,.$dataf$out,sol='num')
        # print(num_soln)
      }

      .$fnames$solver <- 'f_A_r_leaf_analytical'
      .$dataf$out     <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      ana_soln        <- cbind(.$dataf$met,.$dataf$out,sol='ana')
      # print(ana_soln)
      
      odf <- if(ana_only) ana_soln else rbind(num_soln,ana_soln)
      
      if(!ana_only) {
        p1 <- xyplot(A~cc|as.factor(odf$leaf.par),odf,groups=unlist(sol),abline=0,
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
        
        p3 <- xyplot(rs~leaf.ca_conc|as.factor(odf$leaf.par),odf,groups=unlist(sol),abline=0,
                     main=rs,
                     ylab=expression(g[s]*' ['*mol*' '*m^-2*s^-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']')
                     )

        ol <- list(odf,p2,p1,p3)        
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

      # print(p1)
      print(p2)
      ol
    }

    .test_aci_lim <- function(.,rs='f_rs_medlyn2011',et='f_j_farquharwong1984',leaf.par=c(100,1000),leaf.ca_conc=seq(100,1500,50), 
                                     ana_only=F,verbose=F,verbose_loop=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$pars$diag           <- diag
      .$cpars$output        <- 'all_lim'
      
      if(verbose) str.proto(.)
      
      .$fnames$rs           <- rs
      .$fnames$ri           <- 'f_r_zero'
      .$fnames$gas_diff     <- 'f_ficks_ci'
      .$fnames$etrans       <- et
      
      .$dataf     <- list()
      .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))      
      
      .$fnames$Alim <- 'f_lim_farquhar1980'
      .$dataf$out   <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      Fout          <- cbind(.$dataf$met,.$dataf$out,Alim='F1980')
      
      .$fnames$Alim <- 'f_lim_collatz1991'
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
    
    #######################################################################        
    # end object      
})





