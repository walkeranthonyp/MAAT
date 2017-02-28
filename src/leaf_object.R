################################
#
# Leaf object 
# 
# AWalker March 2014
#
################################

library(proto)
library(stringr)
library(plyr)

source('leaf_functions.R')



# LEAF OBJECT
###############################################################################

leaf_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'leaf'
    
    # no expected child objects
    
    # build function
    .build <- function(.) {
      as.proto(.$as.list())
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.,verbose=F) {
      
      if(verbose) .$cpars$verbose
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
      
      
      # calculate state parameters
      # photosynthetic parameters
      .$state_pars$vcmax   <- get(.$fnames$vcmax)(.)
      .$state_pars$jmax    <- get(.$fnames$jmax)(.)
      .$state_pars$tpu     <- get(.$fnames$tpu)(.) # needs work, how is TPU normally set?
      .$state_pars$alpha   <- 0.5 * (1-.$pars$f)
      
      # kinetic pars & temperature dependence
      # - if the solver involves the energy balance all of this crap needs to go in the solver! 
      .$state$leaf_temp    <- .$env$temp               
      .$state_pars$Kc      <- .$pars$atref.Kc * get(.$fnames$Kc_tcor)(.,parlist=list(Ha=.$pars$Ha.Kc))
      .$state_pars$Ko      <- .$pars$atref.Ko * get(.$fnames$Ko_tcor)(.,parlist=list(Ha=.$pars$Ha.Ko)) 
      .$state_pars$Km      <- .$state_pars$Kc * (1+(.$state$oi/.$state_pars$Ko)) 
      .$state_pars$gstar   <- get(.$fnames$gstar)(.) 
      .$state_pars$gamma   <- (-.$state_pars$vcmaxlt * .$state_pars$gstar - .$state$respiration * .$state_pars$Km) / (.$state$respiration - .$state_pars$vcmaxlt)
      .$state_pars$vcmaxlt <- .$state_pars$vcmax * get(.$fnames$vcmax_tcor)(.,parlist=list(Ha=.$pars$Ha.vcmax,Hd=.$pars$Hd.vcmax,Topt=.$pars$Topt.vcmax,q10=.$pars$q10.vcmax,
                                                                                           .$pars$tupp=tupp_cox.vcmax,.$pars$tlow=tlow_cox.vcmax))
      .$state_pars$jmaxlt  <- .$state_pars$jmax  * get(.$fnames$jmax_tcor)(.,parlist=list(Ha=.$pars$Ha.jmax,Hd=.$pars$Hd.jmax,Topt=.$pars$Topt.jmax,q10=.$pars$q10.jmax))
      .$state_pars$tpult   <- .$state_pars$tpu   * get(.$fnames$tpu_tcor)(.,parlist=list(Ha=.$pars$Ha.tpu,Hd=.$pars$Hd.tpu,Topt=.$pars$Topt.tpu,q10=.$pars$q10.tpu))
      
      # conductance/resistance terms
      # - if either of these functions become a function of co2 or assimilation they can be easily moved into the solver
      .$state_pars$rb      <- get(.$fnames$rb)(.)
      .$state_pars$ri      <- get(.$fnames$ri)(.)  
      
      # print state parameters to screen
      if(.$cpars$verbose) {
        print(.$state_pars)
      }  
      
      # calculate physiological state
      # respiration
      .$state$respiration <- get(.$fnames$respiration)(.)
      # if PAR > 0
      if(.$env$par > 0) {
        # diagnostic calculations
        if(.$cpars$diag) {
          # assume infinite conductances to initialise solver
          .$state$cc <- .$state$ci <- .$state$cb <- .$state$ca
          .$state$A_noR      <- f_A_r_leaf_noR(.)
          .$state$transition <- transition_cc(.)
        }
        # run photosynthesis
        # assume infinite conductances to initialise solver
        .$state$cc <- .$state$ci <- .$state$cb <- .$state$ca
        # determine rate limiting step - this is done based on carboxylation, not net assimilation (Gu etal 2010).
        .$state$A       <- get(.$fnames$solver)(.)      
        # assign the limitation state - assumes the minimum is the dominant limiting rate
        .$state$lim     <- c('wc','wj','wp')[which(c(.$state$wc,.$state$wj,.$state$wp)==min(c(.$state$wc,.$state$wj,.$state$wp),na.rm=T))]       
        # after the fact calculations
        .$state$cb      <- f_ficks_ci(.,A=.$state$A,r=.$state_pars$rb,c=.$state$ca)
        if(!grepl('analytical',.$fnames$solver)) .$state_pars$rs <- get(.$fnames$rs)(.) 
        .$state$ci      <- f_ficks_ci(.,A=.$state$A,r=.$state_pars$rs,c=.$state$cb)
        .$state$cc      <- f_ficks_ci(.,A=.$state$A,r=.$state_pars$ri,c=.$state$ci)        
      }
      # if PAR < 0
      else {
        # assume infinite conductances when concentration gradient is small
        # - this ignores the build up of CO2 within the leaf due to respiration and potentially high rs 
        .$state$cc <- .$state$ci <- .$state$cb <- .$state$ca
        .$state$A_noR      <- NA
        .$state$transition <- NA
        .$state$A          <- -.$state$respiration      
        .$state$lim        <- 'night'       
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

    #output processing function
    # -- returns a list of outputs
    output <- function(.){
      if(.$cpars$output=='run') {
        lout <- 
          list(A=.$state$A,cc=.$state$cc,ci=.$state$ci,
               gi=1/.$state_pars$ri,gs=1/.$state_pars$rs,gb=1/.$state_pars$rb,
               respiration=.$state$respiration,lim=.$state$lim)
        
      } else if(.$cpars$output=='all_lim') {
        lout <- 
          list(A=.$state$A,wc=.$state$wc,wj=.$state$wj,wp=.$state$wp,
               cc=.$state$cc,ci=.$state$ci,ca=.$state$ca,
               gi=1/.$state_pars$ri,gs=1/.$state_pars$rs,gb=1/.$state_pars$rb,
               respiration=.$state$respiration,lim=.$state$lim)     
        
      } else if(.$cpars$output=='full') {
        lout <- c(.$state,.$state_pars)
        
      } else if(.$cpars$output=='sphagnum') {
        lout <-
          list(A=.$state$A,cc=.$state$cc,ci=.$state$ci,
               gi=1/.$state_pars$ri,gs=1/.$state_pars$rs,gb=1/.$state_pars$rb,
               respiration=.$state$respiration,lim=.$state$lim,fwdw=.$state$fwdw_ratio) 
        
      }
      
      if(.$cpars$diag) c( lout, list(A_noR=.$state$A_noR,transition=.$state$transition) ) 
      else                lout
      
    }    
    
    
    
    ###########################################################################
    # Variables etc
    
    # function names
    fnames <- list(
      gstar       = 'f_gstar_constref',
      gstar_tcor  = 'f_temp_scalar_quadratic_bf1985',
      tau_tcor    = 'f_temp_scalar_Q10',
      Kc_tcor     = 'f_temp_scalar_Arrhenius',
      Ko_tcor     = 'f_temp_scalar_Arrhenius',
      vcmax       = 'f_vcmax_lin',
      jmax        = 'f_jmax_walker2014',
      tpu         = 'f_constant_tpu',
      vcmax_tcor  = 'f_temp_scalar_modArrhenius',
      jmax_tcor   = 'f_temp_scalar_modArrhenius',
      tpu_tcor    = 'f_temp_scalar_modArrhenius',
      etrans      = 'f_j_harley1992',
      wc          = 'f_wc_farquhar1980',
      wj          = 'f_wj_generic',
      wp          = 'f_wp_vonc2000',            # 'f_wp_collatz1991'
      gas_diff    = 'f_ficks_ci_bound0',
      respiration = 'f_rd_collatz1991',
      fwdw_ratio  = 'f_none',                   # 'f_fwdw_wl_lin' 'f_fwdw_wl_exp'
      cica_ratio  = 'f_cica_constant',             
      ri          = 'f_r_zero',
      rs          = 'f_r_zero',
      rb          = 'f_r_zero',
      solver      = 'f_R_Brent_solver',
      solver_func = 'f_A_r_leaf',
      Alim        = 'f_lim_farquhar1980'
    )
    
    # leaf environment
    env <- list(
      ca_conc   = numeric(0),          # (umol mol-1)
      o2_conc   = 0.21,                # ( mol mol-1)
      par       = numeric(0),          # (umol photons m-2 s-1)
      water_l   = numeric(0),          # (mm) water level relative to hollow surfwce
      sphag_l   = 0,                   # (mm) Sphagnum surface relative to hollow surfwce
      temp      = 25,                  # (oC)
      vpd       = 2,                   # (kPa)
      rh        = 0.8,                 # (unitless - proportion)
      atm_press = 101325               # ( Pa)
      )

    # leaf state
    state <- list(
      #environmental state
      oi = numeric(0),                 # atmospheric & internal O2  (kPa)
      ca = numeric(0),                 # atmospheric CO2            ( Pa)
      cb = numeric(0),                 # boundary layer CO2         ( Pa)
      ci = numeric(0),                 # leaf internal CO2          ( Pa) 
      cc = numeric(0),                 # chloroplast CO2            ( Pa)
      leaf_temp = numeric(0),          # leaf temperature           (oC)
      
      #leaf state - calculated by canopy object so need initialisation
      leafN_area = 2,                  # leaf N per unit area       (g N m-2)
      fwdw_ratio = 5,                  # fresh weight dry weight ratio, used for Sphagnum conductance term 
      
      #calculated state
      J  = numeric(0),                 # electron transport rate                            (umol electrons m-2 s-1) 
      wc = numeric(0),                 # Carboxylaton limited rate of net asssimilation     (umol m-2 s-1)
      wj = numeric(0),                 # light limited rate of carboxylation                (umol m-2 s-1)
      wp = numeric(0),                 # TPU limited rate of carboxylation                  (umol m-2 s-1)
      A            = numeric(0),       # actual rate of carboxylation                       (umol m-2 s-1)
      respiration  = numeric(0),       # actual rate of respiration                         (umol m-2 s-1)
      lim          = character(0),     # flag indicationg limitation state of assimilation, wc = wc limited, wj = wj limited, wp = wp limited

      # diagnostic state
      A_noR        = numeric(0),       # rate of carboxylation assuming zero resistance to CO2 diffusion (umol m-2 s-1)
      transition   = numeric(0)        # cc at the transition point where wc = wj                        (Pa)
    )
    
    # results from solver
    solver_out = NULL
    
    #leaf state parameters (i.e. calculated parameters)
    state_pars <- list(
      vcmax    = numeric(0),   # umol m-2 s-1
      vcmaxlt  = numeric(0),   # umol m-2 s-1
      jmax     = numeric(0),   # umol m-2 s-1
      jmaxlt   = numeric(0),   # umol m-2 s-1
      tpu      = numeric(0),   # umol m-2 s-1
      Kc       = numeric(0),   #  Pa
      Ko       = numeric(0),   # kPa
      Km       = numeric(0),   #  Pa
      gstar    = numeric(0),   #  Pa
      tau      = numeric(0),   #  -
      gamma    = numeric(0),   #  Pa
      rb       = numeric(0),   # m2s mol-1 
      rs       = numeric(0),   # m2s mol-1 
      ri       = numeric(0),   # m2s mol-1     
      alpha    = numeric(0),   # mol electrons mol-1 absorbed photosynthetically active photons
      cica_chi = numeric(0)    # Ci:Ca ratio 
    )
    
    #leaf parameters
    pars   <- list(
      # photosynthetic parameters
      # deprecated    alpha    = 0.24,         # harley 1992 alpha - Williams & Flannagan 1998 use 0.21 but calculate 0.25 
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
      g1_medlyn     = 6,          # Medlyn 2011 gs slope                                   (kPa^0.5)
      g1_leuning    = 10,         # Leuning 1995 gs slope                                  (unitless - likely higher than medlyn and ball g1)
      d0            = 1,          # Leuning 1995 D0                                        (kPa)
      g1_ball       = 6,          # Ball 1987 gs slope                                     (unitless - multiplier on RH as a proportion)
      rs            = 1/0.15,     # stomatal resistance                                    (m2s mol-1 h2o)
      gi            = 0.15,       # mesophyll conductance                                  (molm-2s-1 - expressed in these units for consistency with other conductance terms, often expressed in the literature per unit Pa)
      cica_chi      = 0.7,        # constant Ci:Ca ratio                                   (unitless)
      ri            = 1/0.15,     # mesophyll resistance                                   (m2s mol-1 - expressed in these units for consistency with other resistance terms, often expressed in the literature multiplied by Pa)
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
      reftemp.tau   = 25,         # reference temperature at which tau scalar = 1           (oC)
      atref.rd      = 2,          # rd at ref temp (usually 25oC)    - used to set rd as a parameter                        (umolm-2s-1) 
      atref.vcmax   = 50,         # vcmax at ref temp (usually 25oC) - used to set Vcmax as a parameter instead of an f(N)  (umolm-2s-1) 
      atref.jmax    = 100,        # jmax at ref temp (usually 25oC)  - used to set Jmax as a parameter instead of an f(N)   (umolm-2s-1)
      atref.tpu     = 10,         # tpu at ref temp (usually 25oC)   - used to set TPU as a parameter                       (umolm-2s-1)
      atref.Kc      = 40.49,      # Kc for RuBisCO at ref temp (usually 25oC)               ( Pa)
      atref.Ko      = 27.84,      # Kc for RuBisCO at ref temp (usually 25oC)               (kPa)
      atref.gstar   = 4.325,      # Gamma star at ref temp (usually 25oC), 4.325 is Farquhar & Brooks value converted to Pa (Pa)
      atref.tau     = 2600,       # CO2/O2 specificity ratio at ref temp (usually 25oC), Collatz 1991 (-)
      atref.vomax   = numeric(0),
      Ha.vcmax      = 69830,      # activation energy of Vcmax                              (J mol-1)
      Ha.jmax       = 100280,     # activation energy of Jmax                               (J mol-1)
      Ha.tpu        = 69830,      # activation energy of TPU                                (J mol-1)
      Ha.Kc         = 79430,      # activation energy of Kc                                 (J mol-1)
      Ha.Ko         = 36380,      # activation energy of Ko                                 (J mol-1)
      Ha.gstar      = 37830,      # activation energy of gamma star                         (J mol-1)
      Ha.vomax      = 60110,      # activation energy of Vomax                              (J mol-1)
      Hd.vcmax      = 200000,     # deactivation energy of Vcmax                            (J mol-1)
      Hd.jmax       = 200000,     # deactivation energy of Jmax                             (J mol-1)
      Hd.tpu        = 200000,     # deactivation energy of TPU                              (J mol-1)
      Topt.vcmax    = 27.56,      # temperature optimum of Vcmax                            (oC)
      Topt.jmax     = 19.89,      # temperature optimum of Jmax                             (oC)
      Topt.tpu      = 27.56,      # temperature optimum of TPU                              (oC)
      deltaS.vcmax  = numeric(0), # 
      deltaS.jmax   = numeric(0), #
      deltaS.tpu    = numeric(0), #
      q10.vcmax     = 2,          # Q10 of Vcmax                                            (-)
      q10.jmax      = 2,          # Q10 of Jmax                                             (-)
      q10.tpu       = 2,          # Q10 of TPU                                              (-)
      q10.tau       = 0.57,       # Q10 of tau                                              (-)
      tupp_cox.vcmax= 36,         # upper leaf T for Vcmax temp scaling from Cox 1991       (-)
      tlow_cox.vcmax= 0,          # lower leaf T for Vcmax temp scaling from Cox 1991       (-)
      gstar_jo_a    = numeric(0), # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
      gstar_jo_b    = numeric(0), # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
      gstar_jo_c    = numeric(0), # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
      
      #physical constants
      R   = 8.31446               # molar gas constant                                      (m2 kg s-2 K-1 mol-1  ==  Pa m3 mol-1K-1)
    )
    
    # run control parameters
    cpars <- list(
      verbose       = F,          # write diagnostic output during runtime 
      verbose_loop  = F,          # write diagnostic output on the solver during runtime 
      cverbose      = F,          # write configuration output during runtime 
      diag          = F,          # calculate diagnostic output during runtime and add to output, such as cc transition point and non-stomatal limited assimilation rate 
      output        = 'run'       # type of output from run function
    )
    
    
    
    ###########################################################################
    # Run & configure functions
    
    # initialisation function
    init <- function(.) {
      # when run this function expects this model object to be cloned within the wrapper object, 
      # and for the wrapper object to be named maat
      # the 'init' function resides in the model as initialisation depends on the hierarchical structure of the model (e.g. leaf within canopy )
      
      # this function prefixes the names of a list with 'leaf.' 
      comb_init_list <- function(.,lls) {

        if(sum(is.na(lls))==length(lls)) NA
        else {
          names(lls) <- paste('leaf',names(lls),sep='.') # need to write error message in here
          lls
        }
      }
            
      maat$static$fnames <- comb_init_list(lls=.$init_static$leaf$fnames)
      maat$static$pars   <- comb_init_list(lls=.$init_static$leaf$pars)
      maat$static$env    <- comb_init_list(lls=.$init_static$leaf$env)
      
      maat$vars$fnames   <- comb_init_list(lls=.$init_dynamic$leaf$fnames)
      maat$vars$pars     <- comb_init_list(lls=.$init_dynamic$leaf$pars)
      maat$vars$env      <- comb_init_list(lls=.$init_dynamic$leaf$env)
      
      if(.$wpars$UQ) {
        if(.$wpars$UQtype=='ye') maat$vars$pars_proc <- comb_init_list(lls=.$init_dynamic$leaf$pars_proc)
        if(.$wpars$eval_strings) maat$vars$pars_eval <- comb_init_list(lls=.$init_dynamic$leaf$pars_eval)
      }
    }
    
    configure <- function(.,func,df,o=T){
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
      
#       print(func)
#       print(df)
      
      # name and assign the UQ variables
      uqvars     <- names(df)
      prefix     <- substr(uqvars,1,str_locate(uqvars,'\\.')[,2]-1)
      lapply(uqvars[which(prefix=='leaf')], func, .=., df=df)

#       print(func)
      
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
      .$configure(func='write_env',df=.$dataf$met[l,],F)
      
      # run model
      .$run()              
    }
    
    # variable assignment functions - called from the above configuration function
    write_fnames <- function(.,var,df) .$fnames[which(names(.$fnames)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))] <- df[which(names(df)==var)]
    write_pars   <- function(.,var,df) .$pars[which(names(.$pars)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))]     <- df[which(names(df)==var)]
    write_env    <- function(.,var,df) .$env[which(names(.$env)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))]       <- df[which(names(df)==var)]
    write_state  <- function(.,var,df) .$state[which(names(.$state)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))]   <- df[which(names(df)==var)]
    
    
    
    #######################################################################        
    # Test functions
    # - not copied when the object is cloned

    .test_leaf <- function(.,verbose=T,verbose_loop=T,leaf.par=1000,leaf.ca_conc=300){
      
      if(verbose) {
        str.proto(.)
        print(.$env)
      }
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- 'f_ri_constant'
      .$fnames$solver_func <- 'f_A_r_leaf'
      
      .$env$par     <- leaf.par
      .$env$ca_conc <- leaf.ca_conc
      
      .$run()
    }
    
    .test_aci <- function(.,leaf.par=c(100,1000),leaf.ca_conc=seq(0.1,1200,50), 
                          verbose=F,verbose_loop=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$diag          <- diag

      if(verbose) str.proto(.)
      
      .$fnames$ri          <- 'f_r_zero'
#       .$fnames$rs          <- 'f_ri_constant'
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$pars$output        <- 'all_lim'
      
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
    

    .test_aci_light <- function(.,leaf.par=seq(10,2000,50),leaf.ca_conc=seq(1,1200,50),
                                verbose=F,verbose_loop=F,output=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$diag          <- diag

      if(verbose) str.proto(.)
      
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- 'f_ri_constant'
      .$fnames$solver_func <- 'f_A_r_leaf'
      .$cpars$output       <- 'all_lim'
      
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

    .test_aci_analytical <- function(.,rs='f_r_zero',leaf.par=c(100,1000),leaf.ca_conc=seq(100,1200,50), 
                                     ana_only=F,verbose=F,verbose_loop=F,diag=F) {
      
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$diag          <- diag
      
      if(verbose) str.proto(.)
      
      .$fnames$rs           <- rs
      .$fnames$ri           <- 'f_r_zero'
      .$fnames$gas_diff     <- 'f_ficks_ci'
      .$pars$output         <- 'all_lim'
      
      .$dataf     <- list()
      .$dataf$met <- expand.grid(mget(c('leaf.ca_conc','leaf.par')))      
      
      if(!ana_only) {
        .$fnames$solver <- 'f_R_Brent_solver'
        .$dataf$out <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
        num_soln <- cbind(.$dataf$met,.$dataf$out,sol='num')
        # print(num_soln)
      }

      .$fnames$solver <- 'f_A_r_leaf_analytical'
      .$dataf$out <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      ana_soln <- cbind(.$dataf$met,.$dataf$out,sol='ana')
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
        
        p2 <- xyplot(A~leaf.ca_conc|as.factor(odf$leaf.par),odf,groups=unlist(sol),abline=0,
                     main=rs,
                     ylab=expression('A ['*mu*mol*' '*m^-2*s^-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'),
                     panel=function(subscripts=subscripts,...) {
                       if(diag) {
                         panel.abline(v=odf$transition[subscripts][1])
                         panel.points(y=odf$A_noR[subscripts],x=odf$cc[subscripts],col='black')                       
                       }
                       panel.xyplot(subscripts=subscripts,...)
                     })
        
        p3 <- xyplot(gs~leaf.ca_conc|as.factor(odf$leaf.par),odf,groups=unlist(sol),abline=0,
                     main=rs,
                     ylab=expression(g[s]*' ['*mol*' '*m^-2*s^-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'),
                     )
        
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
        
      }

      # print(p1)
      print(p2)

      list(odf,p2,p1,p3)
    }

    #######################################################################        
    # end object      
})





