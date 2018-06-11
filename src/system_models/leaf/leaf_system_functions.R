################################
#
# Leaf system functions 
# 
# AWalker February 2018
#
################################

source('../environment_functions.R')



################################
# General enzyme kinetic model of photosynthesis
# based on Farquhar et al. 1980 
f_leafsys_enzymek <- function(.) {

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
  # - if the solver involves the energy balance all of this needs to go in the solver
  .$state_pars$Kc      <- .$pars$atref[['Kc']] * get(.$fnames$Kc_tcor)(., var='Kc' )
  .$state_pars$Ko      <- .$pars$atref[['Ko']] * get(.$fnames$Ko_tcor)(., var='Ko' ) 
  .$state_pars$Km      <- .$state_pars$Kc * (1+(.$state$oi/.$state_pars$Ko)) 
  .$state_pars$gstar   <- get(.$fnames$gstar)(.) 
  .$state_pars$vcmaxlt <- .$state_pars$vcmax * get(.$fnames$vcmax_tcor_asc)(., var='vcmax' ) * get(.$fnames$vcmax_tcor_des)(., var='vcmax' )
  .$state_pars$jmaxlt  <- .$state_pars$jmax  * get(.$fnames$jmax_tcor_asc)(., var='jmax' )  * get(.$fnames$jmax_tcor_des)(., var='jmax' )
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
    # run photosynthesis
    # account for decreased respiration in the light
    .$state$respiration  <- get(.$fnames$rl_rd_scalar)(.) * .$state$respiration
    .$state_pars$gamma   <- (-.$state_pars$vcmaxlt * .$state_pars$gstar - .$state$respiration * .$state_pars$Km) / (.$state$respiration - .$state_pars$vcmaxlt)
  
    # diagnostic calculations
    if(.$pars$diag) {
      .$state$A_noR      <- f_A_r_leaf_noR(.)
      .$state$transition <- transition_cc(.)
    }
  
    # calculate assimilation 
    .$state$A       <- get(.$fnames$solver)(.)      
    # assign the limitation state a numerical code - assumes the minimum is the dominant limiting rate
    .$state$lim     <- c(2,3,7)[which(c(.$state$Acg,.$state$Ajg,.$state$Apg)==min(c(.$state$Acg,.$state$Ajg,.$state$Apg),na.rm=T))]       
 
    # after the fact calculations
    if(!grepl('analytical',.$fnames$solver)) {
      # calculate Ag for each limiting process
      .$state$Acg     <- .$state$Acg * .$state$cc 
      .$state$Ajg     <- .$state$Ajg * .$state$cc 
      .$state$Apg     <- .$state$Apg * .$state$cc
 
      # calculate intermediate state variables 
      .$state$cb      <- f_ficks_ci(.,A=.$state$A,r=1.4*.$state_pars$rb,c=.$state$ca)
      .$state_pars$rs <- get(.$fnames$rs)(.) 
      .$state$ci      <- f_ficks_ci(.,A=.$state$A,r=1.6*.$state_pars$rs,c=.$state$cb)
    
 
      # if rs is negative (occurs when A is negative) recalculate with fixed rs at 1/g0 
      if( .$state$A<0 | .$state$cc<0 | .$state_pars$rs<0 ) {

        # perhaps write this flag into the leaf object
        #print(paste('negative values loop'))

        # temporarily rename solver function 
        solver <- .$fnames$solver
        .$fnames$solver <- 'f_A_r0_leaf_analytical_quad'
      
        # calculate assimilation 
        .$state$A       <- get(.$fnames$solver)(.)      
        # assign the limitation state a numerical code - assumes the minimum is the dominant limiting rate
        .$state$lim     <- c(2,3,7)[which(c(.$state$Acg,.$state$Ajg,.$state$Apg)==min(c(.$state$Acg,.$state$Ajg,.$state$Apg),na.rm=T))]       
   
        #.$fnames$solver_func <- solver 
        .$fnames$solver <- solver 
    }}
 
  }
 
  # if PAR <= 0
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
}



### END ###
