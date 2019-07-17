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
f_sys_enzymek <- function(.) {

  # calculate environmental state
  # CO2 partial pressure (Pa)
  .super$state$ca <- .super$env$ca_conc * .super$env$atm_press * 1e-6
  # O2 partial pressure  (kPa)
  .super$state$oi <- .super$env$o2_conc * .super$env$atm_press * 1e-3
  # leaf temperature
  .super$state$leaf_temp <- .super$env$temp
  # RH from VPD
  .super$env$rh   <- f_rh_from_vpd(.)

  # calculate state parameters
  # photosynthetic parameters
  .super$state_pars$vcmax   <- .$vcmax()
  .super$state_pars$jmax    <- .$jmax()
  .super$state_pars$tpu     <- .$tpu()
  .super$state_pars$rd      <- .$rd()
  .super$state_pars$alpha   <- 0.5 * (1-.super$pars$f)

  # kinetic pars & temperature dependence
  # - if the solver involves the energy balance all of this needs to go in the solver
  .super$state_pars$Kc      <- .super$pars$atref[['Kc']] * .[['tcor_asc.Kc']](.,'Kc')
  .super$state_pars$Ko      <- .super$pars$atref[['Ko']] * .[['tcor_asc.Ko']](.,'Ko')
  .super$state_pars$Km      <- .super$state_pars$Kc * (1+(.super$state$oi/.super$state_pars$Ko))
  .super$state_pars$gstar   <- .$gstar()
  .super$state_pars$vcmaxlt <- .super$state_pars$vcmax * .[['tcor_asc.vcmax']](.,'vcmax') * .[['tcor_des.vcmax']](.,'vcmax')
  .super$state_pars$jmaxlt  <- .super$state_pars$jmax  * .[['tcor_asc.jmax']](.,'jmax')   * .[['tcor_des.jmax']](.,'jmax')
  .super$state_pars$tpult   <- .super$state_pars$tpu   * .[['tcor_dep.tpu']](.,'tpu')

  # conductance/resistance terms
  # - if either of these functions become a function of co2 or assimilation they can be easily moved into the solver
  .super$state_pars$rb      <- .$rb()
  .super$state_pars$ri      <- .$ri()

  # print state parameters to screen
  if(.super$cpars$verbose) {
    print(.super$state_pars)
  }

  # calculate physiological state
  # electron transport rate
  .super$state$J <- .$etrans()
  # respiration
  .super$state$rd  <- .super$state_pars$rd * .[['tcor_dep.rd']](.,'rd')

  # if PAR > 0
  if(.super$env$par > 0) {
    # run photosynthesis
    # account for decreased respiration in the light
    .super$state$rd  <- .$rl_rd_scalar() * .super$state$rd
    .super$state_pars$gamma   <- (-.super$state_pars$vcmaxlt * .super$state_pars$gstar - .super$state$rd * .super$state_pars$Km) / (.super$state$rd - .super$state_pars$vcmaxlt)

    # diagnostic calculations
    if(.super$pars$diag) {
      # these need assigning in fns
      .super$state$A_noR      <- .$fns$assim_no_resistance()
      .super$state$transition <- .$fns$transition_cc()
    }

    # calculate assimilation
    if(.$fnames$rs=='f_rs_constant') {
      # temporarily reassign solver function that calculates A when rs is a constant
      #solver <- .$fnames$solver
      .$fns$solver <- f_A_r0_leaf_analytical_quad
    }
    .super$state$A       <- .$solver()

    # assign the limitation state a numerical code - assumes the minimum is the dominant limiting rate
    .super$state$lim     <- c(2,3,7)[which(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg)==min(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg),na.rm=T))]

    # after the fact calculations
    if(!grepl('analytical',.$fnames$solver)) {
      # calculate Ag for each limiting process
      .super$state$Acg     <- .super$state$Acg * .super$state$cc
      .super$state$Ajg     <- .super$state$Ajg * .super$state$cc
      .super$state$Apg     <- .super$state$Apg * .super$state$cc

      # calculate intermediate state variables
      .super$state$cb      <- .$gas_diff(A=.super$state$A, r=1.4*.super$state_pars$rb, c=.super$state$ca )
      .super$state_pars$rs <- .$rs()
      .super$state$ci      <- .$gas_diff(A=.super$state$A, r=1.6*.super$state_pars$rs, c=.super$state$cb)

      # if rs is negative (occurs when A is negative) recalculate with fixed rs at 1/g0
      if(!is.na(.super$state$A)) {
        if( .super$state$A<0 | .super$state$cc<0 | .super$state_pars$rs<0 ) {

          # debugging: temporarily comment out print statement to make reading log file easier
          # perhaps write this flag into the leaf object
          # print(paste('solver returned negative value of A, cc, or rs; recalculate assuming rs = r0'))

          # temporarily reassign solver function
          #solver <- .super$fnames$solver
          #.$solver <- get('f_A_r0_leaf_analytical_quad')
          .$solver <- f_A_r0_leaf_analytical_quad

          # calculate assimilation
          .super$state$A       <- .$solver()
          # assign the limitation state a numerical code - assumes the minimum is the dominant limiting rate
          .super$state$lim     <- c(3,7)[which(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg)==min(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg),na.rm=T))]

          #.super$fnames$solver_func <- solver
          .$solver <- get(.super$fnames$solver)
    }}}

    if(.$fnames$rs=='f_rs_constant') .$fns$solver <- get(.$fnames$solver)
  }

  # if PAR <= 0
  else {
    # assume infinite conductances when concentration gradient is small
    # - this ignores the build up of CO2 within the leaf due to respiration and high rs
    .super$state$cc <- .super$state$ci <- .super$state$cb <- .super$state$ca
    .super$state$A_noR      <- NA
    .super$state$transition <- NA
    .super$state$A          <- -.super$state$rd
    .super$state$lim        <- 0
    .super$state_pars$rs    <- .$rs()
  }
}



### END ###
