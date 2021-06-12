##############################
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
  # VPD from RH (default), or vice versa
  .$rh_or_vpd()
 
  # calculate state parameters
  # photosynthetic parameters
  .super$state_pars$vcmax     <- .$vcmax()
  .super$state_pars$jmax      <- .$jmax()
  .super$state_pars$tpu       <- .$tpu()
  .super$state_pars$k_pepc    <- .$k_pepc()
  .super$state_pars$rd        <- .$rd()
  .super$state_pars$alpha     <- 0.5 * (1-.super$pars$f)

  # kinetic pars & temperature dependence
  # - if the solver involves the energy balance all of this needs to go in the solver
  .super$state_pars$Kc        <- .super$pars$atref[['Kc']] * .[['tcor_asc.Kc']](.,'Kc')
  .super$state_pars$Ko        <- .super$pars$atref[['Ko']] * .[['tcor_asc.Ko']](.,'Ko')
  .super$state_pars$Km        <- .super$state_pars$Kc * (1+(.super$state$oi/.super$state_pars$Ko))
  .super$state_pars$gstar     <- .$gstar()
  .super$state_pars$vcmaxlt   <- .super$state_pars$vcmax  * .[['tcor_asc.vcmax']](.,'vcmax')   * .[['tcor_des.vcmax']](.,'vcmax')
  .super$state_pars$jmaxlt    <- .super$state_pars$jmax   * .[['tcor_asc.jmax']](.,'jmax')     * .[['tcor_des.jmax']](.,'jmax')
  .super$state_pars$tpult     <- .super$state_pars$tpu    * .[['tcor_dep.tpu']](.,'tpu')
  .super$state_pars$k_pepc_lt <- .super$state_pars$k_pepc * .[['tcor_asc.k_pepc']](.,'k_pepc') * .[['tcor_des.k_pepc']](.,'k_pepc')

  # conductance/resistance terms
  # - if either of these functions become a function of co2 or assimilation they can be easily moved into the solver
  .super$state_pars$rb <- .$rb()
  .super$state_pars$ri <- .$ri()

  # print state parameters to screen
  if(.super$cpars$verbose) {
    print(.super$state_pars)
  }

  # calculate physiological state
  # electron transport rate
  .super$state$J  <- .$etrans()
  # respiration
  .super$state$rd <- .super$state_pars$rd * .[['tcor_dep.rd']](.,'rd')

  # if PAR > 0, run photosynthesis
  if(.super$env$par > 0) {

    # account for decreased respiration in the light
    .super$state$rd         <- .$rl_rd() * .super$state$rd
    .super$state_pars$gamma <- (-.super$state_pars$vcmaxlt * .super$state_pars$gstar - .super$state$rd * .super$state_pars$Km) / (.super$state$rd - .super$state_pars$vcmaxlt)

    # calculate assimilation with no resistence term 
    if(.super$cpars$diag) .super$state$A_noR <- .$fns$assim_no_resistance()

    # calculate assimilation
    solve_analytically <- (.$fnames$rs=='f_rs_constant') & (.super$state_pars$ri==0) & (.super$state_pars$rb==0) 
    if(solve_analytically) .$solver <- f_solver_analytical_leaf_quad
    .super$state$A <- .$solver()

    # after the fact calculations
    # APW: need to handle cc calculation for C4 plants
    #      just set to a constant high value somewhere, maybe could add cc for PEPC limitation but maybe not necessary  
    if(!grepl('analytical',.super$fnames$solver)) {

      if(grepl('c3',.super$fnames$assimilation)) {
        # calculate Ag for each limiting process
        # APW: something is not working here, or elsewhere, to calculate these gross rates 
        .super$state$Acg[] <- .super$state$Acg * .super$state$cc
        .super$state$Ajg[] <- .super$state$Ajg * .super$state$cc
        .super$state$Apg[] <- .super$state$Apg * .super$state$cc
      } 

      # calculate intermediate state variables
      .super$state$cb      <- .$gas_diff(A=.super$state$A, r=1.4*.super$state_pars$rb, c=.super$state$ca )
      .super$state_pars$rs <- .$rs()
      .super$state$ci      <- .$gas_diff(A=.super$state$A, r=1.6*.super$state_pars$rs, c=.super$state$cb)

      ### should probably occur outside of !is analytical case statement
      # if rs is negative (occurs when A is negative) recalculate with fixed rs at 1/g0
      if(!is.na(.super$state$A)) {
        if( .super$state$A<0 | .super$state$cc<0 | .super$state_pars$rs<0 ) {
           
          # perhaps write this flag into the leaf object
          #print(paste('solver returned negative value of A, cc, or rs; recalculate assuming rs = r0'))

          # temporarily reassign solver function
          .$solver <- if(grepl('c3',.super$fnames$assimilation)) f_solver_analytical_leaf_quad_r0
                      else                                       f_solver_analytical_leaf_c4_r0

          # calculate assimilation
          .super$state$A <- .$solver()

          # reassign original solver
          .$solver <- get(.super$fnames$solver)
      }}
    }

    # reassign original solver
    if(solve_analytically) .$solver <- get(.super$fnames$solver)

    # assign the limitation state a numerical code - assumes the minimum is the primary limiting rate
    lim              <- min(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg), na.rm=T )
    lim_ss           <- which(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg)==lim)
    .super$state$lim <- c(2,3,7)[lim_ss]

    # diagnostic calculations
    if(.super$cpars$diag & grepl('c3',.super$fnames$assimilation)) { 
      .super$state$transition       <- .$fns$transition_cc()
      .super$state$photorespiration <- .$fns$photorespiration()
    }

  # if PAR <= 0
  } else {
    # assume infinite conductances when concentration gradient is small
    # - this ignores the build up of CO2 within the leaf due to respiration and high rs
    .super$state$cc <- .super$state$ci <- .super$state$cb <- .super$state$ca
    .super$state$A       <- -.super$state$rd
    .super$state$lim     <- 0
    .super$state_pars$rs <- .$rs()
    if(.super$cpars$diag) { 
      .super$state$A_noR            <- NA
      .super$state$transition       <- NA
      .super$state$photorespiration <- NA
    }
  }
}



### END ###
