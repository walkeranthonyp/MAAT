################################
#
# MAAT soil_decomp system representation functions (SRFs) 
# 
# AWalker, Matt Craig, October 2019 
#
################################



################################
# soil_decomp system function 

# a system of n pools, composed of a single dynamic solver function
f_sys_npools <- function(.) {

  # calculate state parameters
  .$calc_state_pars()

  # run solver
  .super$state_pars$solver_out <- 
    tryCatch(
      .$solver(.super$state$cpools, c(0L,1L), .$solver_func ),
      error = function(c) {
        print(c)
        matrix(ncol=.super$pars$n_pools+1)
      })
  
  # assign solved values to state
  .super$state$cpools[,1] <- .super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],2:(.super$pars$n_pools+1)]
}  
  
  
  
# a system of n pools, composed of a single dynamic solver function
# - adds loss fluxes into solve, designed to work with f_solver_func_soilR_outfluxes_lossfluxes 
f_sys_npools_lossfluxes <- function(.) {

  # save current state for balance check  
  .super$state_pars$previous_state <- .super$state$cpools  

  # count number of loss fluxes, starts with 1 for respiration
  n_loss_fluxes <- 1
  if(!is.na(.super$pars$leaching_outflux)) 
    n_loss_fluxes <- n_loss_fluxes + 1
  lfm <- matrix(0, nrow=n_loss_fluxes )
  .super$state_pars$n_loss_fluxes <- n_loss_fluxes 
  .super$state_pars$lfm           <- lfm
 
  # calculate state parameters
  .$calc_state_pars()

  # run solver
  parmslist <- list(steadystate=F)
  .super$state_pars$solver_out <- 
    tryCatch(
      .$solver(rbind(.super$state$cpools,lfm), c(0L,1L), .$solver_func, parms=parmslist ),
      error = function(c) {
        print(c)
        matrix(ncol=.super$pars$n_pools+1)
      })
  
  # assign solved values to state
  .super$state$cpools[,1] <- .super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],2:(.super$pars$n_pools+1)]
  
  # assign losses
  .super$state$leaching    <- 0
  .super$state$respiration <- 
    as.numeric(.super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],.super$pars$n_pools+2]) 
  # APW: need a more robust way to id which pools are the other loss pools
  if(!is.na(.super$pars$leaching_outflux)) 
    .super$state$leaching  <- 
      as.numeric(.super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],.super$pars$n_pools+3]) 
   
  # mass balance check
  .$mass_balance(steadystate=F)
}



# steady-state solver for a system of n pools
f_steadystate_npools <- function(.) {

  # calculate state parameters
  .$calc_state_pars()

  # run solver
  parmslist <- list(steadystate=T)
  .super$state_pars$solver_steadystate_out <- 
    tryCatch(
      .$solver_steadystate(.super$state$cpools, 0L, func=.$solver_func, parms=parmslist ),
      error = function(c) {
        print(c)
        list(y=rep(NA,.super$pars$n_pools))
      })

  # assign solved values to state
  .super$state$cpools[,1] <- .super$state_pars$solver_steadystate_out$y
 
  # calculate respiration and other loss fluxes (e.g. leaching)
  .$calc_loss_fluxes()

  # mass balance check 
  .$mass_balance(steadystate=T)
}



# NULL solver to prevent steady-state initialisation, i.e. to initialise with values specified in input
f_steadystate_null   <- function(.) NULL 



### END ###
