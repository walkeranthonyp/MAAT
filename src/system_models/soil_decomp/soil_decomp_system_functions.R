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

  .super$state_pars$solver_out <- tryCatch(
    .$solver(.super$state$cpools, c(0,1), .$solver_func ), 
    error = function(c) {
      print(c)
      matrix(ncol=.super$pars$n_pools+1)
    })
  .super$state$cpools[,1] <- .super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],2:(.super$pars$n_pools+1)]
}


# steady-state solver for a system of n pools
f_steadystate_npools <- function(.) {

  .super$state_pars$solver_steadystate_out <- tryCatch(
    .$solver_steadystate(.super$state$cpools, 0, func=.$solver_func ),
    error = function(c) {
      print(c)
      list(y=rep(NA,.super$pars$n_pools))
    })
  .super$state$cpools[,1] <- .super$state_pars$solver_steadystate_out$y
}


# NULL solver to prevent steady-state initialisation, i.e. to initialise with values specified in input
f_steadystate_null   <- function(.) NULL 



### END ###
