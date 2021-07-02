################################
#
# MAAT soil_decomp system representation functions (SRFs) 
# 
# AWalker, Matt Craig, October 2019 
#
################################



################################
# soil_decomp system function 

f_sys_npools <- function(.) {
  .super$state_pars$solver_out <- .$solver(.super$state$cpools, c(0,1), .$solver_func )
  .super$state$cpools[,1]      <- .super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],2:(.super$pars$n_pools+1)]
}

f_steadystate_npools <- function(.) {
  .super$state_pars$solver_steadystate_out <- .$solver_steadystate(.super$state$cpools, 0, .$solver_func )
  # APW Matt: this is just copied form above so assumes the output from stode is in the same format as lsoda
  #           obvs that could be wrong so this might need editing
  .super$state$cpools[,1]                  <- 
    .super$state_pars$solver_steadystate_out[dim(.super$state_pars$solver_steadystate_out)[1],2:(.super$pars$n_pools+1)]
}

f_steadystate_null   <- function(.) NULL 


### END ###
