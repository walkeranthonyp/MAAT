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
  .super$state$cpools[,1]      <- .super$state_pars$solver_out[dim(.super$state_pars$solver_out)[1],2:(.super$state$cpools_n+1)]
  #print(.super$state$cpools[,1])
}



### END ###
