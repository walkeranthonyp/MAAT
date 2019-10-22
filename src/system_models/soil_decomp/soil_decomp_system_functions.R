################################
#
# MAAT soil_decomp system representation functions (SRFs) 
# 
# AWalker, Matt Craig, October 2019 
#
################################



################################
# soil_decomp system function one
# - 2 pool microbial model

f_sys_m2pool <- function(.) {

  .super$state$solver_out  <- .$plsoda(.super$state$c_pools, c(0,1), .$solver_func )
  .super$state$c_pools[,1] <- .super$state$solver_out[dim(.super$state$solver_out)[1],2:3]
}



### END ###
