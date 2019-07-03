################################
#
# gwater_rt for MAAT system functions
# from Dai et al 2017 WRR 
#
# AWalker February 2018
#
################################



################################
# system model for simple ground water model (does not include recative transport component in Dai et al 2017)

f_sys_daiye <- function(.) {

  # calculate state parameters
  .super$state_pars$x   <-  seq(0, .super$pars$L, (.super$pars$L-0)/(.super$pars$nx-1) )
  
  # recharge
  .super$state$recharge <- .$rechrg()
  
  # geological transport
  .super$state$h        <- .$geol()

}     



### END ### 
