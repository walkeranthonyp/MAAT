################################
#
# gwater_rt for MAAT system functions
# from Dai et al 2017 WRR 
#
# AWalker February 2018
#
################################



################################
# system model for simple ground water and reactive transport model
# Dai et al 2017
f_gwatersys_daiye <- function(.) {

  # calculate state parameters
  .$state_pars$x <-  seq(0, .$pars$L, (.$pars$L-0)/(.$pars$nx-1) )
  
  # Dai & Ye model 
  # recharge
  .$state$recharge <- get(.$fnames$rechrg)(.)
  #print(get(.$fnames$recharge)(.))
  
  # geological transport
  .$state$h        <- get(.$fnames$geol)(.)

}      
