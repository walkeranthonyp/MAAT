################################
#
# gwater_rt for MAAT object functions
# from Dai et al 2017 WRR 
#
# AWalker November 2017
#
################################

library(proto)
source('gwater_rt_functions.R')
source('gwater_rt_system_functions.R')



# GWATER_RT OBJECT
###############################################################################

# use generic template
setwd('..')
source('generic_model_object.R')
gwater_rt_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('gwater_rt')



# assign object functions
###########################################################################
gwater_rt_object$name              <- 'gwater_rt'
gwater_rt_object$build             <- build
gwater_rt_object$configure         <- configure
gwater_rt_object$configure_sublist <- configure_sublist  
gwater_rt_object$configure_unique  <- NULL  
gwater_rt_object$configure_test    <- configure_test
gwater_rt_object$run_met           <- run_met  



# use default run function
###########################################################################



# assign object variables 
###########################################################################

# function names
####################################
gwater_rt_object$fnames <- list(
  sys    = 'f_sys_daiye',   # system model
  rechrg = 'f_rechrg_lin',  # recharge function
  geol   = 'f_geol_single'  # geology/transport function
)


# environment
####################################
gwater_rt_object$env <- list(
  h1     = 180,               # hydraulic head left  (x=0) boundary condition (m)
  h2     = 100,               # hydraulic head right (x=L) boundary condition (m)
  precip = 1524               # precipitation                                 (mm)
)


# state
####################################
gwater_rt_object$state <- list(
  recharge = numeric(1),      # recharge rate                    (m d-1?)
  h        = numeric(21)      # hydraulic head across the domain (Pa)
)


# state parameters (i.e. calculated parameters)
####################################
gwater_rt_object$state_pars <- list(
  x      = numeric(0)
)


# parameters
####################################
gwater_rt_object$pars   <- list(
  nx  = 21,                   # horizontal discretisation of geological domain   
  K   = 15,                   # hydraulic conductivity for single domain geology 
  K1  = 20,                   # hydraulic conductivity for double domain geology           
  K2  = 10,                   # hydraulic conductivity for double domain geology         
  a   = 3.35,                 # recharge scalar power law (unitless)
  b   = 0.15,                 # recharge scalar linear    (unitless)
  L   = 1e4,                  # horizontal domain length (m)
  out_hsub = 13               # subscript of .$state$h vector for 'slim' output 
)


# run control parameters
####################################
gwater_rt_object$cpars <- list(
  verbose       = F,          # write diagnostic output during runtime 
  cverbose      = F,          # write configuration output during runtime 
  output        = 'run'       # type of output from run function
)



# output functions
#######################################################################        

gwater_rt_object$output <- function(.) {
  if(.$cpars$output=='run') {
    lout <- .$state$h
  } else if(.$cpars$output=='slim') {
    lout <- c(recharge=.$state$recharge,h=.$state$h[.$pars$out_hsub])
  } else stop()
  
  lout
}    



# test functions
#######################################################################        

gwater_rt_object$.test <- function(., verbose=F, 
                                   gwater_rt.precip=1524, gwater_rt.rechrg='f_rechrg_lin', gwater_rt.geol='f_geol_single' ) {
  
  if(verbose) {
    str(.)
    print(.$env)
  }
  
  .$build()

  .$cpars$verbose <- verbose
  .$cpars$output  <-'run'
  
  .$env$precip    <- gwater_rt.precip
  .$fnames$rechrg <- gwater_rt.rechrg
  .$fnames$geol   <- gwater_rt.geol

  .$configure_test()  
  .$run()
}



### END ###
