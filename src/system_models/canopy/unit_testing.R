###########################
#
# Unit testing of the MAAT models and wrapper
#
# AWalker
# May 2015
#
###########################



### Load model scripts 
###############################



# Canopy MAAT
rm(list=ls())
library(lattice)
source('canopy_object.R')

# canopy
canopy_object$fnames$rt <- 'f_rt_beerslaw_goudriaan'
canopy_object$.test()
canopy_object$state$integrated$A
canopy_object$state$vert
canopy_object$state_pars

canopy_object$fnames$can_scale_light <- 'f_rt_goudriaan'
canopy_object$.test()
canopy_object$state$integrated$A
canopy_object$pars
canopy_object$fnames
canopy_object$env

source('canopy_object.R')
canopy_object$fnames$cansys <- 'f_cansys_multilayer'
canopy_object$fnames$rt     <- 'f_rt_goudriaan'
canopy_object$.test()
canopy_object$env

canopy_object$leaf$cpars
canopy_object$leaf$state
canopy_object$leaf$state_pars
canopy_object$leaf$fnames
canopy_object$leaf$env


# canopy ACa curve
source('canopy_object.R')
canopy_object$fnames$cansys <- 'f_cansys_multilayer'
canopy_object$fnames$rt     <- 'f_rt_goudriaan'
canopy_object$.test()
canopy_object$.test(par=1320, ca_conc=400 )
canopy_object$.test_aca()
canopy_object$pars$k_layer  <- 0
canopy_object$.test_aca()
canopy_object$pars$k_layer  <- 0.5
canopy_object$pars$layers   <- 40
canopy_object$.test_aca()
canopy_object$.test_aca(rs='f_rs_medlyn2011')
canopy_object$fnames$k_vcmax <- 'f_k_vcmax_lloyd2012'
canopy_object$.test(par=1320, ca_conc=400 )
canopy_object$.test_aca()
canopy_object$env$par <- 2000
canopy_object$.test_aca()
canopy_object$env
canopy_object$pars
canopy_object$state
canopy_object$state$vert
canopy_object$fnames$rt     <- 'f_rt_beerslaw_goudriaan'
canopy_object$.test_aca()


canopy_object$fnames$cansys <- 'f_cansys_bigleaf_s1992'
canopy_object$.test_aca()



### END ###