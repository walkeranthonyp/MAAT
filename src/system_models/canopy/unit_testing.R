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
canopy_object$fnames$can_scale_light <- 'f_canlight_beerslaw_goudriaan'
canopy_object$.test()
canopy_object$state$integrated$A

canopy_object$pars$layers <- 40
canopy_object$state$vert
canopy_object$state_pars

canopy_object$fnames$can_scale_light <- 'f_canlight_goudriaan'
canopy_object$.test()
canopy_object$state$integrated$A
canopy_object$pars
canopy_object$fnames
canopy_object$env

canopy_object$leaf$cpars
canopy_object$leaf$output
canopy_object$leaf$state
canopy_object$leaf$state_pars
canopy_object$leaf$fnames
canopy_object$leaf$env


# canopy ACa curve
source('canopy_object.R')
canopy_object$.test_aca()
canopy_object$env
canopy_object$state

canopy_object$fnames$rt <- 'f_rt_goudriaan'
canopy_object$.test_aca()
canopy_object$.test_aca(rs='f_rs_medlyn2011')



### END ###