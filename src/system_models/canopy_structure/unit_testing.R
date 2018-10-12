###########################
#
# Unit testing of canopy_structure structure  
#
# AWalker October 2018
#
###########################



### Load model scripts 
###############################



# Canopy structure MAAT
rm(list=ls())
library(lattice)
source('canopy_structure_object.R')

# canopy_structure
# canopy_structure_object$fnames$can_scale_light <- 'f_canlight_beerslaw_goudriaan'
canopy_structure_object$.test()
canopy_structure_object$state$integrated$A

canopy_structure_object$pars$layers <- 40
canopy_structure_object$state$vert
canopy_structure_object$state_pars


# canopy_structure ACa curve
source('canopy_structure_object.R')
canopy_structure_object$.test_aca()
canopy_structure_object$env
canopy_structure_object$state



### END ###
