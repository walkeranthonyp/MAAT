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
canopy_structure_object$.test()
canopy_structure_object$state$integrated$A
canopy_structure_object$state$vert$layer
canopy_structure_object$state
canopy_structure_object$state_pars
canopy_structure_object$pars

canopy_structure_object$pars$layers <- 40



# canopy_structure ACa curve
source('canopy_structure_object.R')
canopy_structure_object$.test_aca()
canopy_structure_object$env
canopy_structure_object$state



### END ###
