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
canopy_object$.test()
canopy_object$state$integrated$A

canopy_object$pars$layers <- 40

canopy_object$pars
canopy_object$fnames

canopy_object$leaf$cpars
canopy_object$leaf$output
canopy_object$leaf$state
canopy_object$leaf$state_pars
canopy_object$leaf$fnames


# canopy ACa curve
canopy_object$.test_aca()



### END ###
