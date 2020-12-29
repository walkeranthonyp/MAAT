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
canopy_object$state$vert
canopy_object$state_pars

canopy_object$.test(canopy.par=1000, canopy.ca_conc=200)
canopy_object$state$integrated$A
canopy_object$pars
canopy_object$fnames
canopy_object$env


# canopy ACa curve
source('canopy_object.R')
canopy_object$.test_aca()
canopy_object$.test_aca(canopy.rt='f_rt_goudriaan')
canopy_object$.test_aca(leaf.rs='f_rs_medlyn2011')
canopy_object$.test_aca(canopy.par=1320)
canopy_object$env
canopy_object$pars
canopy_object$state
canopy_object$state$vert
canopy_object$leaf
canopy_object$leaf$fns
canopy_object$leaf$fns$Ajg



### END ###