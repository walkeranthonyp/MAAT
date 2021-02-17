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
canopy_object$fns$k.vcmax
canopy_object$fns$k.vcmax(var='vcmax')
canopy_object$fns$scale.ca_conc
canopy_object$fns$scale.ca_conc(1, vlist='env', var='ca_conc' )
canopy_object$fns$scale.vcmax(1,var='vcmax')
canopy_object$fns$layer0.vcmax(var='vcmax')
canopy_object$fns$scale.vcmax(1:10,var='vcmax')

canopy_object$.test(canopy.rt='f_rt_goudriaan')
canopy_object$state_pars

source('canopy_object.R') 
canopy_object$.test(
  verbose=F,
  canopy.lai=5,
  canopy.layers=20,
  canopy.rt='f_rt_goudriaan',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$state$integrated$apar
canopy_object$state_pars
canopy_object$state$vert$layer$apar
canopy_object$state$vert$sun$apar
canopy_object$state$vert$shade$apar

source('canopy_object.R') 
canopy_object$.test(
  verbose=F,
  canopy.par=1000,
  canopy.lai=6,
  canopy.layers=10,
  canopy.rt='f_rt_norman'
)
canopy_object$state$integrated$apar
canopy_object$state_pars
canopy_object$state$vert$layer$apar
canopy_object$state$vert$sun$apar
canopy_object$state$vert$shade$apar
canopy_object$env
  
canopy_object$.test(canopy.par=1000, canopy.ca_conc=200)
canopy_object$state$integrated$A

canopy_object$pars
canopy_object$fnames
canopy_object$env

# canopy ACa curve
source('canopy_object.R')
canopy_object$.test_aca()
canopy_object$.test_aca(canopy.rt='f_rt_goudriaan')
canopy_object$.test_aca(canopy.rt='f_rt_norman')
canopy_object$.test_aca(canopy.rt='f_rt_norman',
                        canopy.par=2000, canopy.ca_conc=c(200))
canopy_object$.test_aca(leaf.rs='f_rs_medlyn2011')
canopy_object$.test_aca(canopy.par=1320)
canopy_object$fnames
canopy_object$env
canopy_object$pars
canopy_object$state
canopy_object$state$vert
canopy_object$leaf
canopy_object$leaf$fns
canopy_object$leaf$fns$Ajg



### END ###