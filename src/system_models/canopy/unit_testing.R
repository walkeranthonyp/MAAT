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



### Run test functions
###############################

# canopy default
source('canopy_object.R') 
# model set up
canopy_object$fnames
# run test (verbose)
canopy_object$.test()
# run test (non-verbose)
canopy_object$.test(verbose=F)
# check model data structure
canopy_object$state$integrated$A
canopy_object$state$integrated
canopy_object$state$vert
canopy_object$state_pars


# 2 bigleaf
source('canopy_object.R') 
canopy_object$.test_2bigleaf()
canopy_object$.test_2bigleaf(canopy.can_clump=0.5)
unlist(canopy_object$state$integrated)
unlist(canopy_object$env)
canopy_object$state$vert
canopy_object$state_pars


# test alternate RT functions
source('canopy_object.R') 
canopy_object$.test(canopy.par=1000, canopy.ca_conc=200)
canopy_object$.test(canopy.rt='f_rt_goudriaan')
canopy_object$.test(
  verbose=T,
  canopy.lai=6,
  canopy.layers=20,
  canopy.rt='f_rt_goudriaan',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=T,
  canopy.lai=6,
  canopy.layers=20,
  canopy.rt='f_rt_norman',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.layers=20,
  canopy.rt='f_rt_norman',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=T,
  canopy.lai=6,
  canopy.layers=200,
  canopy.rt='f_rt_norman',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$state$integrated$apar
canopy_object$state_pars
canopy_object$state$vert$layer$apar
canopy_object$state$vert$sun$apar
canopy_object$state$vert$shade$apar

canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.layers=14,
  canopy.rt='f_rt_norman',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0, canopy.can_disc.layer1=6/14,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.layers=20,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0.1, canopy.can_disc.layer1=0.02,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0.5, canopy.can_disc.layer1=0.02,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0.5, canopy.can_disc.layer1=0.5,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0.7, canopy.can_disc.layer1=0.25,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=1.0, canopy.can_disc.layer1=0.2,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.layers=4,
  canopy.rt='f_rt_norman',
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0.05, canopy.can_disc.layer1=0.05,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)
canopy_object$.test(
  verbose=F,
  canopy.lai=6,
  canopy.rt='f_rt_norman',
  canopy.canopy_discretisation='f_canopy_discretisation_fixedbins_exp',
  canopy.k.can_disc=0.0, canopy.can_disc.layer1=6/1000,
  canopy.diffalbedo='f_diffalbedo_goudriaan'
)


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
  

# canopy ACa curves
source('canopy_object.R')
canopy_object$.test_aca()
canopy_object$.test_aca(canopy.rt='f_rt_goudriaan')
canopy_object$.test_aca(canopy.rt='f_rt_norman')
canopy_object$.test_aca(leaf.rs='f_rs_medlyn2011')
canopy_object$.test_aca(canopy.par=1320)
canopy_object$fnames
canopy_object$env
canopy_object$pars
canopy_object$state
canopy_object$state$vert


# check embedded leaf object
canopy_object$leaf
canopy_object$leaf$fns
canopy_object$leaf$fns$Ajg


# check canopy trait/environment scaling functions
canopy_object$fns$k.vcmax
canopy_object$fns$k.vcmax(var='vcmax')
canopy_object$fns$scale.ca_conc
canopy_object$fns$scale.ca_conc(1, vlist='env', var='ca_conc' )
canopy_object$fns$scale.vcmax(1,var='vcmax')
canopy_object$fns$layer0.vcmax(var='vcmax')
canopy_object$fns$scale.vcmax(1:10,var='vcmax')



### END ###