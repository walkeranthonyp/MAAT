###########################
#
# allometry for MAAT object unit testing
# 
# AWalker March 2017
#
###########################



### Load model scripts 
###############################

# allometry
source('allometry_object.R')
allometry_object$.test()
allometry_object$.test(dbh=1)
allometry_object$.test()

allometry_object$.test(verbose=T)

allometry_object$fnames
allometry_object$fns
allometry_object$fns$crown_lai
allometry_object$state
allometry_object$state_pars
allometry_object$pars
allometry_object$env



### END ###