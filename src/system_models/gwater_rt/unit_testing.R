###########################
#
# gwater_rt for MAAT object unit testing
# 
# AWalker March 2017
#
###########################



### Load model scripts 
###############################

# gwater_rt
source('gwater_rt_object.R')
gwater_rt_object$.test(verbose=F)
gwater_rt_object$fnames
gwater_rt_object$state
gwater_rt_object$state_pars
gwater_rt_object$pars
gwater_rt_object$env

gwater_rt_object$.test(gwater_rt.rechrg='f_rechrg_power')
gwater_rt_object$fnames
gwater_rt_object$state

gwater_rt_object$.test(gwater_rt.geol='f_geol_double')
gwater_rt_object$fnames
gwater_rt_object$state

gwater_rt_object$.test(gwater_rt.rechrg='f_rechrg_power', gwater_rt.geol='f_geol_double')
gwater_rt_object$fnames
gwater_rt_object$state



### END ###
