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
gwater_rt_object$.test_gwater_rt(verbose=F,verbose_loop=F)
gwater_rt_object$fnames
gwater_rt_object$state
gwater_rt_object$state_pars
gwater_rt_object$pars
gwater_rt_object$env









