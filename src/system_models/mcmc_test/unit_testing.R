###########################
#
# mcmc_test unit testing
# 
# AWalker July 2018 
#
###########################



### Load model scripts 
###############################

# mcmc_test
source('mcmc_test_object.R')
mcmc_test_object$.test_mixture(verbose=F,verbose_loop=F)

mcmc_test_object$.test_linreg(verbose=F,verbose_loop=F)

mcmc_test_object$fnames
mcmc_test_object$state
mcmc_test_object$state_pars
mcmc_test_object$pars
mcmc_test_object$env



### END ###
