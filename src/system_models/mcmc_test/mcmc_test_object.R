################################
#
# MAAT mcmc_test model object
#
# AWalker November 2017
#
################################

library(proto)
source('mcmc_test_functions.R')
source('mcmc_test_system_functions.R')



# MCMC_TEST OBJECT
###############################################################################

# use generic mcmc_test
setwd('..')
source('generic_model_object.R')
mcmc_test_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('mcmc_test')



# assign object functions
###########################################################################
mcmc_test_object$name <- 'mcmc_test'


mcmc_test_object$configure_unique <- function(., init=F, flist=NULL ) {

  if(any(names(flist)=='reg_func')) {

   .$reg_func_gen_obs <- get(paste0(.$fnames$reg_func,'_gen_obs'), pos=1 )
   .$reg_func_gen_obs()

  }
}


# assign object variables
###########################################################################

# function names
####################################
mcmc_test_object$fnames <- list(
  sys      = 'f_sys_mixture',
  reg_func = 'f_reg_func_linear'
)


# parameters
####################################
mcmc_test_object$pars   <- list(

  # mixture model parameters
  mixture_scale = 1e12,
  height1       = 0.1,
  height2       = 0.3,
  height3       = 0.6,
  mu1           = -8,
  mu2           = 0,
  mu3           = 8,
  sd1           = 1,
  sd2           = 1,
  sd3           = 1,
  proposal1     = 1,
  proposal2     = 1,
  proposal3     = 1,
  proposal4     = 1,

  # linear regresssion model parameters
  obs_error     = 1,
  syn_a_mu      = 2,
  syn_b_mu      = 7,
  syn_a_sd      = 3,
  syn_b_sd      = 2,
  a             = 7,
  b             = 2
)


# environment
####################################
mcmc_test_object$env <- list(
  dummy    = 1,
  linreg_x = 10
)


# state
####################################
mcmc_test_object$state <- list(
  mixture_p = numeric(1), # probability output from mixture model
  linreg_y  = numeric(1)  # y value output from line function
)


# state parameters (i.e. calculated parameters)
####################################
mcmc_test_object$state_pars <- list(
  none      = numeric(0)
)


# run control parameters
####################################
mcmc_test_object$cpars <- list(
  verbose      = F,          # write diagnostic output during runtime
  cverbose     = F,          # write diagnostic output from configure function
  output       = 'mixture'   # type of output from run function
)



# output functions
#######################################################################

f_output_mcmc_test_mixture <- function(.) {
  .$state$mixture_p
}

f_output_mcmc_test_regression <- function(.) {
  .$state$linreg_y
}



# test functions
#######################################################################

mcmc_test_object$.test_mixture <- function(., verbose=F, cverbose=F ) {

  if(verbose) str(.)
  .$cpars$maat_gen_obs <- F
  .$build(mod_out='mixture', switches=c(F,verbose,cverbose) )

  .$fnames$sys <- 'f_sys_mixture'
  .$configure_test()
  .$run()
}

mcmc_test_object$.test_linreg <- function(., verbose=F, cverbose=F ) {

  if(verbose) str(.)
  .$cpars$maat_gen_obs <- F
  .$build(mod_out='regression', switches=c(F,verbose,cverbose) )

  .$fnames$sys <- 'f_sys_regression'
  .$configure_test()
  .$run()
}



### END ###
