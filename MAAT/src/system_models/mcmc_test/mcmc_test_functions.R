################################
#
# MAAT MCMC_test process representation functions (PRFs)
#
# AWalker July 2018
#
################################


### FUNCTIONS
################################

# regression functions
f_reg_func_linear <- function(.) {
  .super$pars$a + .super$pars$b * .super$env$linreg_x
}


# functions to generate synthetic data
f_reg_func_linear_gen_obs <- function(.) {

  # assign met data
  print('',quote=F)
  print('Gen obs:',quote=F)
  #print(.$dataf$met)
  #print(.super$dataf$met)
  .$env$linreg_x     <- .super$dataf$met['mcmc_test.linreg_x',]

  # generate uncertainty on target parameters
  .$pars$a           <- rnorm(length(.$env$linreg_x), .$pars$syn_a_mu, .$pars$syn_a_sd )
  .$pars$b           <- rnorm(length(.$env$linreg_x), .$pars$syn_b_mu, .$pars$syn_b_sd )
  print('parameter a:',quote=F)
  print(.$pars$syn_a_mu)
  print('with error:',quote=F)
  print(.$pars$a)
  print('',quote=F)
  print('parameter b:',quote=F)
  print(.$pars$syn_b_mu)
  print('with error:',quote=F)
  print(.$pars$b)

  # generate synthetic data
  .super$dataf$obs   <- .$fns$reg_func()
  .super$dataf$obsse <- abs(rnorm(length(.$env$linreg_x), mean=0, sd=.$pars$obs_error ))
  print('',quote=F)
  print('x vals:',quote=F)
  print(.$env$linreg_x,quote=F)
  print('',quote=F)
  print('y vals:',quote=F)
  print(.super$dataf$obs,quote=F)
  print('',quote=F)

  # reallocate variables
  .$pars$a           <- numeric(1)
  .$pars$b           <- numeric(1)
  .$env$linreg_x     <- numeric(1)
}



### END ###
