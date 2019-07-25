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
  .$env$linreg_x     <- .super$dataf$met[,'mcmc_test.linreg_x']

  # generate uncertainty on target parameters 
  .$pars$a           <- rnorm(length(.$env$linreg_x), .$pars$syn_a_mu, .$pars$syn_a_sd )
  .$pars$b           <- rnorm(length(.$env$linreg_x), .$pars$syn_b_mu, .$pars$syn_b_sd )

  # generate synthetic data
  .super$dataf$obs   <- .$fns$reg_func()
  .super$dataf$obsse <- abs(rnorm(length(.$env$linreg_x), mean=0, sd=.$pars$obs_error ))

  # reallocate variables
  .$pars$a           <- numeric(1) 
  .$pars$b           <- numeric(1) 
  .$env$linreg_x     <- numeric(1) 
}



### END ###
