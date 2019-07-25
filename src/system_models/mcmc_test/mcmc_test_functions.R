################################
#
# MCMC_test process representation functions (PRFs)
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



### END ###
