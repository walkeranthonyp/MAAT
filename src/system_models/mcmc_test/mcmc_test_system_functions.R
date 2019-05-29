################################
#
# mcmc_test system representation functions (SRFs)
#
# AWalker, Dan Lu, Abbey Johnson
# July 2018
#
################################



################################
# MCMC test system function

f_mcmc_testsys_mixture <- function(.) {

  # Mixture model that combines three normal distributions centered at mu1, mu2, mu3
  # the function evaluates the probability of the proposal vector based on the multi-modal distribution
  x_proposal <- c(.$pars$proposal1, .$pars$proposal2, .$pars$proposal3, .$pars$proposal4 )

  print('proposal = ')
  print(x_proposal)
  print("")

  # calculate proposal probability for each of the three distributions
  p1 <- .$pars$mixture_scale * .$pars$height1 * dnorm(x_proposal, .$pars$mu1, .$pars$sd1 )
  p2 <- .$pars$mixture_scale * .$pars$height2 * dnorm(x_proposal, .$pars$mu2, .$pars$sd2 )
  p3 <- .$pars$mixture_scale * .$pars$height3 * dnorm(x_proposal, .$pars$mu3, .$pars$sd3 )

  # return combined probability
  .$state$mixture_p[] <- sum(prod(p1) + prod(p2) + prod(p3))

  print(paste0('model evalution = ', .$state$mixture_p[]))
  #print(.$state$mixture_p[])
  print('')

  return(.$state$mixture_p[])
}


f_mcmc_testsys_regression <- function(.) {

  #print(c(.$pars$a, .$pars$b, .$env$linreg_x ))

  # line function used as the model in a linear regression
  .$state$linreg_y[] <- get(.$fnames$reg_func)(.)
}



### END ###
