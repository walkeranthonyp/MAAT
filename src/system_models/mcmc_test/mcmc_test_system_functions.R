################################
#
# MAAT mcmc_test system representation functions (SRFs)
#
# AWalker, Dan Lu, Abbey Johnson
# July 2018
#
################################



################################
# MCMC test system function

f_sys_mixture <- function(.) {

  # Mixture model that combines three normal distributions centered at mu1, mu2, mu3
  # the function evaluates the probability of the proposal vector based on the multi-modal distribution
  x_proposal <- c(.super$pars$proposal1, .super$pars$proposal2, .super$pars$proposal3, .super$pars$proposal4 )

  # print('proposal = ')
  # print(x_proposal)

  # calculate proposal probability for each of the three distributions
  p1 <- .super$pars$mixture_scale * dnorm(x_proposal, .super$pars$mu1, .super$pars$sd1 )
  p2 <- .super$pars$mixture_scale * dnorm(x_proposal, .super$pars$mu2, .super$pars$sd2 )
  p3 <- .super$pars$mixture_scale * dnorm(x_proposal, .super$pars$mu3, .super$pars$sd3 )

  # return combined probability
  .super$state$mixture_p[] <- sum(.super$pars$height1 * prod(p1) + .super$pars$height2 * prod(p2) + .super$pars$height3 * prod(p3))

  # print(paste0('model evalution = ', .super$state$mixture_p))

  return(.$state$mixture_p)
}


f_sys_regression <- function(.) {

  #print(c(.$pars$a, .$pars$b, .$env$linreg_x ))

  # line function used as the model in a linear regression
  .super$state$linreg_y[] <- .$reg_func()
}



### END ###
