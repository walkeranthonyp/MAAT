################################
#
# Wrapper functions for MCMC runs
#
# AJohnson, AWalker July 2019
#
################################



# DEMC functions
################################

init_mcmc_demc <- function(.) NULL

# generate proposal using DE-MC algorithm
proposal_generate_mcmc_demc <- function(., j ) {

  # scaling factor
  d          <- dim(.$dataf$pars)[1]
  gamma_star <- 2.38 / sqrt(d + d)

  # b-value should be small compared to width of target distribution; specifies range for drawn "randomization" value
  b_rand  <- 0.01

  uniform_r <- runif(1, min = (-b_rand), max = b_rand)

  # evaluate for each chain
  for (ii in 1:.$dataf$lp) {

    # randomly select two different numbers R1 and R2 unequal to j, from uniform distribution without replacement
    R1 <- 0
    R2 <- 0
    while ((R1 == 0) | (R1 == ii))               R1 <- ceiling(runif(1, min = 0, max = 1) * .$dataf$lp)
    while ((R2 == 0) | (R2 == ii) | (R2 == R1))  R2 <- ceiling(runif(1, min = 0, max = 1) * .$dataf$lp)

    # evaluate for each parameter
    for (jj in 1:d) {

      # generate proposal via Differential Evolution
      .$dataf$pars[jj,ii] <- .$dataf$pars_array[jj,ii,j-1] + gamma_star * (.$dataf$pars_array[jj,R1,j-1] - .$dataf$pars_array[jj,R2,j-1]) + uniform_r

      # call boundary handling function
      .$boundary_handling(ii = ii, jj = jj)

    }
  }
}


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_mcmc_demc <- function(., j, lklihood ) {

  # Metropolis ratio
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[ ,j-1])

  alpha        <- pmin(1, metrop_ratio)

  # evaluate for each chain
  for(ii in 1:.$dataf$lp) {

    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
    accept <- log(alpha[ii]) > log(runif(1, min = 0, max = 1))

    .$dataf$pars_array[,ii,j]   <- if(accept) .$dataf$pars[,ii] else .$dataf$pars_array[,ii,j-1]
    .$dataf$pars_lklihood[ii,j] <- if(accept) lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]

    .$dataf$out_mcmc[ii, , j] <- if(accept | j == 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii, , (j - 1)]

  }
}



# DREAM MCMC functions
################################

# initialisation of DREAM algorithm
init_mcmc_dream <- function(.) {

  # number of parameters being estimated
  .$mcmc$d <- dim(.$dataf$pars)[1]

  # preallocate memory space for algorithmic variables
  .$mcmc$p_state       <- numeric(.$dataf$lp)
  .$mcmc$R             <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$dataf$lp - 1)
  .$mcmc$current_state <- matrix(data = 0, nrow = .$mcmc$d, ncol = .$dataf$lp)
  .$mcmc$draw          <- matrix(data = 0, nrow = .$dataf$lp - 1, ncol = .$dataf$lp)
  .$mcmc$lambda        <- matrix(data = 0, nrow = .$dataf$lp, ncol = 1)
  .$mcmc$jump          <- matrix(data = 0, nrow = .$mcmc$d, ncol = .$dataf$lp)

  # preallocate space for crossover variables
  .$mcmc$t         <- numeric(1)
  .$mcmc$d_star    <- numeric(1)
  .$mcmc$CR_burnin <- numeric(1)
  .$mcmc$sd_state  <- numeric(.$mcmc$d)
  .$mcmc$L         <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$del       <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$p_CR      <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$m         <- numeric(.$wpars$mcmc_chains)

  # if adapting selection of crossover probabilties
  if(.$wpars$mcmc_adapt_pCR) {

    # initialize crossover variables
    .$mcmc$t      <- 1
    .$mcmc$CR     <- 0
    .$mcmc$L[]    <- 0
    .$mcmc$d_star <- .$mcmc$d

    # initial probability of each crossover value
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

  # if not adapting p_CR
  } else {

    # initialize crossover probabilities
    .$mcmc$CR   <- numeric(.$wpars$mcmc_n_CR)
    .$mcmc$CR[] <- 1:.$wpars$mcmc_n_CR / .$wpars$mcmc_n_CR

    # initialize selection probability of crossover values
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

  }

  # burn-in period for adapting crossover selection probabilities
  .$mcmc$CR_burnin <- ceiling(.$wpars$mcmc_CR_burnin * .$wpars$mcmc_maxiter)

  # index of chains for Differential Evolution
  for (ii in 1:.$dataf$lp) .$mcmc$R[ii, ] <- setdiff(1:.$dataf$lp, ii)

}


# generate proposal using DREAM algorithm
proposal_generate_mcmc_dream <- function(., j ) {

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1])

  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[] <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)

  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[] <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)

  # if not adapting crossover values, compute standard deviation of each dimension/parameter
  if(!(.$wpars$mcmc_adapt_pCR)) {

    .$mcmc$sd_state[] <- apply(.$mcmc$current_state, 1, sd)

    # replace any 0's in standard deviation array with 1e-9 to avoid division by 0
    idx <- which(.$mcmc$sd_state == 0)
    .$mcmc$sd_state[idx] <- 1e-9
  }

  # create proposals
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta, 1, replace = T)

    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    if (.$wpars$mcmc_adapt_pCR) {
      # modify each dimension with probability CR each time a proposal vector is generated

      # generate crossover probability
      .$mcmc$m[ii] <- .$generate_CR()

      # calculate jump rate (scaling factor)
      gamma_d <- 2.38 / sqrt(2 * D * .$mcmc$d_star)

      # when gamma = 1, jump between different modes of the posterior (this is approx every 5 iterations with default p_gamma = 0.2)
      gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma))

      # compute jump differential evolution of ii-th chain
      .$mcmc$jump[1:.$mcmc$d, ii] <- .$wpars$mcmc_c_ergod * rnorm(.$mcmc$d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[1:.$mcmc$d, a] - .$mcmc$current_state[1:.$mcmc$d, b]), dim = 1)

      # compute proposal of ii-th chain
      .$dataf$pars[1:.$mcmc$d, ii] <- .$mcmc$current_state[1:.$mcmc$d, ii] + .$mcmc$jump[1:.$mcmc$d, ii]

      # replace each element (jj = 1,...,d) of the proposal with the corresponding current_state element
      #         using a binomial scheme with probability 1 - CR (CR = crossover probability)
      #         when CR = 1, all dimensions are updated jointly and d_star = d
      crossover <- logical(length = .$mcmc$d)
      for (jj in 1:.$mcmc$d) {
        if (runif(1, min = 0, max = 1) <= (1 - .$mcmc$CR)) {
          .$dataf$pars[jj, ii] <- .$mcmc$current_state[jj, ii]
        } else {
          crossover <- T
        }
      }

      # number of dimensions being updated
      .$mcmc$d_star <- length(crossover)

      # numerical check (in case no dimesnions are updated)
      if (.$mcmc$d_star == 0) .$mcmc$d_star <- 1

    } else {

      # select index of crossover value (weighted sample with replacement drawn from multinomial distribution)
      id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)

      # draw d values from uniform distribution between 0 and 1
      zz <- runif(.$mcmc$d)

      # derive subset A of selected dimensions
      A  <- which(zz < .$mcmc$CR[.$mcmc$id])

      #  how many dimensions are sampled (i.e., how many parameters will be updated jointly)
      .$mcmc$d_star <- length(A)

      # numerical check: make sure that A contains at least one value
      if (.$mcmc$d_star == 0) {
        A <- which.min(zz)
        .$mcmc$d_star <- 1
      }

      # calculate jump rate (scaling factor)
      gamma_d <- 2.38 / sqrt(2 * D * .$mcmc$d_star)

      # when gamma = 1, jump between different modes of the posterior (this is approx every 5 iterations with default p_gamma = 0.2)
      gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma))

      # compute jump differential evolution of ii-th chain
      .$mcmc$jump[A, ii] <- .$wpars$mcmc_c_ergod * rnorm(.$mcmc$d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[A, a] - .$mcmc$current_state[A, b]), dim = 1)

      # compute proposal of ii-th chain
      .$dataf$pars[1:.$mcmc$d, ii] <- .$mcmc$current_state[1:.$mcmc$d, ii] + .$mcmc$jump[1:.$mcmc$d, ii]

    }

    # call boundary handling function
    for (jj in 1:.$mcmc$d) .$mcmc_bdry_handling(j=j, ii = ii, jj = jj)

  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_mcmc_dream <- function(., j, lklihood) {

  # likelihood of current state
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ , j-1]

  # iterate through chains
  for (ii in 1:.$dataf$lp) {

    # Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[ii] - .$mcmc$p_state[ii]))

    if (alpha > runif(1, min = 0, max = 1)) {

      # accept proposal
      accept <- TRUE

      .$mcmc$current_state[1:.$mcmc$d, ii] <- .$dataf$pars[1:.$mcmc$d, ii]
      .$mcmc$p_state[ii]                   <- lklihood[ii]

      # append accepted proposal and probability density to storage arrays
      .$dataf$pars_array[1:.$mcmc$d, ii, j] <- .$mcmc$current_state[1:.$mcmc$d, ii]
      .$dataf$pars_lklihood[ii, j]          <- .$mcmc$p_state[ii]

    } else {

      # reject proposal
      accept <- FALSE

      # reset jump back to zero
      .$mcmc$jump[1:.$mcmc$d, ii] <- 0

      # repeat previous accepted proposal and probability density in storage arrays
      .$dataf$pars_array[1:.$mcmc$d, ii, j] <- .$dataf$pars_array[1:.$mcmc$d, ii, j-1]
      .$dataf$pars_lklihood[ii, j]          <- .$dataf$pars_lklihood[ii, j-1]

    }

    # compute squared normalized jumping distance
    if (.$wpars$mcmc_adapt_pCR & (.$mcmc$t < .$mcmc$CR_burnin)) .$calc_del(j = j, ii = ii)

    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
    # write output (model evaluations)
    .$dataf$out_mcmc[ii, , j] <- if(accept | j == 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii, , (j - 1)]

  }

  if (.$wpars$mcmc_adapt_pCR & (.$mcmc$t < .$mcmc$CR_burnin)) {

    # update the selection probability of crossover probabilities/values
    .$adapt_pCR()

    .$mcmc$t <- .$mcmc$t + 1

  } else if (.$wpars$mcmc_adapt_pCR & .$mcmc$t == .$mcmc$CR_burnin) {

    print('Adapted selection probabilities of crossover values = '); print(.$mcmc$p_CR)

    .$mcmc$t <- .$mcmc$t + 1

  } else if (.$wpars$mcmc_adapt_pCR & (.$mcmc$t > .$mcmc$CR_burnin)) {

    .$mcmc$t <- .$mcmc$t + 1

  }
}



# adaptive p_CR functions
################################

# function that generates/updates crossover values based on current probabilities
generate_CR <- function(.) {

  # sample m from numbers 1,...,n_CR using multinomial distribution (with probabilities p_CR)
  m <- sample(1:.$wpars$mcmc_n_CR, size = 1, replace = T, prob = .$mcmc$p_CR)

  # set crossover probability/value
  .$mcmc$CR <- m / .$wpars$mcmc_n_CR

  # index of which crosover probabilities/values are selected
  .$mcmc$L[m] <- .$mcmc$L[m] + 1

  return(m)
}


# function to compute the squared normalized jumping distance
calc_del <- function(., j, ii) {

  # compute standard deviation of each dimension/parameter
  .$mcmc$sd_state[] <- apply(.$dataf$pars_array, 1, sd)

  # replace any 0's in sd vector with 1e-9 to avoid division by 0
  idx                  <- which(.$mcmc$sd_state == 0)
  .$mcmc$sd_state[idx] <- 1e-9

  summation <- sum(((.$dataf$pars_array[1:.$mcmc$d, ii, j] - .$dataf$pars_array[1:.$mcmc$d, ii, j - 1]) / .$mcmc$sd_state)^2)

  .$mcmc$del[.$mcmc$m[ii]] <- .$mcmc$del[.$mcmc$m[ii]] + summation

}


# function that adapts crossover probabilities
adapt_pCR <- function(.) {

  # update the probability of the different crossover values being selected
  for (qq in 1:.$wpars$mcmc_n_CR) {
    # numerical check for divide by zero
    if (.$mcmc$L[qq] != 0) {
      .$mcmc$p_CR[qq] <- .$mcmc$t * .$wpars$mcmc_chains * (.$mcmc$del[qq] / .$mcmc$L[qq]) / sum(.$mcmc$del)
    }
  }

  # normalize
  .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)

}



# prior distribution functions
################################

# initialize chains with uniform distributions
mcmc_prior_uniform <- function(.) {

  # number of Markov chains
  n <- .$wpars$mcmc_chains

  # in case there are nested lists
  unlist(.$dynamic$pars)

  # number of parameters (dimensionality of parameter space)
  d <- length(unlist(.$dynamic$pars, recursive = T)) / 2
  print(names(unlist(.$dynamic$pars, recursive = T)))
  print(names(.$dynamic$pars))
  print(paste0('d = ', d))
  #print(names(.$dynamic$pars))
  print(.$dynamic$pars)

  # replace with periods

  # determine minimums and maximums for parameter values
  #min_vals <- rep(0, d); max_vals <- rep(0, d)
  #for(i in 1:d) {
  #  min_vals[i] <- .$dynamic$pars[[names(.$dynamic$pars)[i]]][['min']]
  #  max_vals[i] <- .$dynamic$pars[[names(.$dynamic$pars)[i]]][['max']]
  #}

  # draw priors from uniform distribution to create pars / proposal matrix
  #.$dataf$pars <- matrix(0, nrow = d, ncol = n)
  #for (jj in 1:d) {
  #  .$dataf$pars[jj, 1:n] <- runif(n, min = min_vals[jj], max = max_vals[jj])
    #nam <- paste0('dimension', jj, sep = '')
    #nam <- runif(n, min = min_vals[jj], max = max_vals[jj])
    #print(nam)
    #assign(nam, runif(n, min = min_vals[jj], max = max_vals[jj]))
  #}

  #attempt <- rbind(c(paste('dimension', 1:d)))
  #print(dimension1)
  #print('attempt')
  #print(attempt)
  #.$dataf$pars <- rbind(dimension1, dimension2, dimension3, dimension4)

  # assign parameter names
  #rownames(.$dataf$pars) <- names(.$dynamic$pars)
  #print(dim(.$dataf$pars))
  #print(.$dataf$pars[1, 1])

  # hard-coding for sake of urgency
  min_vals <- rep(0, 3); max_vals <- rep(0, 3)

  min_vals[1] <- .$dynamic$pars$leaf.atref$vcmax$min
  max_vals[1] <- .$dynamic$pars$leaf.atref$vcmax$max

  min_vals[2] <- .$dynamic$pars$leaf.atref$jmax$min
  max_vals[2] <- .$dynamic$pars$leaf.atref$jmax$max

  min_vals[3] <- .$dynamic$pars$leaf.theta_col_cj$min
  max_vals[3] <- .$dynamic$pars$leaf.theta_col_cj$max

  n <- .$wpars$mcmc_chains
  .$dataf$pars <- matrix(0, nrow = 3, ncol = n)
  for (jj in 1:3) {
    .$dataf$pars[jj, 1:n] <- runif(n, min = min_vals[jj], max = max_vals[jj])
  }

  rownames(.$dataf$pars) <- c('leaf.atref.vcmax', 'leaf.atref.jmax', 'leaf.theta_col_cj')

  # temporary
  .$mcmc$boundary_min <- min_vals
  .$mcmc$boundary_max <- max_vals

  #print('min vals = ')
  #print(min_vals)
  #print('max vals = ')
  #print(max_vals)
  #print('.$dataf$pars')
  #print(.$dataf$pars)
}


# initialize chains with multi-normal distribution (multivariate normal distribution)
mcmc_prior_normal <- function(.) {

  # add check to see if mean and sd are not specified

}


# option for initializing chains when restarting MCMC algorithm
mcmc_prior_none <- function(.) {
  # use the last accepted proposal from previous MCMC run
}


# boundary handling functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {

  #boundary_sample <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text = cs)))

  #min_vals <- rep(0, length(boundary_sample))
  #max_vals <- rep(0, length(boundary_sample))

  #for(jj in 1:length(boundary_sample)) {
  #  min_vals[jj] <- boundary_sample[[names(boundary_sample)[jj]]][['min']]
  #  max_vals[jj] <- boundary_sample[[names(boundary_sample)[jj]]][['max']]
  #}

  # ORIGINAL BOUNDARY HANDLING CODE
  # number of samples to be used in boundary handling
  #n <- 1e4
  #boundary_sample     <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text = cs)))
  #.$mcmc$boundary_min <- unlist(lapply(boundary_sample, min))
  #.$mcmc$boundary_max <- unlist(lapply(boundary_sample, max))
  #rm(boundary_sample)

  # again, temporarily hard-coding for sake of urgency
  #min_vals <- c(100.0, 70.0, 0.9)
  #max_vals <- c(200.0, 300.0, 1.0)

  #.$mcmc$boundary_min <- min_vals
  #.$mcmc$boundary_max <- max_vals

  print(".$mcmc$boundary_min")
  print(.$mcmc$boundary_min)
  print(".$mcmc$boundary_max")
  print(.$mcmc$boundary_max)

  #rm(boundary_sample)
}


# no boundary handling: should be chosen when search space is not theoretically restricted
mcmc_bdry_handling_none <- function(., j, ii, jj) {
  if ((j == .$wpars$mcmc_maxiter) & (ii == .$wpars$mcmc_chains) & (jj == .$mcmc$d)) {
    print('No option was chosen for MCMC boundary handling.')
  }
}


# restrict parameter proposals that are beyond the boundary (reset them back to the min/max bound)
mcmc_bdry_handling_bound <- function(., j, ii, jj) {

  # if outside bound of parameter space, restrict proposal value to corresponding dimension min
  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj]

  # if outside bound of parameter space, restrict proposal value to corresponding dimension max
  else if (.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) .$dataf$pars[jj, ii] <- .$mcmc$boundary_max[jj]
}


# reflect parameter proposals that are beyond the min/max values back accros the boundary
mcmc_bdry_handling_reflect <- function(., j, ii, jj) {

  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) {
    # if outside bound of parameter space, reflect proposal back across minimum boundary
    .$dataf$pars[jj, ii] <- 2 * .$mcmc$boundary_min[jj] - .$dataf$pars[jj, ii]
  } else if (.$dataf$pars[ii, jj] > .$mcmc$boundary_max[jj]) {
    # if outside bound of parameter space, reflect proposal back across maximum boundary
    .$dataf$pars[jj, ii] <- 2 * .$mcmc$boundary_max[jj] - .$dataf$pars[jj, ii]
  }

  # numerical check: see if new reflected proposal value is out of bounds
  if ((.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) | (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj])) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj] + runif(1, min = 0, max = 1) * (.$mcmc$boundary_max[jj] - .$mcmc$boundary_min[jj])
  }
}


# restrict parameter proposals that are beyond the boundary by "folding" them back across the boundary
# this boundary handling approach maintains detailed statistical balance, which is healthy for the MCMC algorithm
mcmc_bdry_handling_fold <- function(., j, ii, jj) {

  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_max[jj] - (.$mcmc$boundary_min[jj] - .$dataf$pars[jj, ii])
  } else if (.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj] + (.$dataf$pars[jj, ii] - .$mcmc$boundary_max[jj])
  }

  # numerical check: see if new "folded" proposal value is out of bounds
  if ((.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) | (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj])) {
    .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj] + runif(1, min = 0, max = 1) * (.$mcmc$boundary_max[jj] - .$mcmc$boundary_min[jj])
  }
}



# convergence diagnostic functions
#####################################

# option for not computing a convergence diagnostic; to be used during post-burn-in MCMC sampling
mcmc_converge_none <- function(., j) {
  if (j == .$wpars$mcmc_maxiter) print('No option was chosen to test for MCMC convergence.')
}


# subroutine calculating the R-statistic of Gelman and Rubin (convergence diagnostic)
mcmc_converge_Gelman_Rubin <- function(., j) {

  # IMPORTANT: beware of pseudo-convergence
  #            R_hat may sometimes be small in early iterations, looking as if it has converged,
  #            so  make sure to run algorithm for an appropriately long number of samples



  # number of parameters/dimensions being sampled
  d <- .$mcmc$d

  # iterate thgouth parameters with jj

  # number of Markov chains
  N <- .$wpars$mcmc_chains

  # iterate through chains with r

  # total number of samples in each chain
  t <- .$wpars$mcmc_maxiter

  # iterate through time samples with i

  # mcmc storage array of samples
  x <- .$dataf$pars_array

  # within-chain variance

  x_bar <- matrix(0, nrow = d, ncol = N)

  # REALLY NOT SURE THAT I'M DOING THESE SUMMATIONS CORRECTLY

  for (jj in 1:d) {
    for (r in 1:N) {
      # summation (probably better way to do this in the future, avoid loops)
      for (i in ceiling(t / 2):t) {
        x_bar[jj, r] <- x_bar[jj, r] + x[jj, r, i]
      }
    }
  }

  x_bar <- (2 / (t - 2)) * x_bar

  W <- rep(0, d)

  for (jj in 1:d) {
    summation <- 0
    # double summation (probably better way to do this in the future, avoid loops)
    for (r in 1:N) {
      for (i in ceiling(t / 2):t) {
        summation <- summation + (x[jj, r, i] - x_bar[jj, r]) * (x[jj, r, i] - x_bar[jj, r])
      }
    }
    W[jj] <- 2 / (N * (t - 2)) * summation
  }

  print('within-chain variance = '); print(W)

  # between chain variance

  B <- rep(0, d)

  x_double_bar <- rep(0, d)
  for (jj in 1:d) {
    # summation (probably better way to do this in the future, avoid loops)
    for (r in 1:N) {
      x_double_bar[jj] <- x_double_bar[jj] + x_bar[jj, r]
    }
  }

  x_double_bar <- (1 / N) * x_double_bar

  for (jj in 1:d) {
    # sumamtion (probably better way to do this in the future, avoid loops)
    summation <- 0
    for (r in 1:N) {
      summation <- summation + (x_bar[jj, r] - x_double_bar[jj]) * (x_bar[jj, r] - x_double_bar[jj])
    }
    B[jj] <- (t / (2 * (N - 1))) * summation
  }

  print('between-chain variance = '); print(B)

  # estimate variance of j-th paraemter of target distribution

  sigma_hat <- rep(0, d)

  for (jj in 1:d) {
    sigma_hat[jj] <- ((t - 2) / t) * W[jj] + (2 / t) * B[jj]
  }

  # R-statistic of Gelman and Rubin
  R_hat <- rep(0, d)

  for (jj in 1:d) {
    R_hat[jj] <- sqrt(((N + 1) / N) * (sigma_hat[jj] / W[jj]) - ((t - 2) / (N * t)))
  }

  # add R_hat to storage array
  counter <- j / .$wpars$mcmc_check_iter
  R_hat_new <- append(R_hat, j, after = 0)
  .$dataf$conv_check[counter, ] <- R_hat_new

  if (j == .$wpars$mcmc_maxiter) print(paste0("At iteration ", j, ", R-statistic of Gelman and Rubin = ")); print(R_hat)
  if (j == .$wpars$mcmc_maxiter) print(.$dataf$conv_check)
}



# likelihood functions
################################

# expects model output to be probability - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {

  # derive log density
  log(.$dataf$out)

  return(log(.$dataf$out))
}


# standard error probability density function with i.i.d. error residuals
f_proposal_lklihood_ssquared <- function(.) {

  # number of measured data points
  obs_n <- length(.$dataf$obs)

  # calculate error residual
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n / 2) * log(SSR)
}


# standard error probability density function with i.i.d. error residuals
# this function incorporates measurement errors (unlike "ssquared" option)
f_proposal_lklihood_ssquared_se <- function(.) {

  # read in measurement error; remove zeros from measurement error
  sspos <- which(.$dataf$obsse > 1e-9)

  # number of measured data points (that do not have zero uncertainty)
  obs_n <- length(sspos)

  # observed error (heteroscedastic and homoscedastic options)
  obsse <- if(.$wpars$mcmc_homosced) rep(mean(.$dataf$obsse[sspos]), obs_n)
           else .$dataf$obsse[sspos]

  # calculate error residual
  error_residual_matrix <- ( t(.$dataf$out)[sspos, ] - .$dataf$obs[sspos] ) / obsse

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n / 2) * log(2 * pi) - sum(log(obsse)) - 0.5 * SSR
}



# outlier handling functions
#####################################

# no outlier handling for Markov chians
mcmc_outlier_none <- function(., j) {
  if (j == .$wpars$mcmc_maxiter) print('No option was chosen to identify outlier Markov chains.')
}


# function that detects and corrects outlier Markov chains using the Inter Quartile-Range (IQR) statistic
mcmc_outlier_iqr <- function(., j) {

  counter <- j / .$wpars$mcmc_check_iter

  # extract last 50% of samples of each chain
  sbst <- .$dataf$pars_lklihood[1:.$wpars$mcmc_chains, (ceiling(j/2)):j]

  # take the mean of the log of the posterior densities and store in omega
  # IMPORANT: all current likelihood function options already return log-likelihood
  #           so it's not necessary to take the log of sbst components here
  #           but this may change in the future with different likelihood functions
  for (ii in 1:.$wpars$mcmc_chains) .$dataf$omega[ii, counter] <- mean(sbst[ii, ])

  # determine upper and lower quantiles of N different chains
  q1 <- quantile(.$dataf$omega[1:.$wpars$mcmc_chains, counter], prob = 0.25, type = 1)
  q3 <- quantile(.$dataf$omega[1:.$wpars$mcmc_chains, counter], prob = 0.75, type = 1)

  # compute IQR statistic
  iqr <- q3 - q1

  # determine which chains are outliers
  outliers <- which(.$dataf$omega[ , counter] < (q1 - 2 * iqr))

  # if outlier chains are detected
  if (length(outliers) > 0) {

    print(paste0('Outlier chain detected. Chain ', outliers, ' at iteration ', j))

    # replace outlier(s) by randomly choosing from the remaining chains
    replace_idx <- rep(0, length(outliers))
    for (qq in 1:length(outliers)) {
      # check: make sure no replacement chains are equal to any other elements of outliers
      while ((replace_idx[qq] == 0) | (replace_idx[qq] %in% outliers)) {
        replace_idx[qq] <- ceiling(runif(1, min = 0, max = 1) * .$wpars$mcmc_chains)
      }
    }

    # check: make sure 2 or more outliers aren't being replaced by the same randomly-chosen chain
    while (any(duplicated(replace_idx))) {
      repeat_idx <- which(duplicated(replace_idx))
      for (qq in 1:length(repeat_idx)) {
          replace_idx[repeat_idx[qq]] <-  ceiling(runif(1, min = 0, max = 1) * .$wpars$mcmc_chains)
      }
    }

    # replace outlier chain(s)
    .$dataf$pars_array[1:.$mcmc$d, outliers, 1:j] <- .$dataf$pars_array[1:.$mcmc$d, replace_idx, 1:j]

    # replace likelihood history for next iqr calculation
    #.$dataf$pars_lklihood[outliers, j] <- .$dataf$pars_lklihood[replace_idx, j]
    .$dataf$pars_lklihood[outliers, 1:j] <- .$dataf$pars_lklihood[replace_idx, 1:j]
  }

  # IMPORTANT: identifying and correcting outliers should only be done during burn-in
  #            because it violates the balance of sampled chains and destroys reversibility
  #            if outlier chain is detected, discard all previous sample history append
  #            and apply another burn-in period before generating posterior moments

  # future work: if outlier is detected, throw out all previous MCMC samples
  #              this is currently done manually in post-processing
  #              maybe automate this moving forward?

}


### END ###
