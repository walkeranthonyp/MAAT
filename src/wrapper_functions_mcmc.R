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
  # APW: commenting out where variables already declared and do do have a varaible extent 
  #.$mcmc$t         <- numeric(1)
  #.$mcmc$d_star    <- numeric(1)
  #.$mcmc$CR_burnin <- numeric(1)
  .$mcmc$sd_state  <- numeric(.$mcmc$d)
  .$mcmc$L         <- numeric(.$wpars$mcmc$n_CR)
  .$mcmc$del       <- numeric(.$wpars$mcmc$n_CR)
  .$mcmc$p_CR      <- numeric(.$wpars$mcmc$n_CR)
  .$mcmc$m         <- numeric(.$wpars$mcmc$chains)

  # if adapting selection of crossover probabilties
  if(.$wpars$mcmc$adapt_pCR) {

    # initialize crossover variables
    .$mcmc$t      <- 1
    .$mcmc$CR     <- 0
    .$mcmc$L[]    <- 0
    .$mcmc$d_star <- .$mcmc$d

    # initial probability of each crossover value
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc$n_CR

  # if not adapting p_CR
  } else {

    # initialize crossover probabilities
    .$mcmc$CR   <- numeric(.$wpars$mcmc$n_CR)
    .$mcmc$CR[] <- 1:.$wpars$mcmc$n_CR / .$wpars$mcmc$n_CR

    # initialize selection probability of crossover values
    .$mcmc$p_CR[] <- 1 / .$wpars$mcmc$n_CR

  }

  # burn-in period for adapting crossover selection probabilities
  #.$mcmc$CR_burnin <- ceiling(.$wpars$mcmc$CR_burnin * .$wpars$mcmc$maxiter)

  # index of chains for Differential Evolution
  for (ii in 1:.$dataf$lp) .$mcmc$R[ii, ] <- setdiff(1:.$dataf$lp, ii)
}


# generate proposal using DREAM algorithm
proposal_generate_mcmc_dream <- function(., j ) {

  # debugging
  #print(paste0('j = ',j))

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1])

  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[] <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)

  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[] <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc$c_rand, .$wpars$mcmc$c_rand), .$dataf$lp)

  # if not adapting crossover values, compute standard deviation of each dimension/parameter
  if(!(.$wpars$mcmc$adapt_pCR)) {

    .$mcmc$sd_state[] <- apply(.$mcmc$current_state, 1, sd)

    # replace any 0's in standard deviation array with 1e-9 to avoid division by 0
    idx <- which(.$mcmc$sd_state == 0)
    .$mcmc$sd_state[idx] <- 1e-9
  }

  # create proposals
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc$delta, 1, replace = T)

    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    if (.$wpars$mcmc$adapt_pCR) {
      # modify each dimension with probability CR each time a proposal vector is generated

      # generate crossover probability
      .$mcmc$m[ii] <- .$generate_CR()

      # calculate jump rate (scaling factor)
      gamma_d <- 2.38 / sqrt(2 * D * .$mcmc$d_star)

      # when gamma = 1, jump between different modes of the posterior (this is approx every 5 iterations with default p_gamma = 0.2)
      gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc$p_gamma, .$wpars$mcmc$p_gamma))

      # compute jump differential evolution of ii-th chain
      .$mcmc$jump[1:.$mcmc$d, ii] <- .$wpars$mcmc$c_ergod * rnorm(.$mcmc$d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[1:.$mcmc$d, a] - .$mcmc$current_state[1:.$mcmc$d, b]), dim = 1)

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
      id <- sample(1:.$wpars$mcmc$n_CR, 1, replace = T, prob = .$mcmc$p_CR)

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

      # when gamma = 1, jump between different modes of the posterior
      # approx every 5 iterations with default p_gamma = 0.2
      gamma <- sample(c(gamma_d, 1), size = 1, replace = T, prob = c(1 - .$wpars$mcmc$p_gamma, .$wpars$mcmc$p_gamma))

      # compute jump differential evolution of ii-th chain
      .$mcmc$jump[A, ii] <- .$wpars$mcmc$c_ergod * rnorm(.$mcmc$d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[A, a] - .$mcmc$current_state[A, b]), dim = 1)

      # compute proposal of ii-th chain
      .$dataf$pars[1:.$mcmc$d, ii] <- .$mcmc$current_state[1:.$mcmc$d, ii] + .$mcmc$jump[1:.$mcmc$d, ii]

    }

    # call boundary handling function
    for (jj in 1:.$mcmc$d) .$mcmc_bdry_handling(j=j, ii = ii, jj = jj)
  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_mcmc_dream <- function(., j, lklihood ) {

  # likelihood of current state
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[, j-1]

  # iterate through chains
  for (ii in 1:.$dataf$lp) {

    # Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[ii] - .$mcmc$p_state[ii] ))

    if(alpha>runif(1, min=0, max=1 )) {

      # accept proposal
      accept <- TRUE

      .$mcmc$current_state[1:.$mcmc$d,ii] <- .$dataf$pars[1:.$mcmc$d,ii]
      .$mcmc$p_state[ii]                  <- lklihood[ii]

      # append accepted proposal and probability density to storage arrays
      .$dataf$pars_array[1:.$mcmc$d,ii,j] <- .$mcmc$current_state[1:.$mcmc$d,ii]
      .$dataf$pars_lklihood[ii,j]         <- .$mcmc$p_state[ii]

    } else {

      # reject proposal
      accept <- FALSE

      # reset jump back to zero
      .$mcmc$jump[1:.$mcmc$d,ii] <- 0

      # repeat previous accepted proposal and probability density in storage arrays
      .$dataf$pars_array[1:.$mcmc$d,ii,j] <- .$dataf$pars_array[1:.$mcmc$d,ii,j-1]
      .$dataf$pars_lklihood[ii,j]         <- .$dataf$pars_lklihood[ii,j-1]

    }

    # compute squared normalized jumping distance
    #if (.$wpars$mcmc$adapt_pCR & (.$mcmc$t < .$mcmc$CR_burnin)) .$calc_del(j = j, ii = ii)
    if(.$wpars$mcmc$adapt_pCR & .$mcmc$CR_burnin) .$calc_del(j=j, ii=ii )

    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted 
    #      as it adds the output from the rejected proposal
    # write output (model evaluations)
    .$dataf$out_mcmc[ii,,j] <- if(accept | j==1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-1)]

  }

  #if (.$wpars$mcmc$adapt_pCR & (.$mcmc$t < .$mcmc$CR_burnin)) {
  if(.$wpars$mcmc$adapt_pCR) {

    # debugging
    #print('if-statement triggered calling .$adapt_pCR()')
    #print('.$mcmc$t'); print(.$mcmc$t)
    #print('.$mcmc$CR_burnin'); print(.$mcmc$CR_burnin)

    #if(.$mcmc$t < .$mcmc$CR_burnin) {
    if(.$mcmc$CR_burnin) { 

      # update the selection probability of crossover probabilities/values
      # APW: update the jump_n selection probabilities
      .$adapt_pCR()

      if(.$mcmc$t==.$wpars$mcmc$CR_burnin) {
        print('Adapted selection probabilities of crossover values = '); print(.$mcmc$p_CR)
        .$mcmc$CR_burnin <- F
      }
    } 
    #.$mcmc$t <- .$mcmc$t + 1

    #} else if (.$wpars$mcmc$adapt_pCR & .$mcmc$t == .$mcmc$CR_burnin) {

    #print('Adapted selection probabilities of crossover values = '); print(.$mcmc$p_CR)

    #.$mcmc$t <- .$mcmc$t + 1

    #} else if (.$wpars$mcmc$adapt_pCR & (.$mcmc$t > .$mcmc$CR_burnin)) {
    #  .$mcmc$t <- .$mcmc$t + 1
    #}

    .$mcmc$t <- .$mcmc$t + 1
  }
}



# adaptive p_CR functions
################################

# function that generates/updates crossover values based on current probabilities
generate_CR <- function(.) {

  # debugging
  #print('in generate_CR function')
  #print('.$wpars$mcmc$n_CR'); print(.$wpars$mcmc$n_CR)
  #print('.$mcmc$p_CR'); print(.$mcmc$p_CR)

  # sample m from numbers 1,...,n_CR using multinomial distribution (with probabilities p_CR)
  m <- sample(1:.$wpars$mcmc$n_CR, size = 1, replace = T, prob = .$mcmc$p_CR)

  # set crossover probability/value
  .$mcmc$CR <- m / .$wpars$mcmc$n_CR

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

  # debugging
  #print('calc_del function being called')
  #print('.$mcmc$sd_state'); print(.$mcmc$sd_state)
  #print('summation = '); print(summation)
  #print('.$dataf$pars_array[,ii,j]'); print(.$dataf$pars_array[,ii,j])
  #print('.$dataf$pars_array[,ii,j-1]'); print(.$dataf$pars_array[,ii,j-1])

}


# function that adapts crossover probabilities
adapt_pCR <- function(.) {

  # debugging
  #print('adapt_pCR function being called')
  #print('.$mcmc$L'); print(.$mcmc$L)
  #print('.$wpars$mcmc$chains'); print(.$wpars$mcmc$chains)
  #print('.$mcmc$del'); print(.$mcmc$del)
  # t is re-initialized to 1 in a restart
  # .$mcmc$ del is zero though

  # update the probability of the different crossover values being selected
  for (qq in 1:.$wpars$mcmc$n_CR) {
    # numerical check for divide by zero
    if (.$mcmc$L[qq] != 0) {
      .$mcmc$p_CR[qq] <- .$mcmc$t * .$wpars$mcmc$chains * (.$mcmc$del[qq] / .$mcmc$L[qq]) / sum(.$mcmc$del)
    }
  }

  # normalize
  .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
}



# prior distribution functions
################################

# initialize chains with uniform distributions
mcmc_prior_uniform <- function(.) {

  # IMPORTANT: when initializing parameters in the init file
  #            it needs to be in the following format:
  #            parameter = 'list(min = value, max = value)'
  #            if there is a nested list
  #            it needs to be in the following format:
  #            list = 'list(
  #                    parameter 1 = list(min = value, max = value),
  #                    parameter 2 = list(min = value, max = value)
  #                   )'

  # future work: re-structure prior functions to be less dependent on the order
  #              of terms listed in the init file
  #              (then can also restructure boundary_handling_set)

  # number of Markov chains
  n <- .$wpars$mcmc$chains

  # number of parameters (dimensionality of parameter space)
  d <- length(unlist(.$dynamic$pars, recursive = T)) / 2

  # determine minimums and maximums for parameters
  dynamic_pars_un <- unlist(.$dynamic$pars)
  vals <- sapply(dynamic_pars_un, function(x) x[[1]])
  max_vals <- vals[seq(2, length(vals), 2)]
  min_vals <- vals[seq(1, length(vals), 2)]

  # draw priors from uniform distribution to create pars / proposal matrix
  .$dataf$pars <- matrix(0, nrow = d, ncol = n)
  for (jj in 1:d) {
    .$dataf$pars[jj, 1:n] <- runif(n, min = min_valsf[jj], max = max_vals[jj])
  }

  # assign parameter names
  row_names <- gsub(pattern = '.min', replacement = '', names(dynamic_pars_un))
  row_names <- row_names[seq(1, length(row_names), 2)]
  rownames(.$dataf$pars) <- row_names
}


# initialize chains with normal distributions
mcmc_prior_normal <- function(.) {

  # IMPORTANT: when initializing parameters in the init file
  #            it needs to be in the following format:
  #            parameter = 'list(min = value, max = value, mean = value, sd = value)'
  #            if there is a nested list
  #            it needs to be in the following format:
  #            list = 'list(
  #                    parameter 1 = list(min = value, max = value, mean = value, sd = value),
  #                    parameter 2 = list(min = value, max = value, mean = value, sd = value)
  #                   )'

  # number of Markov chains
  n <- .$wpars$mcmc$chains

  # number of parameters (dimensionality of parameter space)
  d <- length(unlist(.$dynamic$pars, recursive = T)) / 4

  # determine minimums and maximums for parameters
  dynamic_pars_un <- unlist(.$dynamic$pars)
  vals <- sapply(dynamic_pars_un, function(x) x[[1]])
  mean_vals <- vals[seq(3, length(vals), 4)]
  sd_vals   <- vals[seq(4, length(vals), 4)]

  # draw priors from uniform distribution to create pars / proposal matrix
  .$dataf$pars <- matrix(0, nrow = d, ncol = n)
  for (jj in 1:d) {
    .$dataf$pars[jj, 1:n] <- rnorm(n, mean = mean_vals[jj], sd = sd_vals[jj])
  }

  # assign parameter names
  row_names <- gsub(pattern = '.min', replacement = '', names(dynamic_pars_un))
  row_names <- row_names[seq(1, length(row_names), 4)]
  rownames(.$dataf$pars) <- row_names

  # future work: add check to see if mean and sd are not specified
  #              and if not, make the mean the median of the parameter range
  #              and make it 2-3 standard deviations to the boundary
}


# option for initializing chains when re-starting MCMC algorithm
mcmc_prior_none <- function(.) {

  # use the last accepted proposal from previous MCMC run
  # future work: ideally would not like to hardcode this
  #              currently just copying and pasting

  #.$dataf$pars <- matrix(c(127.1412287, 85.6479380,  0.9263817,
  #                         151.1674128, 267.4227537, 0.9635282,
  #                         153.8731369, 201.6820169, 0.9468014,
  #                         153.7184214, 254.7146791, 0.9295626,
  #                         155.9559049, 251.6358836, 0.9570725,
  #                         140.8245361, 223.8698299, 0.9495525,
  #                         133.0800970, 249.0815415, 0.9322175),
  #                        ncol = .$wpars$mcmc$chains)

  #rownames(.$dataf$pars) <- c('leaf.atref.vcmax', 'leaf.atref.jmax', 'leaf.theta_col_cj')
  print('')
  print('dataf$pars from MCMC run:')
  print(.$dataf$pars)
}



# boundary handling functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {

  dynamic_pars_un <- unlist(.$dynamic$pars)
  vals <- sapply(dynamic_pars_un, function(x) x[[1]])

  if (.$wpars$mcmc$prior == 'uniform') {
    min_vals <- vals[seq(1, length(vals), 2)]
    max_vals <- vals[seq(2, length(vals), 2)]
  }

  if (.$wpars$mcmc$prior == 'normal') {
    min_vals  <- vals[seq(1, length(vals), 4)]
    max_vals  <- vals[seq(2, length(vals), 4)]
  }

  if (.$wpars$mcmc$prior == 'none') {
    # future work: need to figure this out
    #              but this depends on the restart process
  }

  .$mcmc$boundary_min <- min_vals
  .$mcmc$boundary_max <- max_vals

  # if prior initialized with normal distribution, do boundary check
  if(.$wpars$mcmc$prior == 'normal') {
    d <- length(unlist(.$dynamic$pars, recursive = T)) / 4
    for (ii in 1:.$wpars$mcmc$chains) {
      for (jj in  1:d) {
        .$mcmc_bdry_handling(j = 1, ii = ii, jj = jj)
      }
    }
  }
}


# no boundary handling: should be chosen when search space is not theoretically restricted
mcmc_bdry_handling_none <- function(., j, ii, jj) {
  if ((j == .$wpars$mcmc$maxiter) & (ii == .$wpars$mcmc$chains) & (jj == .$mcmc$d)) {
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
mcmc_converge_none <- function(.,j) {
  if(j==.$wpars$mcmc$maxiter) print('No option was chosen to test for MCMC convergence.')
}


# subroutine calculating the R-statistic of Gelman and Rubin (convergence diagnostic)
mcmc_converge_Gelman_Rubin <- function(.,j) {

  # effective number of iterations since burn-in began
  # ALJ: this needs to equivalent to the total number of "usable" time samples in each chain?
  iter_effective <- j-.$mcmc$j_start_burnin+1
  half_effective <- iter_effective/2

  # within-chain variance
  #W <- numeric(.$mcmc$d)
  #x_bar <- matrix(0, nrow=.$mcmc$d, ncol=.$wpars$mcmc$chains)
  x_bar     <- (2/(iter_effective-2)) * apply(.$dataf$pars_array[,,half_effective:iter_effective], 1:2, sum)
  summation <- numeric(.$mcmc$d)
  for (jj in 1:.$mcmc$d) summation[jj] <- sum((.$dataf$pars_array[jj,,half_effective:iter_effective]-x_bar[jj,])^2)
  W         <- 2/(.$wpars$mcmc$chains*(iter_effective-2))*summation

  # between-chain variance
  #B <- numeric(.$mcmc$d)
  #x_double_bar <- numeric(.$mcmc$d)
  x_double_bar <- (1/.$wpars$mcmc$chains) * apply(x_bar[,], 1, sum)
  summation    <- apply(((x_bar[,]-x_double_bar[])^2), 1, sum)
  B            <- (iter_effective/(2*(.$wpars$mcmc$chains-1)))*summation

  # estimate variance of jjth parameter of target distribution
  sigma_hat <- ((iter_effective-2)/iter_effective)*W + (2/iter_effective)*B

  # R-statistic of Gelman and Rubin
  R_hat <- sqrt(((.$wpars$mcmc$chains+1)/.$wpars$mcmc$chains)*(sigma_hat/W) - ((iter_effective-2)/(.$wpars$mcmc$chains*iter_effective)))

  # append corresponding effective iteration number to R_hat vector
  R_hat_new <- append(R_hat, iter_effective, after=0 )
  .$dataf$conv_check[,.$wpars$mcmc$check_ss] <- R_hat_new

  if (j==.$wpars$mcmc$maxiter) {
    print(paste0("At iteration ", j, ", R-statistic of Gelman and Rubin = ")); print(R_hat)
    print('Convergence criterion:'); print(.$dataf$conv_check)
  }
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
  obsse <- if(.$wpars$mcmc$homosced) rep(mean(.$dataf$obsse[sspos]), obs_n)
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

# no outlier handling for Markov chains
mcmc_outlier_none <- function(.,j) {
  if(j==.$wpars$mcmc$maxiter) print('No option was chosen to identify outlier Markov chains.')
}


# function that detects and corrects outlier Markov chains using the Inter Quartile-Range (IQR) statistic
mcmc_outlier_iqr <- function(.,j) {

  #counter <- j/.$wpars$mcmc$check_iter

  #.$mcmc$j_burnin50 <-
  #   if(.$mcmc$outlier_detected|!.$wpars$parsinit_read) .$mcmc$j_start_burnin + ceiling(j/2)
  #   else                                                     j - ceiling((j+.$wpars$mcmc$start_iter-1)/2)

  # extract last 50% of samples of each chain
  #sbst <- .$dataf$pars_lklihood[,.$mcmc$j_burnin50:j]

  # take the mean of the log of the posterior densities and store in omega
  # IMPORTANT: all current likelihood function options already return log-likelihood
  #            so it's not necessary to take the log of sbst components here
  #            but this may change in the future with different likelihood functions
  #for (ii in 1:.$wpars$mcmc$chains) .$dataf$omega[ii,.$wpars$mcmc$check_ss] <- mean(sbst[ii, ])
  .$dataf$omega[,.$wpars$mcmc$check_ss] <- apply(.$dataf$pars_lklihood[,.$mcmc$j_burnin50:j], 1, mean)

  # determine upper and lower quantiles of N different chains
  q1 <- quantile(.$dataf$omega[1:.$wpars$mcmc$chains,.$wpars$mcmc$check_ss], prob = 0.25, type = 1)
  q3 <- quantile(.$dataf$omega[1:.$wpars$mcmc$chains,.$wpars$mcmc$check_ss], prob = 0.75, type = 1)

  # compute IQR statistic
  iqr <- q3-q1

  # determine which chains are outliers
  outliers <- which(.$dataf$omega[,.$wpars$mcmc$check_ss]<(q1-2*iqr))

  # if outlier chains are detected
  if (length(outliers) > 0) {

    print(paste0('Outlier chain detected. Chain ', outliers, ' at iteration ', j))

    # replace outlier(s) by randomly choosing from the remaining chains
    replace_idx <- rep(0, length(outliers))
    for (qq in 1:length(outliers)) {
      # check: make sure no replacement chains are equal to any other elements of outliers
      while ((replace_idx[qq] == 0) | (replace_idx[qq] %in% outliers)) {
        replace_idx[qq] <- ceiling(runif(1, min = 0, max = 1) * .$wpars$mcmc$chains)
      }
    }

    # check: make sure 2 or more outliers aren't being replaced by the same randomly-chosen chain
    while (any(duplicated(replace_idx))) {
      repeat_idx <- which(duplicated(replace_idx))
      for (qq in 1:length(repeat_idx)) {
          replace_idx[repeat_idx[qq]] <-  ceiling(runif(1, min = 0, max = 1) * .$wpars$mcmc$chains)
      }
    }

    # replace outlier chain(s)
    #.$dataf$pars_array[1:.$mcmc$d, outliers, 1:j] <- .$dataf$pars_array[1:.$mcmc$d, replace_idx, 1:j]
    .$dataf$pars_array[1:.$mcmc$d,outliers,j] <- .$dataf$pars_array[1:.$mcmc$d,replace_idx,j]

    # replace likelihood history for next iqr calculation
    #.$dataf$pars_lklihood[outliers, 1:j] <- .$dataf$pars_lklihood[replace_idx, 1:j]
    .$dataf$pars_lklihood[outliers,j] <- .$dataf$pars_lklihood[replace_idx,j]

    # restart burn-in
    .$mcmc$outlier_detected <- T
    .$mcmc$j_start_burnin   <- j
  }

  # IMPORTANT: identifying and correcting outliers should only be done during burn-in
  #            because it violates the balance of sampled chains and destroys reversibility
  #            if outlier chain is detected, discard all previous sample history, then append
  #            and apply another burn-in period before generating posterior moments
}



### END ###
