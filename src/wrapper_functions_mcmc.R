################################
#
# Wrapper functions for MCMC runs
#
# AJohnson, AWalker July 2019
#
################################



# prior distribution functions
################################

# initialize chains with uniform distributions
# APW: code can be simplified I think, fix, add in code snippets that are more consistent with the existing method 
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

  #print('')
  #print(.$dynamic$pars)
  #print('')
  
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
    .$dataf$pars[jj, 1:n] <- runif(n, min = min_vals[jj], max = max_vals[jj])
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

  if(.$wpars$mcmc$prior=='uniform') {
    min_vals <- vals[seq(1, length(vals), 2)]
    max_vals <- vals[seq(2, length(vals), 2)]
  }

  if(.$wpars$mcmc$prior=='normal') {
    min_vals  <- vals[seq(1, length(vals), 4)]
    max_vals  <- vals[seq(2, length(vals), 4)]
  }

  if(.$wpars$mcmc$prior=='none') {
    # future work: need to figure this out
    #              but this depends on the restart process
  }

  .$mcmc$boundary_min <- min_vals
  .$mcmc$boundary_max <- max_vals

  # if prior initialized with normal distribution, do boundary check
  if(.$wpars$mcmc$prior=='normal') {
    d <- length(unlist(.$dynamic$pars, recursive = T)) / 4
    for(ii in 1:.$wpars$mcmc$chains) {
      for(jj in  1:d) {
        .$mcmc_bdry_handling(j=1, ii=ii, jj=jj )
      }
    }
  }
}


# no boundary handling: should be chosen when search space is not theoretically restricted
mcmc_bdry_handling_none <- function(., j, ii, jj) {
  if ((j == .$wpars$mcmc$maxiter) & (ii == .$wpars$mcmc$chains) & (jj == .$mcmc$pars_n)) {
    print('No option was chosen for MCMC boundary handling.')
  }
}


# restrict parameter proposals that are beyond the boundary (reset them back to the min/max bound)
mcmc_bdry_handling_bound <- function(., j, ii, jj ) {

  # if outside bound of parameter space, restrict proposal value to corresponding dimension min
  if (.$dataf$pars[jj, ii] < .$mcmc$boundary_min[jj]) .$dataf$pars[jj, ii] <- .$mcmc$boundary_min[jj]

  # if outside bound of parameter space, restrict proposal value to corresponding dimension max
  else if (.$dataf$pars[jj, ii] > .$mcmc$boundary_max[jj]) .$dataf$pars[jj, ii] <- .$mcmc$boundary_max[jj]
}


# reflect parameter proposals that are beyond the min/max values back accros the boundary
mcmc_bdry_handling_reflect <- function(., j, ii, jj ) {

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
mcmc_bdry_handling_fold <- function(., j, ii, jj ) {

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



# likelihood functions
################################

# expects model output to be probability
# - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {

  # return log density
  log(.$dataf$out)
}


# standard error probability density function with i.i.d. error residuals
f_proposal_lklihood_ssquared <- function(.) {

  # number of measured data points
  obs_n <- length(.$dataf$obs)

  # calculate error residual
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2) )

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2) * log(SSR)
}


# standard error probability density function with i.i.d. error residuals
# - incorporates measurement errors (unlike "ssquared" option)
f_proposal_lklihood_ssquared_se <- function(.) {

  # read in measurement error; remove zeros from measurement error
  sspos <- which(.$dataf$obsse>1e-9)

  # number of measured data points (that do not have zero uncertainty)
  obs_n <- length(sspos)

  # observed error (heteroscedastic and homoscedastic options)
  obsse <- if(.$wpars$mcmc$homosced) rep(mean(.$dataf$obsse[sspos]), obs_n)
           else                      .$dataf$obsse[sspos]

  # calculate error residual
  error_residual_matrix <- ( t(.$dataf$out)[sspos,]-.$dataf$obs[sspos] ) / obsse

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2) )

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2) * log(2*pi) - sum(log(obsse)) - 0.5*SSR
}



# DREAM MCMC functions
################################

# initialisation of DREAM algorithm
init_mcmc_dream <- function(.) {

  # number of parameters being estimated
  .$mcmc$pars_n   <- dim(.$dataf$pars)[1]
  .$mcmc$sd_state <- numeric(.$mcmc$pars_n)

  # preallocate memory space for algorithmic variables
  # APW: note lp and chains are the same
  .$mcmc$current_state <- matrix(0, nrow=.$mcmc$pars_n, ncol=.$wpars$mcmc$chains )
  .$mcmc$jump          <- matrix(0, nrow=.$mcmc$pars_n, ncol=.$wpars$mcmc$chains )
  .$mcmc$lambda        <- numeric(.$wpars$mcmc$chains)

  # initialise crossover variables if not a restart
  if(!.$wpars$parsinit_read) {
    .$mcmc$CR_counter      <- numeric(.$wpars$mcmc$n_CR)
    .$mcmc$jump_delta_norm <- numeric(.$wpars$mcmc$n_CR)
    .$mcmc$p_CR            <- numeric(.$wpars$mcmc$n_CR)
    # initial probability of each crossover value
    .$mcmc$p_CR[]          <- 1/.$wpars$mcmc$n_CR
    .$mcmc$CR              <- numeric(.$wpars$mcmc$chains)
  }
}


# generate proposal using DREAM algorithm (Vrugt et al. 2011)
# - gamma = gamma in V2011
# - lambda = e in V2011
# - chain_pairs_n = delta in V2011
# - length(jump_pars_ss) = dprime in V2011
# - CR = m in V2011
# - CR/nCR = CR in V2011
# - p_CR = p_m in V2011
proposal_generate_mcmc_dream <- function(.,j) {

  # debugging
  #print(paste0('j = ',j))

  # initialise
  # - continuous uniform random values between -c_rand and c_rand
  .$mcmc$jump[]          <- 0
  .$mcmc$current_state[] <- .$dataf$pars_array[,,j-1]
  .$mcmc$lambda[]        <- runif(.$dataf$lp, -.$wpars$mcmc$c_rand, .$wpars$mcmc$c_rand )

  # if adapting crossover values, compute standard deviation of each parameter/dimension
  if(.$wpars$mcmc$adapt_pCR) {
    .$mcmc$sd_state[] <- apply(.$mcmc$current_state, 1, sd )
    .$mcmc$sd_state[.$mcmc$sd_state==0] <- 1e-9
  }

  # create proposals for each chain
  for (ii in 1:.$dataf$lp) {

    # determine chain pairs used to calculate each jump
    # APW: why choose 1 value w replacement? Maybe this should be outside of the chain loop? Given it's inside the loop it will be with replacement
    chain_pairs_n  <- sample(1:.$wpars$mcmc$chain_delta, 1, T )
    chain_pairs_ss <- t(sapply(1:chain_pairs_n, function(v) sample((1:.$dataf$lp)[-ii],2,F) ))

    # select crossover value
    # - weighted sample from multinomial distribution
    # - replacement relevant if this gets moved outside of chain loop 
    .$mcmc$CR[ii]  <- sample(1:.$wpars$mcmc$n_CR, 1, T, .$mcmc$p_CR )

    # determine which parameters will "crossover" (i.e. how many dimensions are sampled/updated jointly)
    zz             <- runif(.$mcmc$pars_n)
    jump_pars_ss   <- which(zz < (.$mcmc$CR[ii]/.$wpars$mcmc$n_CR) )
    if(length(jump_pars_ss)==0) jump_pars_ss <- which.min(zz)

    # jump rate / scaling factor
    gamma_d        <- 2.38 / sqrt(2*chain_pairs_n*length(jump_pars_ss))
    gamma          <- sample(c(gamma_d,1), 1, T, c(1-.$wpars$mcmc$p_gamma, .$wpars$mcmc$p_gamma ))

    # compute 'jump' for params to be updated/crossover (differential evolution)
    chain_diff     <- apply(.$mcmc$current_state[jump_pars_ss,chain_pairs_ss[,1],drop=F], 1, sum ) - 
                      apply(.$mcmc$current_state[jump_pars_ss,chain_pairs_ss[,2],drop=F], 1, sum )
    .$mcmc$jump[jump_pars_ss,ii] <- .$wpars$mcmc$c_ergod*rnorm(length(jump_pars_ss)) + (1+.$mcmc$lambda[ii])*gamma*chain_diff 
    .$dataf$pars[,ii]            <- .$mcmc$current_state[,ii] + .$mcmc$jump[,ii]

    # boundary handling
    # APW: prob can happen outside of chain loop 
    for(jj in 1:.$mcmc$pars_n) .$mcmc_bdry_handling(j=j, ii=ii, jj=jj )

  # chain loop
  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_mcmc_dream <- function(.,j) {

  # likelihoods of proposed and current states
  prop_lklihood <- .$proposal_lklihood()
  curr_lklihood <- .$dataf$pars_lklihood[,j-1]

  # iterate through chains
  for(ii in 1:.$dataf$lp) {

    # Metropolis acceptance probability
    alpha  <- min(1, exp(prop_lklihood[ii]-curr_lklihood[ii]) )
    accept <- alpha>runif(1)

    # APW: also not strictly related to specific acceptance function
    if(accept) {
      .$dataf$pars_array[,ii,j]   <- .$dataf$pars[,ii]
      .$dataf$pars_lklihood[ii,j] <-  prop_lklihood[ii]
      .$dataf$out_mcmc[,ii,j]     <- .$dataf$out[ii,] 
    } else {
      .$dataf$pars_array[,ii,j]   <- .$dataf$pars_array[,ii,j-1]
      .$dataf$pars_lklihood[ii,j] <- curr_lklihood[ii]
      .$dataf$out_mcmc[,ii,j]     <- .$dataf$out_mcmc[,ii,j-1]
    }

    # APW: not directly related to acceptance, could move to DREAM run function, fix
    # compute squared normalized jumping distance
    # - count selected crosover values
    # APW: potentially inefficient, does this need calculating for each chain or can it be simultaneous?, fix
    # - jump_delta_norm = Delta_m V2011
    # - CR_counter = L_m V2011
    if(.$wpars$mcmc$adapt_pCR) {
      summation <- sum(((.$dataf$pars_array[,ii,j] - .$mcmc$current_state[,ii]) / .$mcmc$sd_state)^2 )
      .$mcmc$jump_delta_norm[.$mcmc$CR[ii]] <- .$mcmc$jump_delta_norm[.$mcmc$CR[ii]] + summation
      .$mcmc$CR_counter[.$mcmc$CR[ii]]      <- .$mcmc$CR_counter[.$mcmc$CR[ii]] + 1
    } 
  }

  # APW: not directly related to acceptance, could move to DREAM run function, fix
  # update the selection probability of crossover probabilities/values
  # debugging
  if(.$wpars$mcmc$adapt_pCR) {
    print('')
    print('.$mcmc$CR'); print(.$mcmc$CR)
    print('.$mcmc$CR_counter'); print(.$mcmc$CR_counter)
    print('.$mcmc$sd_state'); print(.$mcmc$sd_state)
    print('.$mcmc$jump_delta_norm'); print(.$mcmc$jump_delta_norm)
  }
  if(.$mcmc$adapt_pCR) {

    .$mcmc$p_CR[] <- .$mcmc$j_true*.$wpars$mcmc$chains * (.$mcmc$jump_delta_norm/.$mcmc$CR_counter) / sum(.$mcmc$jump_delta_norm)
    .$mcmc$p_CR[] <- .$mcmc$p_CR/sum(.$mcmc$p_CR)
  
    # debugging
    print('')
    print('adapt_pCR calculation')
    print('.$mcmc$j_true'); print(.$mcmc$j_true)
    print('.$mcmc$CR_counter'); print(.$mcmc$CR_counter)
    print('.$mcmc$jump_delta_norm'); print(.$mcmc$jump_delta_norm)
    print('.$mcmc$p_CR'); print(.$mcmc$p_CR)

    if(.$mcmc$j_true==.$wpars$mcmc$CR_burnin) {
      print('',quote=F)
      print('',quote=F)
      print('Adapted selection probabilities of crossover values:',quote=F); print(.$mcmc$p_CR,quote=F)
      .$wpars$mcmc$adapt_pCR[] <- .$mcmc$adapt_pCR[] <- F
    }
  }
}



# outlier handling functions
#####################################

# IMPORTANT: identifying and correcting outliers should only be done during burn-in
#            because it violates the balance of sampled chains and destroys reversibility
#            if outlier chain is detected, discard all previous sample history, then append
#            and apply another burn-in period before generating posterior moments

# no outlier handling for Markov chains
mcmc_outlier_none <- function(.,j) {
  if(j==.$wpars$mcmc$maxiter) print('No option was chosen to identify outlier Markov chains.')
}


# function that detects and corrects outlier Markov chains using the Inter Quartile-Range (IQR) statistic
mcmc_outlier_iqr <- function(.,j) {

#  print('')
#  print('jstartburnin,jb50, j, check_ss:')
#  print(c(.$mcmc$j_start_burnin,.$mcmc$j_burnin50,j,.$mcmc$check_ss))

  # identify outlier chains 
  # - based on IQR across chains mean log posterior densities of last 50 % of burnin
  # IMPORTANT: all current likelihood function options return log-likelihood
  #            so it's not necessary to take the log lklihood here
  #            but may change in future with different likelihood functions
  .$dataf$omega[,.$mcmc$check_ss] <- apply(.$dataf$pars_lklihood[,.$mcmc$j_burnin50:j], 1, mean )
  q1q3     <- quantile(.$dataf$omega[1:.$wpars$mcmc$chains,.$mcmc$check_ss], prob=c(0.25,0.75), type=1 )
  iqr      <- q1q3[2]-q1q3[1]
  outliers <- which(.$dataf$omega[,.$mcmc$check_ss] < (q1q3[1]-2*iqr))

  # if outlier chains are detected
  # APW: could move this to main DREAM run functon as it is generic
  if (length(outliers)>0) {

    print('',quote=F)
    print(paste('Outlier chain(s) detected. Chain(s):', outliers, 'at iteration:', .$mcmc$j_true ), quote=F )

    # replace outlier chain(s) & likelihood history for next iqr calculation
    replace_ss <- sample((1:.$wpars$mcmc$chains)[-outliers], length(outliers) )
    .$dataf$pars_array[1:.$mcmc$pars_n,outliers,j] <- .$dataf$pars_array[1:.$mcmc$pars_n,replace_ss,j]
    .$dataf$pars_lklihood[outliers,j]         <- .$dataf$pars_lklihood[replace_ss,j]

    # restart burn-in
    .$mcmc$outlier_detected <- T
    .$mcmc$j_start_burnin   <- j + 1
    .$mcmc$j_burnin50       <- j
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
  iter_effective <- 2*(j-.$mcmc$j_burnin50)
  half_effective <- ceiling(iter_effective/2)

  if(iter_effective>0) {
    # within-chain variance of parameters
    x_bar        <- (2/(iter_effective-2)) * apply(.$dataf$pars_array[,,.$mcmc$j_burnin50:j], 1:2, sum )
    summation    <- numeric(.$mcmc$pars_n)
    for(jj in 1:.$mcmc$pars_n) summation[jj] <- sum((.$dataf$pars_array[jj,,.$mcmc$j_burnin50:j] - x_bar[jj,])^2 )
    W            <- 2/(.$wpars$mcmc$chains*(iter_effective-2))*summation

    # between-chain variance of parameters
    x_double_bar <- apply(x_bar, 1, sum ) / .$wpars$mcmc$chains 
    summation    <- apply((x_bar-x_double_bar)^2, 1, sum )
    B            <- (iter_effective/(2*(.$wpars$mcmc$chains-1)))*summation
  
    # parameter variance
    sigma_hat    <- ((iter_effective-2)/iter_effective)*W + (2/iter_effective)*B
  
    # R-statistic of Gelman and Rubin
    R_hat        <- sqrt(((.$wpars$mcmc$chains+1)/.$wpars$mcmc$chains)*(sigma_hat/W) - ((iter_effective-2)/(.$wpars$mcmc$chains*iter_effective)))
    R_hat_new    <- c(.$mcmc$j_true, j, iter_effective, R_hat )

  } else {
    R_hat_new    <- c(.$mcmc$j_true, j, 0, rep(NA,.$mcmc$pars_n) )
    R_hat        <- 'NA, outlier detected on final iteration'
  }

  # APW: could move this to main DREAM run functon as it is generic
  .$dataf$conv_check[,.$mcmc$check_ss] <- R_hat_new

  #print('',quote=F)
  #print(.$mcmc$t,quote=F)
  #print(.$mcmc$j_true,quote=F)

  if (j==.$wpars$mcmc$maxiter) {
    print('',quote=F)
    print('',quote=F)
    print(paste("At (final) iteration:", .$mcmc$j_true, ", R-statistic of Gelman and Rubin ="), quote=F )
    print(R_hat,quote=F)
    print('',quote=F)
    #print('Convergence criterion:',quote=F)
    #print(.$dataf$conv_check,quote=F)
    #print('',quote=F)
  }
}



# DEMC functions
################################

init_mcmc_demc <- function(.) NULL

# generate proposal using DE-MC algorithm
proposal_generate_mcmc_demc <- function(., j ) {

  # scaling factor
  # APW: can be calculated once I think, fix
  d          <- dim(.$dataf$pars)[1]
  gamma_star <- 2.38 / sqrt(d + d)

  # b-value should be small compared to width of target distribution; specifies range for drawn "randomization" value
  b_rand  <- 0.01
  uniform_r <- runif(1, min=(-b_rand), max=b_rand)

  # evaluate for each chain
  for(ii in 1:.$dataf$lp) {

    # randomly select two different numbers R1 and R2 unequal to j, from uniform distribution without replacement
    chain_pair <- sample((1:.$dataf$lp)[-.$dataf$lp], 2, F )

    # evaluate for each parameter
    for(jj in 1:d) {

      # generate proposal via Differential Evolution
      .$dataf$pars[jj,ii] <- .$dataf$pars_array[jj,ii,j-1] + uniform_r + 
        gamma_star*( .$dataf$pars_array[jj,chain_pair[1],j-1] - .$dataf$pars_array[jj,chain_pair[2],j-1] )

      # boundary handling 
      .$boundary_handling(ii=ii, jj=jj )
    }
  }
}


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_mcmc_demc <- function(., j, lklihood ) {

  # Metropolis ratio
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[,j-1])
  alpha        <- pmin(1, metrop_ratio)

  # evaluate for each chain
  for(ii in 1:.$dataf$lp) {
    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
    accept <- log(alpha[ii]) > log(runif(1, min = 0, max = 1))
    .$dataf$pars_array[,ii,j]   <- if(accept)          .$dataf$pars[,ii] else .$dataf$pars_array[,ii,j-1]
    .$dataf$pars_lklihood[ii,j] <- if(accept)          lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]
    .$dataf$out_mcmc[ii,,j]     <- if(accept | j == 1) .$dataf$out[ii,]  else .$dataf$out_mcmc[ii,,(j-1)]
  }
}



### END ###
