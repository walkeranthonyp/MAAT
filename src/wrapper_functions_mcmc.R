################################
#
# Wrapper functions for MCMC runs
#
# AJohnson, AWalker July 2019
#
################################



# boundary handling functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {
  # number of samples to be used in boundary handling
  n <- 1000
  boundary_sample     <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text=cs)) )
  .$mcmc$boundary_min <- unlist(lapply(boundary_sample,min))
  .$mcmc$boundary_max <- unlist(lapply(boundary_sample,max))
  rm(boundary_sample)
}



# no boundary handling: should be chosen when search space is not theoretically restricted
boundary_handling_none <- function(., j, ii, jj ) {
  if((j==.$wpars$mcmc$maxiter) & (ii==.$wpars$mcmc$chains) & (jj==.$mcmc$pars_n)) {
    print('No option was chosen for MCMC boundary handling.')
  }
}


# restrict parameter proposals that are beyond the boundary to the boundary 
boundary_handling_bound <- function(.) {

  pm   <- .$dataf$pars
  bmax <- .$mcmc$boundary_max 
  bmin <- .$mcmc$boundary_min

  .$dataf$pars[] <- 
    t(sapply(1:dim(pm)[1], function(p, pm2=pm[p,], min2=bmin[p], max2=bmax[p])
       ifelse(pm2<min2, min2, 
       ifelse(pm2>max2, max2, 
              pm2) 
       )))
}

# reflect parameter proposals beyond the min/max values back across the boundary
boundary_handling_reflect <- function(.) {

  pm   <- .$dataf$pars
  bmax <- .$mcmc$boundary_max 
  bmin <- .$mcmc$boundary_min

  # vectorised ifelse
  pm3 <- 
    t(sapply(1:dim(pm)[1], function(p, pm2=pm[p,], min2=bmin[p], max2=bmax[p])
       ifelse(pm2<min2, 2*min2 - pm2, 
       ifelse(pm2>max2, 2*max2 - pm2, 
              pm2) 
       )))

  # redo in case new reflected proposal value is out of bounds 
  .$dataf$pars[] <- 
    t(sapply(1:dim(pm)[1], function(p, pm2=pm3[p,], min2=bmin[p], max2=bmax[p])
       ifelse(pm2<min2, 2*min2 - pm2, 
       ifelse(pm2>max2, 2*max2 - pm2, 
              pm2) 
       )))
}

# restrict parameter proposals beyond the boundary by "folding" them back across the opposite boundary
boundary_handling_fold <- function(.) {

  pm   <- .$dataf$pars
  bmax <- .$mcmc$boundary_max 
  bmin <- .$mcmc$boundary_min

  # vectorised ifelse
  pm3 <- 
    t(sapply(1:dim(pm)[1], function(p, pm2=pm[p,], min2=bmin[p], max2=bmax[p])
       ifelse(pm2<min2, max2 - (min2-pm2), 
       ifelse(pm2>max2, min2 + (pm2-max2), 
              pm2) 
       )))

  # redo in case new "folded" proposal value is out of bounds 
  .$dataf$pars[] <- 
    t(sapply(1:dim(pm)[1], function(p, pm2=pm3[p,], min2=bmin[p], max2=bmax[p])
       ifelse(pm2<min2, max2 - (min2-pm2), 
       ifelse(pm2>max2, min2 + (pm2-max2), 
              pm2) 
       )))
}



# likelihood functions
################################

# expects model output to be probability
# - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {
  log(.$dataf$out)
}


# standard error probability density function with i.i.d. error residuals
f_proposal_lklihood_ssquared <- function(.) {

  obs_n                 <- length(.$dataf$obs)
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs
  SSR                   <- apply(error_residual_matrix, 2, function(v) sum(v^2) )

  # log-likelihood vector corresponding to each chain/column in .$dataf$pars matrix
  -(obs_n/2) * log(SSR)
}


# standard error probability density function with i.i.d. error residuals
# - incorporates measurement errors (unlike "ssquared" option)
f_proposal_lklihood_ssquared_se <- function(.) {

  # remove obs with zero measurement error
  sspos <- which(.$dataf$obsse>1e-9)
  obs_n <- length(sspos)

  # observed error (heteroscedastic and homoscedastic options)
  obsse <- if(.$wpars$mcmc$homosced) rep(mean(.$dataf$obsse[sspos]), obs_n)
           else                      .$dataf$obsse[sspos]

  # calculate error residual
  error_residual_matrix <- ( t(.$dataf$out)[sspos,]-.$dataf$obs[sspos] ) / obsse
  SSR                   <- apply(error_residual_matrix, 2, function(v) sum(v^2) )

  # log-likelihood vector corresponding to each chain/column in .$dataf$pars matrix
  -(obs_n/2) * log(2*pi) - sum(log(obsse)) - 0.5*SSR
}



# Proposal generation functions
################################

# generate proposal using DREAM algorithm (Vrugt et al. 2011)
proposal_generate_mcmc_dream <- function(.,j) {
  # - gamma = gamma in V2011
  # - lambda = e in V2011
  # - chain_pairs_n = delta in V2011
  # - length(jump_pars_ss) = dprime in V2011
  # - CR = m in V2011
  # - CR/nCR = CR in V2011
  # - p_CR = p_m in V2011

  # debugging
  #print(paste0('j = ',j))

  # initialise
  .$mcmc$jump[]          <- 0
  .$mcmc$current_state[] <- .$dataf$pars_array[,,j-1]

  # if adapting crossover values, compute standard deviation of each parameter/dimension
  if(.$wpars$mcmc$adapt_pCR) {
    .$mcmc$sd_state[] <- apply(.$mcmc$current_state, 1, sd )
    .$mcmc$sd_state[.$mcmc$sd_state==0] <- 1e-9
  }

  # chain pairs for jump on each chain
  chain_pairs_n <- sample(1:.$wpars$mcmc$chain_delta, .$wpars$mcmc$chains, T )

  # create proposals for each chain
  for(ii in 1:.$dataf$lp) {

    # determine chain pairs used to calculate each jump
    chain_pairs_ss <- t(sapply(1:chain_pairs_n[ii], function(v) sample((1:.$dataf$lp)[-ii],2,F) ))
    #chain_pairs_ss <- matrix((1:.$wpars$mcmc$chains)[-ii][1:(2*chain_pairs_n[ii])], chain_pairs_n[ii] )  

    # select crossover value
    # - weighted sample from multinomial distribution
    # - replacement relevant if this gets moved outside of chain loop 
    .$mcmc$CR[ii]  <- sample(1:.$wpars$mcmc$n_CR, 1, T, .$mcmc$p_CR )

    # determine which parameters will "crossover" (i.e. how many dimensions are sampled/updated jointly)
    zz             <- runif(.$mcmc$pars_n)
    jump_pars_ss   <- which(zz < (.$mcmc$CR[ii]/.$wpars$mcmc$n_CR) )
    if(length(jump_pars_ss)==0) jump_pars_ss <- which.min(zz)

    # jump rate / scaling factor
    # APW: for online code, when gamma = 1 CR=nCR i.e. all dimensions update
    gamma_d        <- 2.38 / sqrt(2*chain_pairs_n[ii]*length(jump_pars_ss))
    gamma          <- sample(c(gamma_d,1), 1, T, c(1-.$wpars$mcmc$p_gamma, .$wpars$mcmc$p_gamma ))
    if(gamma==1) { jump_pars_ss <- 1:.$mcmc$pars_n; .$mcmc$CR[ii] <- NA }

    # compute 'jump' for params to be updated/crossover (differential evolution)
    lambda         <- runif(length(jump_pars_ss), -.$wpars$mcmc$c_rand, .$wpars$mcmc$c_rand )
    chain_diff     <- apply(.$mcmc$current_state[jump_pars_ss,chain_pairs_ss[,1],drop=F], 1, sum ) - 
                      apply(.$mcmc$current_state[jump_pars_ss,chain_pairs_ss[,2],drop=F], 1, sum )
    .$mcmc$jump[jump_pars_ss,ii] <- .$wpars$mcmc$c_ergod*rnorm(length(jump_pars_ss)) + (1+lambda)*gamma*chain_diff 
    .$dataf$pars[,ii]            <- .$mcmc$current_state[,ii] + .$mcmc$jump[,ii]

  # chain loop
  }
}


# generate proposal using DREAM-ZS algorithm (Vrugt et al. 2016)
proposal_generate_mcmc_dreamzs <- function(.,j) {
  # - gamma = gamma in V2011
  # - lambda = e in V2011
  # - chain_pairs_n = delta in V2011
  # - length(jump_pars_ss) = dprime in V2011
  # - CR = m in V2011
  # - CR/nCR = CR in V2011
  # - p_CR = p_m in V2011

  # debugging
  #print(paste0('j = ',j))

  # initialise
  .$mcmc$jump[]          <- 0
  .$mcmc$current_state[] <- .$dataf$pars_array[,,j-1]

  # if adapting crossover values, compute standard deviation of each parameter/dimension
  if(.$wpars$mcmc$adapt_pCR) {
    .$mcmc$sd_state[] <- apply(.$mcmc$current_state, 1, sd )
    .$mcmc$sd_state[.$mcmc$sd_state==0] <- 1e-9
  }

  # determine past states used to calculate each jump
  past_states_n      <- 2*.$wpars$mcmc$chain_delta*.$wpars$mcmc$chains
  past_states_ss     <- sample(1:.$dataf$lps, past_states_n, F ) 
  past_states_sample <- .$dataf$past_states[,past_states_ss]
  chain_delta        <- sample(1:.$wpars$mcmc$chain_delta, .$wpars$mcmc$chains, T )
  #past_ss            <- numeric(4)

  # create proposals for each chain
  if(runif(1)>.$wpars$mcmc$psnooker) {
    for(ii in 1:.$dataf$lp) {

      # select crossover value
      # - weighted sample from multinomial distribution
      # - replacement relevant if this gets moved outside of chain loop 
      .$mcmc$CR[ii]  <- sample(1:.$wpars$mcmc$n_CR, 1, T, .$mcmc$p_CR )
  
      # determine which parameters will "crossover" (i.e. how many dimensions are sampled/updated jointly)
      zz             <- runif(.$mcmc$pars_n)
      jump_pars_ss   <- which(zz < (.$mcmc$CR[ii]/.$wpars$mcmc$n_CR) )
      if(length(jump_pars_ss)==0) jump_pars_ss <- which.min(zz)
  
      # subscripts for past_states_sample
      #past_ss[1]     <- past_ss[4] + 1
      #past_ss[2]     <- past_ss[1] + chain_delta[ii] - 1
      #past_ss[3]     <- past_ss[2] + 1
      #past_ss[4]     <- past_ss[3] + chain_delta[ii] - 1
      past_ss        <- t(sapply(1:chain_delta[ii], function(v) sample(1:past_states_n,2,F) ))
      #past_ss        <- t(sapply(1:chain_delta[ii], function(v) sample(1:.$dataf$lps,2,F) ))
  
      # jump rate / scaling factor
      gamma_d        <- 2.38 / sqrt(2*chain_delta[ii]*length(jump_pars_ss))
      gamma          <- sample(c(gamma_d,1), 1, T, c(1-.$wpars$mcmc$p_gamma, .$wpars$mcmc$p_gamma ))
      if(gamma==1) { jump_pars_ss <- 1:.$mcmc$pars_n; .$mcmc$CR[ii] <- NA }
  
      # compute 'jump' for params to be updated/crossover (differential evolution)
      lambda         <- runif(length(jump_pars_ss), -.$wpars$mcmc$c_rand, .$wpars$mcmc$c_rand )
      #chain_diff     <- apply(past_states_sample[jump_pars_ss, past_ss[1]:past_ss[2], drop=F], 1, sum ) - 
      #                  apply(past_states_sample[jump_pars_ss, past_ss[3]:past_ss[4], drop=F], 1, sum )
      chain_diff     <- apply(past_states_sample[jump_pars_ss, past_ss[,1], drop=F], 1, sum ) - 
                        apply(past_states_sample[jump_pars_ss, past_ss[,2], drop=F], 1, sum )
      .$mcmc$jump[jump_pars_ss,ii] <- .$wpars$mcmc$c_ergod*rnorm(length(jump_pars_ss)) + (1+lambda)*gamma*chain_diff 
      .$dataf$pars[,ii]            <- .$mcmc$current_state[,ii] + .$mcmc$jump[,ii]
    }
   
  # snooker update, orthogonal rather than parallel jump
  } else {

    draw <- matrix(1:past_states_n,2) 
    for(ii in 1:.$dataf$lp) {

      # select 3 past states 
      #past_ss           <- c(draw[,ii], sample((1:past_states_n)[-draw[,ii]], 1 ) )
      past_ss           <- sample(1:past_states_n, 3 ) 
      #past_ss           <- sample(1:.$dataf$lps, 3 ) 
  
      # jump rate / scaling factor
      gamma             <- 1.2 + runif(1) 
      .$mcmc$CR[ii]     <- NA
      
      # compute 'jump' for all parameters
      # - taken from https://github.com/Zaijab/DREAM/blob/master/functions/offde.m 
      Fv                <- t(.$mcmc$current_state[,ii] - past_states_sample[,past_ss[3]])
      D                 <- pmax(Fv %*% t(Fv), 1e-20 )
      chain_diff        <- Fv * as.numeric(sum((past_states_sample[,past_ss[1]] - past_states_sample[,past_ss[2]]) * Fv) / D )
      .$mcmc$jump[,ii]  <- gamma*chain_diff 
      .$dataf$pars[,ii] <- .$mcmc$current_state[,ii] + .$mcmc$jump[,ii]
    }
  }
}


# Adapt the probability of selecting a specific 'crossover' value
# - i.e. fraction of parameters (on average) that are updated for each proposal
adapt_pCR <- function(.) {

  .$mcmc$p_CR[] <- .$mcmc$j_true*.$wpars$mcmc$chains * (.$mcmc$jump_delta_norm/.$mcmc$CR_counter) / sum(.$mcmc$jump_delta_norm)
  .$mcmc$p_CR[] <- .$mcmc$p_CR/sum(.$mcmc$p_CR)

#   # debugging
#   print('');  print('adapt_pCR calculation')
#   print('.$mcmc$CR'); print(.$mcmc$CR)
#   print('.$mcmc$CR_counter'); print(.$mcmc$CR_counter)
#   print('.$mcmc$sd_state'); print(.$mcmc$sd_state)
#   print('.$mcmc$jump_delta_norm'); print(.$mcmc$jump_delta_norm)
#   print('')
#   print('.$mcmc$j_true'); print(.$mcmc$j_true)
#   print('.$mcmc$CR_counter'); print(.$mcmc$CR_counter)
#   print('.$mcmc$jump_delta_norm'); print(.$mcmc$jump_delta_norm)
#   print('.$mcmc$p_CR'); print(.$mcmc$p_CR)

  if(.$mcmc$j_true==.$wpars$mcmc$CR_burnin) {
    print('',quote=F); print('',quote=F)
    print(paste0('Adapted selection probabilities of crossover values, at iteration, ',.$mcmc$j_true,':'),quote=F); print(.$mcmc$p_CR,quote=F)
    .$wpars$mcmc$adapt_pCR[] <- .$mcmc$adapt_pCR[] <- F
  }
}



# proposal acceptance functions
#####################################

# Metropolis acceptance probability
# - assumes likelihood is log-likelihood
proposal_accept_mcmc_dream <- function(., prop_lklihood, curr_lklihood ) {
  alpha <- exp(prop_lklihood-curr_lklihood)
  alpha > runif(.$dataf$lp)
}
proposal_accept_mcmc_dreamzs <- proposal_accept_mcmc_dream 



# outlier detection functions
#####################################

# no outlier handling for Markov chains
mcmc_outlier_none <- function(.,j) {
  if(j==.$wpars$mcmc$maxiter) print('No option was chosen to identify outlier Markov chains.')
  numeric(0)
}


# function that detects outlier Markov chains using the Inter Quartile-Range (IQR) statistic (Vrugt et al 2011)
# - based on IQR across chains mean log posterior densities of last 50 % of burnin
# - all current likelihood functions return log-likelihood, so log likelihood not taken here
mcmc_outlier_iqr <- function(.,j) {

  .$dataf$omega[,.$mcmc$check_ss] <- apply(.$dataf$pars_lklihood[,.$mcmc$j_burnin50:j], 1, mean )
  q1q3     <- quantile(.$dataf$omega[1:.$wpars$mcmc$chains,.$mcmc$check_ss], prob=c(0.25,0.75), type=1 )
  iqr      <- q1q3[2]-q1q3[1]
  which(.$dataf$omega[,.$mcmc$check_ss] < (q1q3[1]-2*iqr))
}


# handle detected outliers
# - burnin restarted
mcmc_outlier_handling <- function(., outliers, j ) {
  
  print('',quote=F)
  print(paste('Outlier chain(s) detected. Chain(s):', outliers, 'at iteration:', .$mcmc$j_true ), quote=F )
  
  # replace outlier chain(s) & likelihood history for next iqr calculation
  replace_ss <- sample((1:.$wpars$mcmc$chains)[-outliers], length(outliers) )
  .$dataf$pars_array[1:.$mcmc$pars_n,outliers,j] <- .$dataf$pars_array[1:.$mcmc$pars_n,replace_ss,j]
  .$dataf$pars_lklihood[outliers,j]              <- .$dataf$pars_lklihood[replace_ss,j]
  
  # restart burn-in
  .$mcmc$outlier_detected <- T
  .$mcmc$j_start_burnin   <- j + 1
  .$mcmc$j_burnin50       <- j
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
    R_hat        <- sqrt( ((.$wpars$mcmc$chains+1)/.$wpars$mcmc$chains)*(sigma_hat/W) - ((iter_effective-2)/(.$wpars$mcmc$chains*iter_effective)) )
    R_hat        <- c(iter_effective, R_hat )

  } else {
    R_hat        <- c(0, rep(NA,.$mcmc$pars_n) )
  }

  R_hat
}


# final iteration handline
mcmc_handle_iter_final <- function(., R_hat ) {

  names(R_hat)[1] <- 'iterations since outlier detection'
  print('',quote=F); print('',quote=F)
  print(paste("At (final) iteration:", .$mcmc$j_true, ", R-statistic of Gelman and Rubin:"), quote=F )
  print(R_hat,quote=F); print('',quote=F)

  R_hat <- R_hat[2:length(R_hat)]
  if((sum(R_hat<1.2)) == length(R_hat)) print('ALL PARAMETERS CONVERGED.',quote=F)
  else                                  print('NON-CONVERGENCE, RESTART RECOMMENDED.',quote=F)
  print('',quote=F); print('',quote=F); print('',quote=F)
} 



# DEMC functions
# APW: not currrently in use as during development it became clear that this algorithm is not computationally parallel 
################################

#init_mcmc_demc <- function(.) NULL
#
## generate proposal using DE-MC algorithm
#proposal_generate_mcmc_demc <- function(., j ) {
#
#  # scaling factor
#  # APW: can be calculated once I think, fix
#  d          <- dim(.$dataf$pars)[1]
#  gamma_star <- 2.38 / sqrt(d + d)
#
#  # b-value should be small compared to width of target distribution; specifies range for drawn "randomization" value
#  b_rand  <- 0.01
#  uniform_r <- runif(1, min=(-b_rand), max=b_rand)
#
#  # evaluate for each chain
#  for(ii in 1:.$dataf$lp) {
#
#    # randomly select two different numbers R1 and R2 unequal to j, from uniform distribution without replacement
#    chain_pair <- sample((1:.$dataf$lp)[-.$dataf$lp], 2, F )
#
#    # evaluate for each parameter
#    for(jj in 1:d) {
#
#      # generate proposal via Differential Evolution
#      .$dataf$pars[jj,ii] <- .$dataf$pars_array[jj,ii,j-1] + uniform_r + 
#        gamma_star*( .$dataf$pars_array[jj,chain_pair[1],j-1] - .$dataf$pars_array[jj,chain_pair[2],j-1] )
#
#      # boundary handling 
#      .$boundary_handling(ii=ii, jj=jj )
#    }
#  }
#}
#
#
## calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
## - this is the same as DREAM
#proposal_accept_mcmc_demc <- function(., j, lklihood ) {
#
#  # Metropolis ratio
#  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[,j-1])
#  alpha        <- pmin(1, metrop_ratio)
#
#  # evaluate for each chain
#  for(ii in 1:.$dataf$lp) {
#    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
#    accept <- log(alpha[ii]) > log(runif(1, min = 0, max = 1))
#    .$dataf$pars_array[,ii,j]   <- if(accept)          .$dataf$pars[,ii] else .$dataf$pars_array[,ii,j-1]
#    .$dataf$pars_lklihood[ii,j] <- if(accept)          lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]
#    .$dataf$out_mcmc[ii,,j]     <- if(accept | j == 1) .$dataf$out[ii,]  else .$dataf$out_mcmc[ii,,(j-1)]
#  }
#}



### END ###
