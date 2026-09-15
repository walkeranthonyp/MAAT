################################
#
# MAAT soil_decomp generalisation functions 
# - based entirely on functions in SoilR package functions modifed to work with proto objects and MAAT    
# 
# Carlos Sierra, Markus Mueller (SoilR developers) 
# AWalker, Matt Craig, October 2019 onwards; Hannah DeHetre, August 2026  
#
################################



# MODIFIED SoilR FUNCTIONS
################################

# decomp matrix, single column, rows = n_pools
f_DotO <- function(., C, t ) {
  
  # search fns proto object for functions named starting with 'decomp.' and use to make dnames 
  #dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  
  # get integer id's of the decomp function for each pool
  # - could just use 1:n_pools for id but this reads directly
  #id     <- sub('decomp.d', '', dnames )
  
  # call functions and create decomp matrix 
  m <- matrix(ncol=1, nrow=.super$pars$n_pools )
  #print(m)
  for(i in 1:nrow(m)) { 
    #print(i); print(.[[paste0('decomp.d',i)]]) 
    m[i,] <- .[[paste0('decomp.d',i)]](C=C, t=t, i=i )
    #print(m)
  }
  m
}


# function that calls outflux functions multiplied by env scalars 
# - assigns output to state outflux list 
f_outfluxes <- function(., C, t ) {

  # cycle through pools to assign correct pool to outflux 
  for(i in 1:.super$pars$n_pools) { 
    #print(i); print(.[[paste0('outflux.of',i)]]); print(.[[paste0('tcor.t',i)]]); print(.[[paste0('wcor.w',i)]]) 

    # cycle through outfluxes associated with each pool 
    for(of in unlist(.super$pars[[paste0('decomp_outflux',i)]]) ) {
      #print('f_outfluxes, i, of:')
      #print(c(i,of)) 
      .super$state$outflux[[paste0('of',of)]] <- 
        .[[paste0('outflux.of',of)]](C=C, t=t, i=i, of=of ) * 
        .[[paste0('tcor.tcor',of)]](i=of) * .[[paste0('wcor.wcor',of)]](i=of) * .[[paste0('scor.scor',of)]](i=of) 
  }}
}


# transfer matrix, square, n_pools extent
f_transfermatrix <- function(., C, t ) {
  
  #print('')
  #print('Transfer matrix generation:')
  #print('**********************')
  
  # search fns proto object for functions named starting with 'transfer.' and use to make tnames 
  tnames <- grep('transfer\\.', names(.), value=T )
  #print(as.matrix(tnames)) 
 
  # get integer id's of the transfer functions
  id     <- sub('transfer.t', '', tnames )

  # remove transfers that are not required, i.e. are 0 or from/to pools > n_pools
  id     <- id[!is.na(.super$fnames$transfer[paste0('t',id)])]
  #print(as.matrix(id))  
  
  # get integer id's of the from and to pools of the transfers 
  idm    <- apply(as.matrix(id,nrow=1), 1, function(i) as.numeric(unlist(strsplit(i,'_to_'))) )
  idm    <- apply(idm,1,function(v) v[v<=.super$pars$n_pools] )
  
  # call functions and create transfer matrix 
  m      <- -1 * diag(nrow=.super$pars$n_pools)
  #print(idm)
  #print(dim(idm)[1])
  #print('')
  for(i in 1:dim(idm)[1]) {
    ss <- idm[i,]
    #print(''); print(paste('transfer', i, '... from pool',ss[1],'to',ss[2],':')) 
    #print(.super$fnames$transfer[[paste0('t',ss[1],'_to_',ss[2])]])
    #print(.[[paste0('transfer.t',ss[1],'_to_',ss[2])]])
    #print(.[[paste0('transfer.t',ss[1],'_to_',ss[2])]](.=.,C=C, t=t, from=ss[1], to=ss[2] )) 
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',ss[1],'_to_',ss[2])]](.=., C=C, t=t, from=ss[1], to=ss[2] )
  }
  #print(m)
  m
}
  
  
# deSolve/lsoda style functions to solve the SoilR style function
# - parms is a dummy argument to work with lsoda
f_solver_func_soilR <- function(., t, y, parms) {
  #print(y)
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t) 
  #print(list(as.vector(YD)))
  list(as.vector(YD))
}


# as above but calculates individual fluxes out of pools prior to SoilR function execution
f_solver_func_soilR_outfluxes <- function(., t, y, parms) {
  #print(y)
  .$outfluxes(y,t)
  #print(unlist(.$state$outflux)) 
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t)
  #print(.$transfermatrix(y,t)); print(.$DotO(y,t)); print(.$input(t)) 
  #print(YD) 
  #print(unlist(.super$state$outflux)) 
  #print(list(as.vector(YD)))
  if(is.na(sum(as.vector(YD)))) stop('NAs in SoilR solver function.')
  list(as.vector(YD))
}


# as above but calculates respiration and other loss fluxes
# - by temporarily assigning them a pool to accumulate the loss fluxes 
f_solver_func_soilR_outfluxes_lossfluxes <- function(., t, y, parms ) {
  
  # calculate outfluxes
  .$outfluxes(y,t)
 
  # calculate soilR terms
  tm    <- .$transfermatrix(y,t)
  do    <- .$DotO(y,t) 
  input <- .$input(t)
  if(.super$cpars$verbose) {
    print('')
    print('Dynamic solver function')
    print('initial state:'); print(y)
    print('outfluxes:'); print(unlist(.super$state$outflux)) 
    print('initial soilR terms:'); print(input); print(do); print(tm)
  }

  # incorporate loss fluxes in soilR terms
  if(!parms$steadystate) { 
    n_loss_fluxes <- .super$state_pars$n_loss_fluxes 
    lfm           <- .super$state_pars$lfm 
    # - simple add zeros for do and input 
    do    <- rbind(do, lfm )
    input <- rbind(input, lfm ) 
    # - add rows and cols to transfer matrix
    tm_add_cols    <- matrix(0, nrow=.super$pars$n_pools+n_loss_fluxes , ncol=n_loss_fluxes ) 
    pool_loss_frac <- abs(apply(tm, 2, sum ))
    tm_add_rows    <- matrix(pool_loss_frac, nrow=1 )
    if(n_loss_fluxes!=1) {
      tm_add_rows0   <- tm_add_rows 
      tm_add_rows0[] <- 0 
      for(i in 1:.super$pars$n_pools) { 
        flux_from_pool_ids <- unlist(.super$pars[[paste0('decomp_outflux',i)]])
        for(of in flux_from_pool_ids) { 
          if(of==.super$pars$leaching_outflux) {
            target_flux_prop <- 
              if(length(flux_from_pool_ids)==1) 1
              else {
                fluxes_from_pool <- unlist(.super$state$outflux)[flux_from_pool_ids] 
                fluxes_from_pool[which(flux_from_pool_ids==of)] / sum(fluxes_from_pool) 
              } 
            #print(fluxes_from_pool)
            #print(target_flux_prop)
            tm_add_rows      <- rbind(tm_add_rows,tm_add_rows0) 
            tm_add_rows[1,i] <- tm_add_rows[1,i] - target_flux_prop 
            tm_add_rows[dim(tm_add_rows)[1],i] <- target_flux_prop 
          }  
      }}
    } 
    tm_add_rows[tm_add_rows<1e-15] <- 0 
    tm <- rbind(tm,tm_add_rows )
    tm <- cbind(tm,tm_add_cols )
    #print('Augmented soilR terms:')
    #print(tm); print(do); print(input) 
  }
 
  # cacluate soilR function 
  YD = tm %*% do + input 
  #print(YD) 
  #print(list(as.vector(YD)))
  if(is.na(sum(as.vector(YD)))) stop('NAs in SoilR solver function.')
  list(as.vector(YD))
}



# calculates instantaneous loss fluxes after steadystate solver has run
# - good for steady-state mass balance check but not dynamic which needs integrated loss fluxes 
f_calc_loss_fluxes <- function(.) {

  # calculate loss fractions for each pool from transfer matrix
  # - 1,1 arguments to tm and ,1 doto are dummy args
  .super$state_pars$transfer_matrix    <- .$transfermatrix(1,1)
  .super$state_pars$total_fluxes       <- as.numeric(.$DotO(.super$state$cpools[,1], 1 ))
  .super$state_pars$flux_matrix        <- t(.super$state_pars$total_fluxes*t(.super$state_pars$transfer_matrix))
  pool_loss_frac                       <- abs(apply(.super$state_pars$transfer_matrix, 2, sum ))
  pool_loss_frac[pool_loss_frac<1e-15] <- 0 

  # sum pool fluxes multiplied by loss fractions
  pool_loss_total <- pool_loss_frac * .super$state_pars$total_fluxes

  # assign losses
  .super$state$leaching    <- 0
  .super$state$respiration <- sum(pool_loss_total) 
  if(!is.na(.super$pars$leaching_outflux)) {
    .super$state$leaching    <- .super$state$outflux[[.super$pars$leaching_outflux]]
    .super$state$respiration <- .super$state$respiration - .super$state$leaching 
  } 

  # pool imbalance
  # - at steady state pool inputs and output should sum to zero
  .super$state_pars$pool_balance <- apply(cbind(.$input(),.super$state_pars$flux_matrix), 1, sum )

  if(.super$cpars$verbose) {
    print('') 
    print('Calculate loss fluxes:') 
    print('transfer matrix:')
    print(.super$state_pars$transfer_matrix)
    print('pool loss frac.:')
    print(pool_loss_frac) 
    print('pool loss abs.:')
    print(pool_loss_total)
    print('respiration:')
    print(.super$state$respiration) 
    print('leaching:')
    print(.super$state$leaching) 
  }
}


f_mass_balance <- function(., steadystate=T ) {

  input    <- sum(.$input()) 
  output   <- .super$state$respiration + .super$state$leaching 
  previous <- sum(.super$state_pars$previous_state)
  current  <- sum(.super$state$cpools)

  # at steady state inputs = outputs
  mass_imbalance <- .super$state_pars$mass_imbalance <-  
    if(steadystate) input - output
    else            input - output - (current - previous)
  # error tolerance determined as proportion of input sum
  error <- abs(mass_imbalance)/input > .super$pars$error_tolerance

  if(!is.na(error)) {
    if(error|.super$cpars$verbose) {
      print('') 
      print('Calculate mass balance:') 
      print('input:') 
      print(input) 
      print('output:') 
      print(output) 
      print('current state sum:') 
      print(current) 
      print('previous state sum:') 
      print(previous) 
      print('delta state sum:') 
      print(current-previous) 
      print('mass imbalance (absolute):') 
      print(mass_imbalance) 
      print('mass imbalance (proportion of input):') 
      print(mass_imbalance/input) 
    }
    if(error) stop(paste('ERROR:: mass imbalance (should be zero) = ', mass_imbalance )) 
  }
}


### END ###
