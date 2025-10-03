################################
#
# MAAT soil_decomp generalisation functions 
# - based entirely on functions in SoilR package functions modifed to work with proto objects and MAAT    
# 
# Carlos Sierra, Markus Mueller (SoilR developers) 
# AWalker, Matt Craig, October 2019 
#
################################



# MODIFIED SoilR FUNCTIONS
################################

# input matrix, single column, rows = cpools_n
# - this is where inputs would be divided among pools
# - APW I'm not sure this needs to be part of the solver
f_input <- function(., t ) {
  # why is .$env$litter not .super$env$litter??
  .$env$litter * matrix(unlist(.super$pars$input_coefs)[1:.super$pars$n_pools], ncol=1 )
}


f_input_clm5 <- function(., t ) {
  m <- matrix(0, nrow=.super$pars$n_pools, ncol=1 )
  f_litter_not_cwd <- 1 - .super$pars$input_coefs[[1]]
  m[1,] <- .super$pars$input_coefs[[1]]
  m[2,] <- f_litter_not_cwd * .super$pars$input_coefs[[2]]
  m[3,] <- f_litter_not_cwd * .super$pars$input_coefs[[3]]
  m[4,] <- f_litter_not_cwd * .super$pars$input_coefs[[4]]
  
  # print(.$env$litter * m)
  .$env$litter * m
}


## decomp matrix, single column, rows = n_pools
f_DotO <- function(., C, t ) {
  
  # search fns proto object for functions named starting with 'decomp.' and use to make dnames 
  #dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  
  # get integer id's of the decomp function for each pool
  # - could just use 1:n_pools for id but this reads directly
  #id     <- sub('decomp.d', '', dnames )
  
  # call functions and create decomp matrix 
  m      <- matrix(ncol=1, nrow=.super$pars$n_pools )
  #for(i in id) { 
  for(i in 1:nrow(m)) { 
    #print(i); print(.[[paste0('decomp.d',i)]]) 
    #m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C, t=t, i=as.numeric(i) )
    m[i,] <- .[[paste0('decomp.d',i)]](C=C, t=t, i=i )
  }
  m
}


f_outfluxes <- function(., C, t ) {
  # call outflux functions and assign to state outflux list 
  for(i in 1:.super$pars$n_pools) { 
    #print(i); print(.[[paste0('outflux.of',i)]]); print(.[[paste0('tcor.t',i)]]); print(.[[paste0('wcor.w',i)]]) 
    #ofv[i] <- .[[paste0('outflux.of',i)]](C=C, t=t, i=i ) * .[[paste0('tcor.t',i)]]() * .[[paste0('wcor.w',i)]]()
    
    for(of in unlist(.super$pars[[paste0('decomp_outflux',i)]]) ) {
      #print(i); print(of)
      .super$state$outflux[[paste0('of',of)]] <- 
        .[[paste0('outflux.of',of)]](C=C, t=t, i=i, of=of ) * .[[paste0('tcor.t',of)]](i=of) * .[[paste0('wcor.w',of)]](i=of) 
  }}
}


## decomp matrix with individual out flux terms calculated, single column, rows = n_pools
f_DotO_outflux <- function(., C, t ) {
 
  # APW: this prob needs to be called from the SoilR function before the main operation so that transfer functions can access the results 
  # call outflux functions and create outflux vector
  # APW: this needs to be in state list to allow access from decomp and transfer functions  
  for(i in 1:length(.super$pars$n_outfluxes)) { 
    #print(i); print(.[[paste0('outflux.of',i)]]); print(.[[paste0('tcor.t',i)]]); print(.[[paste0('wcor.w',i)]]) 
    #ofv[i] <- .[[paste0('outflux.of',i)]](C=C, t=t, i=i ) * .[[paste0('tcor.t',i)]]() * .[[paste0('wcor.w',i)]]() 
    .super$state$outflux[[paste0('of',i)]] <- .[[paste0('outflux.of',i)]](C=C, t=t, i=i ) * .[[paste0('tcor.t',i)]]() * .[[paste0('wcor.w',i)]]() 
  }


  # - the following two steps are not needed, think they're a hangover from pre n_pools days 
  # - could just use 1:n_pools for id but this reads directly from dnames
  # search fns proto object for functions named starting with 'decomp.' and use to make dnames 
  #dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  
  # get integer id's of the decomp function for each pool
  #id     <- sub('decomp.d', '', dnames )
  
  # call functions and create decomp matrix 
  m      <- matrix(ncol=1, nrow=.super$pars$n_pools )
  for(i in 1:nrow(m)) { 
    #print(i); print(.[[paste0('decomp.d',i)]]) 
    m[i,] <- .[[paste0('decomp.d',i)]](C=C, t=t, i=i )
  }
  m
}


## transfer matrix, square, n_pools extent
f_transfermatrix <- function(., C, t ) {
  
  # search fns proto object for functions named starting with 'transfer.' and use to make tnames 
  tnames <- grep('transfer\\.', names(.), value=T )
  
  # get integer id's of the transfer functions
  id     <- sub('transfer.t', '', tnames )
  
  # remove transfers that are not required, i.e. are 0 or from/to pools > n_pools
  id     <- id[!is.na(.super$fnames$transfer[paste0('t',id)])]
  
  # get integer id's of the from and to pools of the transfers 
  idm    <- apply(as.matrix(id,nrow=1), 1, function(i) as.numeric(unlist(strsplit(i,'_to_'))) )
  idm    <- apply(idm,1,function(v) v[v<=.super$pars$n_pools] )
  
  # call functions and create transfer matrix 
  m      <- -1 * diag(nrow=.super$pars$n_pools)
  #print(idm)
  #print(dim(idm)[1])
  #print(m)
  for(i in 1:dim(idm)[1]) {
    ss <- idm[i,]
    #print(i); print(ss); print(.[[paste0('transfer.t',ss[1],'_to_',ss[2])]]) 
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',ss[1],'_to_',ss[2])]](.=.,C=C, t=t, from=ss[1], to=ss[2] )
  }
  m
}
  
  
# deSolve/lsoda style functions to solve the SoilR style function
# - parms is a dummy argument to work with lsoda
f_solver_func_soilR <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t) 
  list(as.vector(YD))
}


# as above but allows calculation of individual fluxes out of pools prior to SoilR function execution
f_solver_func_soilR_outfluxes <- function(., t, y, parms) {
  print(y)
  .$outfluxes(y,t)
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t)
  #print(.$state$outflux); print(.$transfermatrix(y,t)); print(.$DotO(y,t)); print(.$input(t)) 
  #print(YD) 
  print(list(as.vector(YD)))
  list(as.vector(YD))
}



### END ###
