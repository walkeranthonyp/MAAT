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
  print(m)
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
  print(unlist(.$state$outflux)) 
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t)
  #print(.$transfermatrix(y,t)); print(.$DotO(y,t)); print(.$input(t)) 
  #print(YD) 
  #print(unlist(.super$state$outflux)) 
  #print(list(as.vector(YD)))
  if(is.na(sum(as.vector(YD)))) stop('NAs in SoilR solver function.')
  list(as.vector(YD))
}



### END ###
