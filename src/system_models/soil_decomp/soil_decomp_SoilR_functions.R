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
# - APW I don't think this needs to be part of the solver
f_input <- function(., t ) {
  #why is .$env$litter not .super$env$litter??
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
# - APW seems a little convoluted, why can't just use 1:n_pools for the id
f_DotO <- function(., C, t ) {
  
  # search fns proto object for functions named starting with 'decomp.' and use to make dnames 
  dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  
  # get integer id's of the decomp function for each pool
  id     <- sub('decomp.d', '', dnames )
  
  # call functions and create decomp matrix 
  m      <- matrix(ncol=1, nrow=.super$pars$n_pools )
  for(i in id) { 
    #print(i); print(.[[paste0('decomp.d',i)]]) 
    m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C, t=t, i=as.numeric(i) )
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
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',ss[1],'_to_',ss[2])]](C=C, t=t, from=ss[1], to=ss[2] )
  }
  m
}
  
  
# deSolve/lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func_soilR <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t) 
  list(as.vector(YD))
}



### END ###