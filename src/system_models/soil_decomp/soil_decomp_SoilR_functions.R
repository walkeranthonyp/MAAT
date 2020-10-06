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
f_input <- function(., t ) {
  
  .$env$litter * matrix(unlist(.super$pars$input_coefs)[1:.super$pars$n_pools], ncol=1 )
}


## decomp matrix, single column, rows = n_pools
f_DotO <- function(., C, t ) { 
  dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  id     <- sub('decomp.d', '', dnames )
  m      <- matrix(ncol=1, nrow=.super$pars$n_pools )
  # call functions and create decomp matrix 
  for(i in id) { m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t,i=as.numeric(i))
    #print(id); print(.[[paste0('decomp.d',i)]]) 
  }
  m
}


## transfer matrix, square, n_pools extent
f_transfermatrix <- function(., C, t ) {
  tnames <- grep('transfer\\.', names(.), value=T )
  id     <- sub('transfer.t', '', tnames )
  # remove transfers that are not required, i.e. are 0 or from/to pools > n_pools
  id     <- id[!is.na(.super$fnames$transfer[paste0('t',id)])]
  idm    <- apply(as.matrix(id,nrow=1), 1, function(i) as.numeric(unlist(strsplit(i,'_to_'))) )
  idm    <- apply(idm,1,function(v) v[v<=.super$pars$n_pools] )
  # call functions and create transfer matrix 
  m      <- -1 * diag(nrow=.super$pars$n_pools)
  #print(idm)
  #print(dim(idm)[1])
  #print(m)
  for(i in 1:dim(idm)[1]) {
    ss <- idm[i,]
    #print(i) 
    #print(ss)
    #print(.[[paste0('transfer.t',ss[1],'_to_',ss[2])]]) 
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',ss[1],'_to_',ss[2])]](C=C,t=t,from=ss[1],to=ss[2])
  }
  m
}
  
  
# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t) 
  list(as.vector(YD))
}

# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func_corpse <- function(., t, y, parms) {
  dCs <- .$input(t)[[1]] + .$desorp.ds5(t=t,C=y,i=5) - .$sorp.s1(t=t,C=y,i=1)*.$scor(.) - .$decomp.d1(t = t, C=y, i=1)*.$wcor(.)*.$tcor.t1(i=1)
  dCr <- .$input(t)[[2]] + .$desorp.ds6(t=t,C=y,i=6) - .$sorp.s2(t=t,C=y,i=2)*.$scor(.) - .$decomp.d2(t = t, C=y, i=2)*.$wcor(.)*.$tcor.t2(i=2)
  dCn <- .$decomp.d4(t = t, C=y, i=4)*.super$pars$cue[[4]] + .$desorp.ds7(t=t,C=y,i=7) - .$sorp.s3(t=t,C=y,i=3)*.$scor(.) - .$decomp.d3(t = t, C=y, i=3)*.$wcor(.)*.$tcor.t3(i=3)
  dM <- .$decomp.d1(t = t, C=y, i=1)*.$wcor(.)*.$tcor.t1(i=1)*.super$pars$cue[[1]]+.$decomp.d2(t = t, C=y, i=2)*.$wcor(.)*.$tcor.t2(i=2)*.super$pars$cue[[2]]+.$decomp.d3(t = t, C=y, i=3)*.$wcor(.)*.$tcor.t3(i=3)*.super$pars$cue[[3]] - .$decomp.d4(t = t, C=y, i=4)
  dPs <- .$sorp.s1(t=t,C=y,i=1)*.$scor(.) - .$desorp.ds5(t=t,C=y,i=5)
  dPr <- .$sorp.s2(t=t,C=y,i=2)*.$scor(.) - .$desorp.ds6(t=t,C=y,i=6)
  dPn <- .$sorp.s3(t=t,C=y,i=3)*.$scor(.) - .$desorp.ds7(t=t,C=y,i=7)
  list(c(dCs, dCr, dCn, dM, dPs, dPr, dPn))
}


### END ###
