################################
#
# MAAT soil_decomp generalisation functions 
# - based entirely on functions in SoilR package functions modifed to work with proto objects and MAAT    
# 
# Carlos Sierra, Matthais Mueller (SoilR developers) 
# AWalker, Matt Craig, October 2019 
#
################################



# MODIFIED SoilR FUNCTIONS
################################

# input rate matrix, single column, rows = cpools_n
# - this is where inputs would be divided among pools
f_inputrates <- function(., t ) {
  matrix(ncol = 1, c(.$env$litter, rep(0,.super$state$cpools_n-1)) )
}


# decomp vector, single column, rows = cpools_n
f_DotO <- function(., C, t ) { 
  dnames <- grep('decomp\\.', names(.), value=T )
  id     <- sub('decomp.d', '', dnames)
  m      <- matrix(ncol=1, nrow=.super$state$cpools_n )
  for(i in id) m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t)
  m
}

# decomp vector, single column, rows = cpools_n
f_DotOi <- function(., C, t ) { 
  dnames <- grep('decomp\\.', names(.), value=T )
  id     <- sub('decomp.d', '', dnames)
  m      <- matrix(ncol=1, nrow=.super$state$cpools_n )
  for(i in id) m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t,i=as.numeric(i))
  m
}


# transfer matrix, square, cpools_n extent
f_transfermatrix <- function(., C, t ) {
  tnames <- grep('transfer\\.', names(.), value=T )
  id     <- sub('transfer.t', '', tnames)
  #print(id)
  m      <- -1 * diag(nrow=.super$state$cpools_n)
  #print(m)
  #print(.super$state)
  for(i in id) m[matrix(rev(as.numeric(unlist(strsplit(i,'_to_')))),nrow=1)] <- .[[paste0('transfer.t',i)]](C=C,t=t)
  m
}

# transfer matrix, square, cpools_n extent
f_transfermatrixi <- function(., C, t ) {
  tnames <- grep('transfer\\.', names(.), value=T )
  id     <- sub('transfer.t', '', tnames)
  #print(id)
  m      <- -1 * diag(nrow=.super$state$cpools_n)
  #print(m)
  #print(.super$state)
  for(i in id) {
    ss <- as.numeric(unlist(strsplit(i,'_to_'))) 
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',i)]](C=C,t=t,from=ss[1],to=ss[2])
  }
  m
}
  
  
# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$inputrates(t) 
  list(as.vector(YD))
}



### END ###
