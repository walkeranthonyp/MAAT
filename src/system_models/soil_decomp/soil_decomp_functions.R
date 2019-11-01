################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# AWalker, Matt Craig, October 2019 
# Carlos Sierra, Matthais Mueller (SoilR developers) 
#
################################


### FUNCTIONS
################################

# decay functions
f_decomp_MM_microbe <- function(.,C,t) { (.super$pars$Af*.super$pars$ks*C[2] / (.super$pars$Km + C[2])) * C[1]  }
f_decomp_lin        <- function(.,C,t) { (.super$pars$kb) * C[2] }

# transfer functions
f_transfer_resploss <- function(.,C,t) {1 - .super$pars$r }
f_transfer_all      <- function(.,C,t) {1}



### MODIFIED SoilR FUNCTIONS
################################

# input rate matrix, single column, rows = length(cpools)
# - this is where inputs would be divided among pools
f_inputrates <- function(.,t) {
  matrix(ncol = 1, 
         c(.$env$litter, 0))
}


# decomp vector 
f_DotO  <- function(.,C,t) { 
  dnames <- grep('decomp\\.', names(.), value=T )
  id     <- sub('decomp.d', '', dnames)
  mat    <- matrix(ncol=1,nrow=length(dnames))
  for(i in id) mat[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t)
  mat
}


# transfer matrix
f_transfermatrix <- function(., C, t ) {
  dnames <- grep('decomp\\.', names(.), value=T )
  tnames <- grep('transfer\\.', names(.), value=T )
  id     <- sub('transfer.t', '', tnames)
  m      <- -1 * diag(nrow=length(dnames))
  for(i in id) m[matrix(rev(as.numeric(unlist(strsplit(i,'_to_')))),nrow=1)] <- .[[paste0('transfer.t',i)]](C=C,t=t)
  m
}

  
# function to solve
f_DotC <- function(., C, t) {
  .$transfermatrix(C,t) %*% .$DotO(C,t) + .$inputrates(t) 
}


# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$DotC(y,t)
  list(as.vector(YD))
}



### END ###
