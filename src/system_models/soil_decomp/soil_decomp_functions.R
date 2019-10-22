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
f_c1dec <- function(.,C,t) { (.super$pars$Af*.super$pars$ks*C[2] / (.super$pars$Km + C[2])) * C[1]  }
f_c2dec <- function(.,C,t) { (.super$pars$kb) * C[2] }

# transfer functions
f_eff12 <- function(.,C,t) {1 - .super$pars$r }
f_eff21 <- function(.,C,t) {1}



### MODIFIED SoilR FUNCTIONS
################################

# input rate matrix, single column, rows = length(cpools)
# - this is where inputs would be divided among pools
f_inputrates <- function(.,t) {
  matrix(ncol = 1, 
         c(.$env$litter, 0))
}

# decomp operator
#f_DotO  <- function(.,C,t) { 
#  matrix(
#    ncol=1, 
#    c(
#      .$c1dec(C,t), 
#      .$c2dec(C,t) 
#     ))
#}
f_DotO  <- function(.,C,t) { 
  dnames <- grep('decomposition.', names(.), value=T )
  id     <- sub('decomposition.d.', '', dnames)
  mat    <- matrix(ncol=1,nrow=length(dnames))
  for(i in id) mat[as.numeric(i),] <- .[[paste0('decomposition.d.',i)]](C=C,t=t)
  mat
}

# transfer matrix
f_transfermatrix <- function(., C, t ) {
  dnames <- grep('decomposition.', names(.), value=T )
  tnames <- grep('transition.', names(.), value=T )
  id     <- sub('transition.t.', '', tnames) 
  m      <- -1 * diag(nrow=length(dnames))
  for(i in id) m[matrix(rev(as.numeric(unlist(strsplit(i,'_to_')))),nrow=1)] <- .[[paste0('transition.t.',i)]](C=C,t=t)
  m
}

  
# function to solve
f_DotC <- function(., C, t) {
  .$transfermatrix(C, t) %*% .$DotO(C, t) + .$inputrates(t) 
}


# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$DotC(y, t)
  list(as.vector(YD))
}



### END ###
