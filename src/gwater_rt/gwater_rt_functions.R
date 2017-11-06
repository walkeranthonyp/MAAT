################################
#
# gwater_rt for MAAT object functions
# 
# AWalker March 2017
#
################################


f_none <- function(.) {
  NA
}

# robust numerical solution to the quadratic
# taken from Numerical Recipes
quad_sol <- function(a,b,c,out='lower') {

  q     <- -0.5 * ( b + sign(b)*(b^2 - 4*a*c)^0.5 )
  roots <- c( q/a , c/q )
  
  if(out=='lower')      min(roots,na.rm=T) 
  else if(out=='upper') max(roots,na.rm=T)
  else roots 
}

# function that mimics some of the behaviour of the matlab 'diag' function
diag_m <- function(v,k=0) {
  if(k==0) diag(v)
  else if(k>0) {
    rbind(cbind(matrix(0,length(v),k),diag(v)), matrix(0,k,length(v)+k) )
  } else {
    k <- abs(k)
    rbind(matrix(0,k,length(v)+k), cbind(diag(v), matrix(0,length(v),k) ))
  }  
}



### MODEL FUNCTIONS
################################

# recharge
# power law
f_rechrg_power <- function(.) {
  (.$pars$a*(.$env$precip - 355.6)^0.5) * 1e-3/365
}

# linear
f_rechrg_lin <- function(.) {
  (.$pars$b*(.$env$precip - 399.8))     * 1e-3/365
}


# transport over geological domain

# single layer model
f_trans_single <- function(.) {
  vapply(1:.$pars$nx, function(jj,.) sqrt(.$env$h1^2 - (.$env$h1^2-.$env$h2^2)*.$state_pars$x[jj]/.$pars$L + 
                                            .$state$recharge*(.$pars$L-.$state_pars$x[jj])*.$state_pars$x[jj]/.$pars$K), 1, .=. )
}

# double layer model
f_trans_double <- function(.) {
  h11 <-  Double_layer_model_h(.)
  h11[seq(1,101,5)]
}

Double_layer_model_h <- function(.) {
  # Domain information
  nx    <- 101 # for some reason the horizontal discretisation is finer for this model
  nx1   <- 70
  h     <- numeric(nx)
  A     <- B  <- x2 <- numeric(nx-2)
  x1    <- x3 <- numeric(nx-3)
  delta <- 100
  
  # calculations
  B[1]        <- -2*delta^2*.$state$recharge - .$pars$K1*.$env$h1^2
  B[nx-2]     <- -2*delta^2*.$state$recharge - .$pars$K2*.$env$h2^2
  B[2:(nx-3)] <- -2*delta^2*.$state$recharge
  
  x1[1:(nx1-2)]      <- .$pars$K1
  x1[(nx1-1):(nx-3)] <- .$pars$K2
  x2[1:(nx1-2)]      <- -2*.$pars$K1
  x2[nx1:(nx-2)]     <- -2*.$pars$K2
  x2[nx1-1]          <- -(.$pars$K1 + .$pars$K2)
  x3[1:(nx1-2)]      <- .$pars$K1
  x3[(nx1-1):(nx-3)] <- .$pars$K2
  
  A                  <- diag_m(x1,1) + diag_m(x2) + diag_m(x3,-1)
  ht                 <- B %*% solve(A)
  htt                <- ht^0.5
  h[1]               <- .$env$h1
  h[nx]              <- .$env$h2
  h[2:(nx-1)]        <- htt
  h
}
















