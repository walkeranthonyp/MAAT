################################
#
# gwater_rt for MAAT object functions
# 
# AWalker March 2017
#
################################

source('../generic_mathematical_functions.R')

f_none <- function(.) {
  NA
}



### MODEL FUNCTIONS
################################

# recharge
# power law (R1 in Dai et al 2017)
f_rechrg_power <- function(.) {
  (5.04*.super$pars$a*(.super$env$precip - 355.6)^0.5) * 1e-3 / 365
}

# linear (R2 in Dai et al 2017)
f_rechrg_lin <- function(.) {
  (.super$pars$b*(.super$env$precip - 399.8)) * 1e-3 / 365
}


# transport over geological domain
# single layer model
f_geol_single <- function(.) {
  vapply(1:.super$pars$nx, function(jj,.) sqrt(.super$env$h1^2 - (.super$env$h1^2-.super$env$h2^2)*.super$state_pars$x[jj]/.super$pars$L + 
                                            .super$state$recharge*(.super$pars$L-.super$state_pars$x[jj])*.super$state_pars$x[jj]/.super$pars$K), 1, .=. )
}

# double layer model
f_geol_double <- function(.) {
  h11 <-  f_double_layer_model_h(.)
  h11[seq(1,101,5)]
}

f_double_layer_model_h <- function(.) {
  # Domain information
  nx    <- 101 # for some reason the horizontal discretisation is finer for this model
  nx1   <- 70
  h     <- numeric(nx)
  A     <- B  <- x2 <- numeric(nx-2)
  x1    <- x3 <- numeric(nx-3)
  delta <- 100
  
  # calculations
  B[1]        <- -2*delta^2*.super$state$recharge - .super$pars$K1*.super$env$h1^2
  B[nx-2]     <- -2*delta^2*.super$state$recharge - .super$pars$K2*.super$env$h2^2
  B[2:(nx-3)] <- -2*delta^2*.super$state$recharge
  
  x1[1:(nx1-2)]      <- .super$pars$K1
  x1[(nx1-1):(nx-3)] <- .super$pars$K2
  x2[1:(nx1-2)]      <- -2*.super$pars$K1
  x2[nx1:(nx-2)]     <- -2*.super$pars$K2
  x2[nx1-1]          <- -(.super$pars$K1 + .super$pars$K2)
  x3[1:(nx1-2)]      <- .super$pars$K1
  x3[(nx1-1):(nx-3)] <- .super$pars$K2
  
  A                  <- diag_m(x1,1) + diag_m(x2) + diag_m(x3,-1)
  ht                 <- B %*% solve(A)
  htt                <- ht^0.5
  h[1]               <- .super$env$h1
  h[nx]              <- .super$env$h2
  h[2:(nx-1)]        <- htt
  h
}



### END ###
