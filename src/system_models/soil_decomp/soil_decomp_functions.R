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
###################

f_decomp_MM_microbe <- function(.,C,t) { (.super$pars$Af*.super$pars$ks*C[2] / (.super$pars$Km + C[2])) * C[1]  }
f_decomp_lin        <- function(.,C,t) { (.super$pars$kb) * C[2] }

# 3-pool model
f_decomp_MM_microbe1 <- function(.,C,t) (.super$pars$vmax1*C[2]*C[1]) / (.super$pars$km1+C[1]) 
f_decomp_MM_microbe3 <- function(.,C,t) (.super$pars$vmax3*C[2]*C[3]) / (.super$pars$km3+C[3])

# Microbial turnover hypotheses
# H1a) no density dependent microbial turnover                
f_decomp_lin2 <- function(.,C,t) C[2]*.super$pars$k23                                              
                                                            
# H1b) density-dependent turnover as in Georgiou et al. 2017  
f_decomp_dd2  <- function(.,C,t) (C[2]^.super$pars$beta) * .super$pars$k23                                  


                                                           
# transfer functions
###################

f_transfer_resploss <- function(.,C,t) {1 - .super$pars$r }
f_transfer_all      <- function(.,C,t) {1}


f_transfer_zero21 <- function(.,C,t) 0
f_transfer_zero13 <- function(.,C,t) 0
f_transfer_zero31 <- function(.,C,t) 0

# Microbial growth hypotheses (eff13 should also change)
# H2a) density dependent CUE                          
f_transfer_cue_dd12 <- function(.,C,t) .super$pars$cuec1*(1-(C[2]/.super$pars$mbcmax))                  
                                                    
# H2b) no density dependent CUE (i.e. constant CUE)   
f_transfer_cue12 <- function(.,C,t) .super$pars$cuec1                                         
                                                    

# MAOM saturation hypotheses
# H3a) no saturation of MAOM pool                     
f_transfer_cue23 <- function(.,C,t) .super$pars$h                                             
                                                    
# H3b saturation of MAOM pool                         
f_transfer_sat23 <- function(.,C,t) .super$pars$h * (1-C[3]/.super$pars$maommax)                     
                                                    

# Microbial growth hypotheses (eff12 should also change)
# H2a) density dependent CUE                          
f_transfer_cue_dd32 <- function(.,C,t) .super$pars$cuec3*(1-(C[2]/.super$pars$mbcmax))                  
                                                    
# H2b) no density dependent CUE (i.e. constant CUE)   
f_transfer_cue32 <- function(.,C,t) .super$pars$cuec3                                         



# MODIFIED SoilR FUNCTIONS
################################

# input rate matrix, single column, rows = cpools_n
# - this is where inputs would be divided among pools
f_inputrates <- function(.,t) {
  matrix(ncol = 1, c(.$env$litter, rep(0,.super$state$cpools_n-1)) )
}


# decomp vector, single column, rows = cpools_n
f_DotO  <- function(.,C,t) { 
  dnames <- grep('decomp\\.', names(.), value=T )
  id     <- sub('decomp.d', '', dnames)
  m      <- matrix(ncol=1, nrow=.super$state$cpools_n )
  for(i in id) m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t)
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

  
## function to solve
#f_DotC <- function(., C, t) {
#  .$transfermatrix(C,t) %*% .$DotO(C,t) + .$inputrates(t) 
#}


# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$inputrates(t) 
  #YD = .$DotC(y,t)
  list(as.vector(YD))
}



### END ###
