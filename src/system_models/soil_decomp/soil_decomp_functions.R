################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# Matt Craig, AWalker October 2019 
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



### END ###
