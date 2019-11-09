################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# Matt Craig, AWalker October 2019 
# Carlos Sierra, Markus Mueller (SoilR developers) 
#
################################



### FUNCTIONS
################################


# decay functions
###################

# linear decomp
f_decomp_lin        <- function(.,C,t,i) C[i]*.super$pars$k[[i]]                                              

# non-linear Michaelis-Menten decomp, a function of microbial biomass
f_decomp_MM_microbe <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[2]*C[i]) / (.super$pars$km[[i]]+C[i])

# density-dependent turnover as in Georgiou et al. 2017  
f_decomp_dd         <- function(.,C,t,i) (C[i]^.super$pars$beta) * .super$pars$k[[i]]                                  


# 3-pool model
#f_decomp_MM_microbe1 <- function(.,C,t) (.super$pars$vmax1*C[2]*C[1]) / (.super$pars$km1+C[1]) 
#f_decomp_MM_microbe3 <- function(.,C,t) (.super$pars$vmax3*C[2]*C[3]) / (.super$pars$km3+C[3])

# Microbial turnover hypotheses
# H1a) no density dependent microbial turnover                
#f_decomp_lin2 <- function(.,C,t) C[2]*.super$pars$k23                                              
                                                            
# H1b) density-dependent turnover as in Georgiou et al. 2017  
#f_decomp_dd2  <- function(.,C,t) (C[2]^.super$pars$beta) * .super$pars$k23                                  


                                                           
# transfer functions
###################

# transfer all or nothing
f_transfer_all  <- function(.,C,t,from,to) 1
f_transfer_zero <- function(.,C,t,from,to) 0

# CUE or carbon transfer efficiency sets transfer from one pool to another
f_transfer_cue  <- function(.,C,t,from,to) .super$pars$cue[[from]]    

# CUE / transfer effiiciency sets transfer subject to a maximum pool size 
# - can be used both for saturating MAOM pool and density dependent microbial growth efficiency
f_transfer_cue_sat <- function(.,C,t,from,to) .super$pars$cue[[from]] * (1-C[to]/.super$pars$poolmax[[to]])                  

# Microbial growth hypotheses (eff13 should also change)
# H2a) density dependent CUE                          
#f_transfer_cue_dd12 <- function(.,C,t) .super$pars$cuec1*(1-(C[2]/.super$pars$mbcmax))                  
                                                    
# H2b) no density dependent CUE (i.e. constant CUE)   
#f_transfer_cue12 <- function(.,C,t) .super$pars$cuec1                                         
                                                    
# MAOM saturation hypotheses
# H3a) no saturation of MAOM pool                     
#f_transfer_cue23 <- function(.,C,t) .super$pars$h                                             
                                                    
# H3b saturation of MAOM pool                         
#f_transfer_sat23 <- function(.,C,t) .super$pars$h * (1-C[3]/.super$pars$maommax)                     
                                                    
# Microbial growth hypotheses (eff12 should also change)
# H2a) density dependent CUE                          
#f_transfer_cue_dd32 <- function(.,C,t) .super$pars$cuec3*(1-(C[2]/.super$pars$mbcmax))                  
                                                    
# H2b) no density dependent CUE (i.e. constant CUE)   
#f_transfer_cue32 <- function(.,C,t) .super$pars$cuec3                                         



### END ###
