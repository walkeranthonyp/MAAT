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

# linear decomp, Oleson 1963
f_decomp_lin        <- function(.,C,t,i) C[i]*.super$pars$k[[i]]                                              

# non-linear Michaelis-Menten decomp, a function of microbial biomass
f_decomp_MM_microbe <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[2]*C[i]) / (.super$pars$km[[i]]+C[i])

# reverse Michaelis-Menten decomp
f_decomp_RMM_microbe <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[2]*C[i]) / (.super$pars$km[[i]]+C[2])

# density-dependent turnover, Georgiou et al. 2017  
f_decomp_dd         <- function(.,C,t,i) (C[i]^.super$pars$beta) * .super$pars$k[[i]]                                  


# transfer functions
###################

# transfer all or nothing
f_transfer_all  <- function(.,C,t,from,to) 1
f_transfer_zero <- function(.,C,t,from,to) 0

# CUE or carbon transfer efficiency sets transfer from one pool to another
f_transfer_cue  <- function(.,C,t,from,to) .super$pars$cue[[from]]    

# CUE / transfer efficiency sets transfer subject to a maximum pool size 
# - can be used both for saturating MAOM pool and density dependent microbial growth efficiency
f_transfer_cue_sat <- function(.,C,t,from,to) .super$pars$cue[[from]] * (1-C[to]/.super$pars$poolmax[[to]])



### END ###
