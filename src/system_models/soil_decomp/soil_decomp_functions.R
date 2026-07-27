################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# Matt Craig, AWalker October 2019 
#
################################

source('soil_decomp_parameter_functions.R')
source('soil_decomp_temperature_functions.R')
source('soil_decomp_water_functions.R')



# FUNCTIONS
################################

# litter input functions 
################################

# divides inputs among pools
# returns input matrix, single column, rows = n_pools

f_input <- function(., t ) {
  # APW: why is .$env$litter not .super$env$litter??
  .$env$litter * matrix(unlist(.super$pars$input_coef)[1:.super$pars$n_pools], ncol=1 )
}

# APW: need to make CWD litter input?
f_input_clm5 <- function(., t ) {
  m <- matrix(0, nrow=.super$pars$n_pools, ncol=1 )
  f_litter_not_cwd <- 1 - .super$pars$input_coef[[1]]
  m[1,] <- .super$pars$input_coef[[1]]
  m[2,] <- f_litter_not_cwd * .super$pars$input_coef[[2]]
  m[3,] <- f_litter_not_cwd * .super$pars$input_coef[[3]]
  m[4,] <- f_litter_not_cwd * .super$pars$input_coef[[4]]
  
  # print(.$env$litter * m)
  .$env$litter * m
}


# function for converting ANPP to litter inputs from MIMICS
# APW: some of these paramters could be calculated in calc_pars
# APW: I'd like to make this calculate litter and put the in the env list, and then do the partitioning as per the other models 
# APW: input a function of depth
f_input_mimics <- function(., t ) {
  # convert anpp to litter inputs
  EST_LIT_in <- .super$env$anpp / (365*24) # gC/m2/h (from gC/m2/y)
  EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h)
  input      <- .super$env$litter <- EST_LIT/.super$env$depth 
 
  # partitioning inputs
  i_LITm <- input * .super$state_pars$quality_mod     * (1-.super$pars$input_coef[[5]])
  i_LITs <- input * (1-.super$state_pars$quality_mod) * (1-.super$pars$input_coef[[6]])
  i_SOMp <- input * .super$state_pars$quality_mod     * .super$pars$input_coef[[5]]
  i_SOMc <- input * (1-.super$state_pars$quality_mod) * .super$pars$input_coef[[6]]

  # input matrix
  matrix(c(i_LITm, i_LITs, 0, 0, i_SOMp, i_SOMc, 0 ), ncol=1 ) 
}



# decomp functions
# APW: we might want a broader class of functions called flux functions given we have uptake, sorption, and other processes 
#      that cause a flux from one pool to another 
################################

# APW: not sure this is needed 
f_decomp_none <- function(., ... ) 0


# generic functions that tie outfluxes to decomp pools
#####

# function where there are multiple fluxes from a single pool
f_decomp_fluxsum <- function(., i, ... )
  sum(unlist(.super$state$outflux[unlist(.super$pars[[paste0('decomp_outflux',i)]]) ]))


# for outfluxes structure when all pools only have one flux -- limited use
f_decomp_outflux <- function(., i, ... ) .super$state$outflux[[i]]
 


# specific decomp/outflux functions
#####

# linear decomp, Oleson 1963
# APW: not sure ... is needed
f_decomp_lin <- function(., C, t, i, of=i, ... ) 
  C[i]*.super$state_pars$k[[of]]


# non-linear Michaelis-Menten decomp 
# - typically function of microbial biomass: C[cat_pool]
# - rate saturates as a function of decomp/substrate pool mass 
f_decomp_mm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]] )
  .super$state_pars$vmax[[of]]*C[cat_pool]*C[i] / (.super$state_pars$km[[of]] + C[i])


# reverse Michaelis-Menten decomp
# - rate saturates as a function of cat_pool mass 
# APW: not sure ... is needed
f_decomp_rmm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... )  
  .super$state_pars$vmax[[of]]*C[cat_pool]*C[i] / (.super$state_pars$km[[of]] + C[cat_pool])


# double Michaelis-Menten decomp
# - rate saturates as a function of both cat_pool and decomp/substrate pool mass 
f_decomp_dmm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]] )
  .super$state_pars$vmax[[of]] * (C[i]/(.super$state_pars$km[[of]] + C[i])) * (C[cat_pool]/(.super$state_pars$km2[[of]]+C[cat_pool]))


# CORPSE reverse M-M with half-saturation a function of the microbial biomass pool relative to the sum of the unprotected pools
# APW: need a method to make indexing of pools flexible 
f_decomp_rmm_sulman <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... ) 
  (.super$state_pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[cat_pool] + .super$state_pars$km[[of]]*(C[1]+C[2]+C[3]))


# as above but protected pools as multiplier on Km
# APW: need a way to pass multiple pools
f_decomp_rmm_sulman_protected <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... ) 
  (.super$state_pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[cat_pool] + .super$state_pars$km[[of]]*(C[5]+C[6]+C[7]))


# michaelis-menten version of above decomp function
# APW: f_decomp_rmm_sulman_moore was this type of func 
f_decomp_mm_sulman <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... ) 
  (.super$state_pars$vmax[[i]]*C[cat_pool]*C[i]) / (C[i] + .super$state_pars$km[[i]]*(C[1]+C[2]+C[3]))


# density-dependent turnover, often for microbes, e.g. Georgiou et al. 2017  
f_decomp_dd_georgiou <- function(., C, t, i, of=i ) 
  (C[i]^.super$pars$beta) * .super$state_pars$k[[of]] 


# CORPSE microbial turnover
# APW: CAUTION: used to be divide by k as k was really a mean residence time 
f_micturn_sulman <- function(., C, t, i, of=i ) 
  (C[i] - .super$pars$minmic * (C[1]+C[2]+C[3])) * .super$state_pars$k[[of]] 


# as above but with density dependence
f_micturn_sulman_dd <- function(., C, t, i, of=i ) 
  (C[i]^.super$pars$beta - .super$pars$minmic * (C[1]+C[2]+C[3])) * .super$state_pars$k[[of]] 


# saturating sorption
f_sorp_sat <- function(., C, t, i, of=i, sat_pool=.super$pars$sat_pool )
  C[i]*.super$state_pars$k[[of]] * (1-C[sat_pool]/.super$state_pars$poolmax[[sat_pool]])


# desorption from a saturating pool
# - i in this case should be the same as sat_pool in the above
# - because this is desorption from the pool that saturates
f_desorp_sat <- function(., C, t, i, of=i, ... ) 
  .super$state_pars$k[[of]] * (C[i]/.super$state_pars$poolmax[[i]])



# transfer functions
################################

# transfer all or nothing
f_transfer_all  <- function(., ... ) 1
f_transfer_zero <- function(., ... ) 0 # APW: not needed

# CUE or carbon transfer efficiency sets transfer from one pool to another
# APW: these may be obsolete but can keep for now for back compatibility
f_transfer_cue            <- function(., C, t, from, to ) .super$state_pars$cue[[from]]    
f_transfer_cue2           <- function(., C, t, from, to ) .super$state_pars$cue2[[from]]    
# transfers remainder of cue function to another pool instead of CO2 loss
f_transfer_cue_remainder  <- function(., C, t, from, to ) 1 - .super$state_pars$cue[[from]]    
f_transfer_cue2_remainder <- function(., C, t, from, to ) 1 - .super$state_pars$cue[[from]] - .super$state_pars$cue2[[from]]    


# CUE or carbon transfer efficiency sets transfer from one pool to another
# - compatible with outfluxes method when only one flux from a pool
# - obtains appropriate flux id (of) for the pool id (from)
# APW: is there a way to tidy these up and call from wrapper functions as is done below? 
f_transfer_fluxid_cue <- function(., C, t, from, to ) {
  # get outflux index
  of <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  .super$state_pars$cue[[of]]
}

f_transfer_fluxid_cue2 <- function(., C, t, from, to ) {
  # get outflux index
  of <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  .super$state_pars$cue2[[of]]
}

f_transfer_fluxid_cue_remainder <- function(., C, t, from, to ) {
  # get outflux index
  of <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  1 - .super$state_pars$cue[[of]]
}

f_transfer_fluxid_cue2_remainder <- function(., C, t, from, to ) {
  # get outflux index
  of <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  1 - .super$state_pars$cue[[of]] - .super$state_pars$cue2[[of]]
}

# as above but calculates remainder from CUE & CUE2, then applies CUE3 
f_transfer_fluxid_cue2_remainder_cue3 <- function(., C, t, from, to ) {
  # get outflux index
  of   <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  rem1 <- 1 - .super$state_pars$cue[[of]] - .super$state_pars$cue2[[of]]
  .super$state_pars$cue3[[of]] * rem1
}

# caluclates remainder of above
f_transfer_fluxid_cue2_remainder_cue3_remainder <- function(., C, t, from, to ) {
  # get outflux index
  of   <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  rem1 <- 1 - .super$state_pars$cue[[of]] - .super$state_pars$cue2[[of]]
  (1 - .super$state_pars$cue3[[of]]) * rem1
}



# transfer function that when there are multiple out fluxes, calculates the proportion of each flux of the total flux
# - not applied directly, called by below wrapper functions
# - can apply CUE
f_transfer_fluxsum_prop <- function(., C, t, from, to, f, cue=F ) {
  #print("transfer fluxsum generic")
  #print(.); print(C); print(t); print(to); print(from)
  
  # get outflux index
  of <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[f]

  # if CUE get CUE
  cue_val <- if(cue) .super$state_pars$cue[[of]] else 1    
 
  total_flux <- .[[paste0('decomp.d',from)]](C=C, t=t, i=from )
  # when all fluxes are zero these functions create a NaN
  if(total_flux==0) {
    return(0)
  } else {
    cue_val * .super$state$outflux[[of]] / total_flux
  }
}


# wrapper functions for above
# APW: it's probably possible to combine these with above and create a case where total flux is not calculated if of extent is 1
f_transfer_fluxsum_prop_one       <- function(., ... ) .$transfer_fluxsum_prop(f=1, ... ) 
f_transfer_fluxsum_prop_two       <- function(., ... ) .$transfer_fluxsum_prop(f=2, ... ) 
f_transfer_fluxsum_prop_three     <- function(., ... ) .$transfer_fluxsum_prop(f=3, ... ) 
f_transfer_fluxsum_prop_one_cue   <- function(., ... ) .$transfer_fluxsum_prop(f=1, cue=T, ... ) 
f_transfer_fluxsum_prop_two_cue   <- function(., ... ) .$transfer_fluxsum_prop(f=2, cue=T, ... ) 
f_transfer_fluxsum_prop_three_cue <- function(., ... ) .$transfer_fluxsum_prop(f=3, cue=T, ... ) 


# wrapper of above, for MEND where there are two fluxes going to one pool and both have the same CUE
# APW: could also modify above to handle two outfluxes
f_transfer_fluxsum_prop_one_two_cue <- function(., ... ) 
  .$transfer_fluxsum_prop(f=1, cue=T, ... ) + .$transfer_fluxsum_prop(f=2, cue=T, ... )  

f_transfer_fluxsum_prop_two_three_cue <- function(., ... ) 
  .$transfer_fluxsum_prop(f=2, cue=T, ... ) + .$transfer_fluxsum_prop(f=3, cue=T, ... )  


# CUE / transfer efficiency sets transfer subject to a maximum pool size 
# - can be used both for saturating MAOM pool and density dependent microbial growth efficiency
# APW: should this include a min 0 function in the final term?
# APW: instead of modifying flux by pool size, this modifies CUE, needs work to be compatible
f_transfer_cue_sat <- function(.,C,t,from,to) 
  .super$pars$cue[[from]] * (1-C[to]/.super$pars$poolmax[[to]])



# APW: decomp functions not yet converted to outfluxes methods and state_pars lists 
###############################

# alternative CORPSE hypotheses
#########

## adding a decomp function that saturates and returns excess to unprotected pool
## APW: can make part of this a poolmax function and combine with f_sorp_sat, but need a way to specify more than one sat pool
#f_decomp_lin_sat_corpse        <- function(.,C,t,i) {
#  #match millennialv2 for ESA2022 sims
#  if(.super$env$clay == 27){ 
#    protected_max = 6.837
#  } else protected_max = 10.32
#  # protected_max = .super$pars$poolmax[[5]] + .super$pars$poolmax[[6]] + .super$pars$poolmax[[7]] #this is how to set it generally
#    C[i]*.super$pars$k[[i]]* (1-(C[5] + C[6] + C[7])/protected_max)
#}



# MILLENNIAL-specific functions
#########

# aggregate formation
f_aggform_abramoff <- function(.,C,t,i, agg_pool = 5){
  if(i==1){
     # print(C[i])
     # print(.super$pars$millennial[['Vpa']])
     # print(.super$pars$millennial[['Kpa']])
     # print(.super$pars$poolmax[[5]])
    ( .super$pars$millennial[['Vpa']] * C[i]) / (.super$pars$millennial[['Kpa']] + C[i])*(1 - C[agg_pool] / .super$pars$poolmax[[5]])
  }
  else if(i==3){
    ( .super$pars$millennial[['Vma']] * C[i]) / (.super$pars$millennial[['Kma']] + C[i])*(1 - C[agg_pool] / .super$pars$poolmax[[5]])  
  }
}

# sorption of DOC to soil minerals
f_docsorp_abramoff <- function(.,C,t,i){
  Qmax = (.super$env$BD * 10^(.super$pars$millennial[['c1']]*log10(.super$env$clay) + .super$pars$millennial[['c2']]))/1000
  Kdm = 10^(-.186*.super$env$pH - .216)
  C[i] * ((Kdm*Qmax*C[i])/(1+Kdm*C[i]) - C[3]) / Qmax
}



# MEND 2020 specific functions
#########

# MEND losses of microbial biomass toward CO2 (Fr) and toward other pools (Fe)
# MEND O2
# APW: similar to MEND 2013 respiraton, but different so keeping for now 
f_decomp_mbc_mend <- function(.,C,t,i) {
  ( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
    C[2]*.super$pars$mr #Fe
}

# MEND losses of dissolved organic matter toward mbc uptake (Fu) and adsorption to mineral surfaces (Fa)
# MEND O5
f_decomp_doc_mend <- function(.,C,t,i) {
  C[5]*((.super$pars$vmax[[i]]+.super$pars$mr)/.super$pars$cue[[i]])*(C[2]/(.super$pars$km[[i]]+C[5])) +  #Fu
  C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]])             #Fa
}


# MIMICS-specific functions
#########
# alt mimics function
# APW: this can be converted if needed
# APW: can be broken into two fluxes that have a different cat pool but share the same parameters 
#f_decomp_rmm_twocat <- function(.,C,t,i, cat1 = 3, cat2 = 4) 
#  (.super$pars$vmax[[i]]*(C[cat1]+C[cat2])*C[i]) / ((C[cat1]+C[cat2]) + .super$pars$km[[i]])



### END ###
