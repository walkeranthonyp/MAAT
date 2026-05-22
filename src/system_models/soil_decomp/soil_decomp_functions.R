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
  # why is .$env$litter not .super$env$litter??
  .$env$litter * matrix(unlist(.super$pars$input_coef)[1:.super$pars$n_pools], ncol=1 )
}


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
f_input_mimics <- function(., t ) {
  # convert anpp to litter inputs
  EST_LIT_in <- .super$env$anpp / (365*24) # gC/m2/h (from gC/m2/y)
  EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h)
  input      <- .super$env$litter <- EST_LIT/.super$env$depth 
 
  # partitioning inputs
  i_LITm <- input * .super$state_pars$fmet     * (1-.super$pars$mimics[['fi_LITm']])
  i_LITs <- input * (1-.super$state_pars$fmet) * (1-.super$pars$mimics[['fi_LITs']])
  i_SOMp <- input * .super$state_pars$fmet     * .super$pars$mimics[['fi_LITm']]
  i_SOMc <- input * (1-.super$state_pars$fmet) * .super$pars$mimics[['fi_LITs']]

  # input matrix
  matrix(c(i_LITm, i_LITs, 0, 0, i_SOMp, i_SOMc, 0 ), ncol=1 ) 
}



# decomp functions
# APW: we might want a broader class of functions called flux functions given we have uptake, sorption, and other processes 
#      that cause a flux from one pool to another 
################################

# generic functions

# APW: both are not needed
f_decomp_none <- function(., ... ) 0
f_zero        <- function(., ... ) 0
# APW: where is this used? 
f_identity    <- function(., C, t, i ) 1


# function where there are multiple fluxes from a single pool
f_decomp_fluxsum <- function(., i, ... )
  sum(unlist(.super$state$outflux[unlist(.super$pars[[paste0('decomp_outflux',i)]]) ]))


# for outfluxes structure when all pools only have one flux
f_decomp_outflux <- function(., i, ... ) .super$state$outflux[[i]]
 

# specific decomp functions

# linear decomp, Oleson 1963
# APW: not sure ... is needed
f_decomp_lin <- function(., C, t, i, of=i, ... ) 
  C[i]*.super$state_pars$k[[of]]


# density-dependent decomp, used for density-dependent turnover as in Georgiou et al. 2017  
f_decomp_dd_georgiou <- function(., C, t, i, of=i ) 
  (C[i]^.super$pars$beta) * .super$state_pars$k[[of]] 


# non-linear Michaelis-Menten decomp 
# - typically function of microbial biomass: C[cat_pool]
# - rate saturates as a function of decomp/substrate pool mass 
f_decomp_mm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]] )
  .super$state_pars$vmax[[of]]*C[cat_pool]*C[i] / (.super$pars$km[[of]] + C[i])

f_decomp_MM_microbe <- f_decomp_mm 


# reverse Michaelis-Menten decomp
# - rate saturates as a function of cat_pool mass 
# APW: not sure ... is needed
f_decomp_rmm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... )  
  .super$state_pars$vmax[[of]]*C[cat_pool]*C[i] / (.super$state_pars$km[[of]] + C[cat_pool])


# APW: not converted to new structure as yet 
# double Michaelis-Menten decomp
# - rate saturates as a function of cat_pool and decomp/substrate pool mass 
f_decomp_dmm <- function(., C, t, i )
  .super$pars$vmax[[i]] * (C[1]/(.super$pars$km[[1]] + C[1])) * (C[2]/(.super$pars$rkm[[1]]+C[2]))


# saturating sorption
f_sorp_sat <- function(., C, t, i, of=i, sat_pool=.super$pars$sat_pool )
  C[i]*.super$state_pars$k[[of]] * (1-C[sat_pool]/.super$state_pars$poolmax[[sat_pool]])


# desorption from a saturating pool
# - i in this case should be the same as sat_pool in the above
# - because this is desorption from the pool that saturates
f_desorp_sat <- function(., C, t, i, of=i, ... ) 
  .super$pars$k[[of]] * C[i] / .super$pars$poolmax[[i]]






# CORPSE-specific functions
#########
f_decomp_rmm_sulman <- function(.,C,t,i) 
  (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[4] + .super$pars$km[[i]]*(C[1]+C[2]+C[3]))

f_decomp_rmm_sulman_gen <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... ) 
  #(.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[cat_pool] + .super$pars$km[[of]]*(C[1]+C[2]+C[3]))
  (.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[cat_pool] + .super$state_pars$km[[of]]*(C[1]+C[2]+C[3]))

f_micturn_sulman <- function(., C, t, i, of=i ) 
  #(C[i] - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[of]] 
  (C[i] - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$state_pars$k[[of]] 


f_micturn_sulman_mdd <- function(.,C,t,i) 
  (C[i]^.super$pars$beta - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[i]] 




# alternative CORPSE hypotheses
#########

# michaelis-menten version of decomp function
# would need to change km parameterization
f_decomp_mm_sulman <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[i] + .super$pars$km[[i]]*(C[1]+C[2]+C[3]))

# CORPSE density dependence version of microbial turnover
f_micturn_sulman_dd <- function(.,C,t,i) {
  (C[i]^.super$pars$beta - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[i]] 
  # ((C[i] - .super$pars$minmic * (C[1]+C[2]+C[3])))^.super$pars$beta/.super$pars$k[[i]]
  # C[i]^.super$pars$beta/.super$pars$k[[i]] 
}

# adding a decomp function that saturates and returns excess to unprotected pool
f_decomp_lin_sat_corpse        <- function(.,C,t,i) {
  #match millennialv2 for ESA2022 sims
  if(.super$env$clay == 27){ 
    protected_max = 6.837
  } else protected_max = 10.32
  # protected_max = .super$pars$poolmax[[5]] + .super$pars$poolmax[[6]] + .super$pars$poolmax[[7]] #this is how to set it generally
    C[i]*.super$pars$k[[i]]* (1-(C[5] + C[6] + C[7])/protected_max)
}

# Improvement to CORPSE described in Moore et al. 2019, scales microbial biomass to each specific substrate type
# This type of function would be more general to protected pools. 
f_decomp_rmm_sulman_moore <- function(.,C,t,i) .super$pars$vmax[[i]]*C[i] * C[4] / (C[i] + .super$pars$km[[i]]*(C[1]+C[2]+C[3]))

# f_desorp_lin_sat_corpse <- function(.,C,t,i) {
#   C[i]*.super$pars$k[[i]]
# }

# decay of protected pool
f_decomp_rmm_sulman_protected <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[4] + .super$pars$km[[i]]*(C[5]+C[6]+C[7]))



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

# uptake of DOC by micropbes
# reverse michaelis-menten equation so could potentially use a f_decomp_rmm function instead
f_docuptake_abramoff <- function(., C, t, i )
  .super$pars$millennial[['Vdm']] * C[i] * C[2]/(C[2]+.super$pars$millennial[['Kdb']])


# linear
f_uptake_lin <- function(., C, t, i, cat )
  .super$pars$vmax[[i]] * C[i]




# MIMICS-specific functions
#########
# environmental controls are embedded,
# would be ideal to pull these out
# for now we could define all the parameters with if...then statements in the solver function


# alt mimics function
f_decomp_rmm_twocat <- function(.,C,t,i, cat1 = 3, cat2 = 4) 
  (.super$pars$vmax[[i]]*(C[cat1]+C[cat2])*C[i]) / ((C[cat1]+C[cat2]) + .super$pars$km[[i]])



# MILLENNIALV2-specific functions 
#########

# cat is dummy
# APW: can we use ... ? cat seems unnecessary, can we remove?
#f_desorp_millennialv2 <- function(., C, t, i, k=1, cat=NULL ) { 
#  k * C[i] / .super$pars$poolmax[[i]]
# APW: can be modified to pull k from list using i/of indexing
f_desorp_millennialv2 <- function(., C, t, i, of=i, k=.super$pars$millennialV2[['kld']]/1000, ... )  
  k*(C[i]/.super$pars$poolmax[[i]])
  #if(C[i]>0) k * C[i] / .super$pars$poolmax[[i]] else 0 
f_desorp_millennialv2 <- function(., C, t, i, of=i, ... )
  .super$state_pars$k[[of]]*(C[i]/.super$state_pars$poolmax[[i]])


# cat is dummy
# APW: can we use ... ?
#f_desorp_millennialv2_nosat <- function(., C, t, i, k=1, cat=NULL ) 
f_desorp_millennialv2_nosat <- function(., C, t, i, ... ) 
  .super$pars$millennialV2[['kld']] * C[i]



# transfer functions
################################

# transfer all or nothing
f_transfer_all  <- function(., ... ) 1
f_transfer_zero <- function(., ... ) 0

# CUE or carbon transfer efficiency sets transfer from one pool to another
# MEND15
f_transfer_cue       <- function(., C, t, from, to ) .super$pars$cue[[from]]    
f_transfer_cue2      <- function(., C, t, from, to ) .super$pars$cue2[[from]]    
f_transfer_cue_resp  <- function(., C, t, from, to ) 1 - .super$pars$cue[[from]]    
f_transfer_cue       <- function(., C, t, from, to ) .super$state_pars$cue[[from]]    
f_transfer_cue2      <- function(., C, t, from, to ) .super$state_pars$cue2[[from]]    
f_transfer_cue_resp  <- function(., C, t, from, to ) 1 - .super$state_pars$cue[[from]]    

# transfers remainder of cue function to another pool instead of CO2
# MEND13
f_transfer_cue_remainder <- function(.,C,t,from,to) (1-.super$pars$cue[[from]])
f_transfer_cue_remainder <- function(.,C,t,from,to) (1-.super$state_pars$cue[[from]])
f_transfer_cue_remainder <- function(.,C,t,from,to) (1-.super$state_pars$cue[[from]])

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

f_transfer_fluxid_cue2_remainder_cue3 <- function(., C, t, from, to ) {
  # get outflux index
  of   <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  rem1 <- 1 - .super$state_pars$cue[[of]] - .super$state_pars$cue2[[of]]
  .super$state_pars$cue3[[of]] * rem1
}

f_transfer_fluxid_cue2_remainder_cue3_remainder <- function(., C, t, from, to ) {
  # get outflux index
  of   <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[1]
  rem1 <- 1 - .super$state_pars$cue[[of]] - .super$state_pars$cue2[[of]]
  (1 - .super$state_pars$cue3[[of]]) * rem1
}



# transfer functions that when there are multiple out fluxes, calculate the proportion of each flux of the total flux
# APW: could turn this into a standard calculation in the SoilR function that calculates proportion for all fluxes 
# APW problem: when all fluxes are zero these functions create a NaN
f_transfer_fluxsum_prop <- function(., C, t, from, to, f, cue=F ) {
  #print("transfer fluxsum generic")
  #print(.); print(C); print(t); print(to); print(from)
  
  # get outflux index
  of <- unlist(.super$pars[[paste0('decomp_outflux',from)]])[f]

  # if CUE get CUE
  #cue_val <- if(cue) .super$pars$cue[[of]] else 1    
  cue_val <- if(cue) .super$state_pars$cue[[of]] else 1    
 
  total_flux <- .[[paste0('decomp.d',from)]](C=C, t=t, i=from )
  if(total_flux==0) {
    return(0)
  } else {
    cue_val * .super$state$outflux[[of]] / total_flux
  }
}

f_transfer_fluxsum_prop_one   <- function(., ... ) .$transfer_fluxsum_prop(f=1, ... ) 
f_transfer_fluxsum_prop_two   <- function(., ... ) .$transfer_fluxsum_prop(f=2, ... ) 
f_transfer_fluxsum_prop_three <- function(., ... ) .$transfer_fluxsum_prop(f=3, ... ) 
f_transfer_fluxsum_prop_one_cue   <- function(., ... ) .$transfer_fluxsum_prop(f=1, cue=T, ... ) 
f_transfer_fluxsum_prop_two_cue   <- function(., ... ) .$transfer_fluxsum_prop(f=2, cue=T, ... ) 
f_transfer_fluxsum_prop_three_cue <- function(., ... ) .$transfer_fluxsum_prop(f=3, cue=T, ... ) 

f_transfer_fluxsum_prop_one_two_cue <- function(., ... ) 
  .$transfer_fluxsum_prop(f=1, cue=T, ... ) + .$transfer_fluxsum_prop(f=2, cue=T, ... )  


# APW: this is an issue that needs solved
#      we don't need n fluxes functions for each transfer coefficient
#      can something like above be done so it's just one function with n fluxes wrappers
# APW: also this includes a temperature function to calculate CUE, could be a state_pars type calculation   
#f_transfer_fluxsum_prop_three_cue <- function(., from, ... ) 
#  (.super$pars$cue[[from]] - .super$pars$millennialV2[['cue_t']] * (.super$env$temp - .super$pars$millennialV2[['Taeref']])) * 
#  .$transfer_fluxsum_prop(f=3, from=from, ... ) 

# CUE / transfer efficiency sets transfer subject to a maximum pool size 
# - can be used both for saturating MAOM pool and density dependent microbial growth efficiency
# APW: should this include a min 0 function in the final term?
f_transfer_cue_sat <- function(.,C,t,from,to) 
  .super$pars$cue[[from]] * (1-C[to]/.super$pars$poolmax[[to]])



###########################
# functions not yet updated (or maybe even needed) for new code

# MEND 2020? specific functions

# non-linear Michaelis-Menten decomp of POM, a function of POM-specific enzymes
# MEND O1
# APW: why is this needed instead of f_decomp_mm with cat=6?
f_decomp_MM_enzpom <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[6]*C[i]) / (.super$pars$km[[i]]+C[i])

# MEND losses of microbial biomass toward CO2 (Fr) and toward other pools (Fe)
# MEND O2
# APW: seems to include CO2 transfer coeffcient, why can't that be specified in transfer matrix?
f_decomp_mbc_mend <- function(.,C,t,i) {
  ( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
    C[2]*.super$pars$mr #Fe
}

# non-linear Michaelis-Menten decomp of MAOM, a function of MAOM-specific enzymes
# MEND O3
# APW: why is this needed instead of f_decomp_mm with cat=7?
f_decomp_MM_enzmaom <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[7]*C[i]) / (.super$pars$km[[i]]+C[7])

# density-dependent decomp, as used in MEND for desorption from the Q pool
# MEND O4
f_decomp_dd_mend <- function(.,C,t,i) .super$pars$k[[i]]*(C[i]/.super$pars$poolmax[[i]])

# MEND losses of dissolved organic matter toward mbc uptake (Fu) and adsorption to mineral surfaces (Fa)
# MEND O5
f_decomp_doc_mend <- function(.,C,t,i) {
  C[5]*((.super$pars$vmax[[i]]+.super$pars$mr)/.super$pars$cue[[i]])*(C[2]/(.super$pars$km[[i]]+C[5])) +  #Fu
  C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]])             #Fa
}


# APW: not sure where these MEND-speciifc transfer functions have been used
# transfer from mbc to pom in mend: (1-gd)(1-pe)(fe/f_decomp_mbc_mend)
f_transfer_mend21 <- function(.,C,t,from,to){
  (1-.super$pars$cue[[from]])*
  (1-(.super$pars$pep+.super$pars$pem))*
  ((C[2]*.super$pars$mr)/
     (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
     C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend25 <- function(.,C,t,from,to){
  (.super$pars$cue[[from]])*
    (1-(.super$pars$pep+.super$pars$pem))*
    ((C[2]*.super$pars$mr)/
       (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
       C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend26 <- function(.,C,t,from,to) {
  .super$pars$pep * 
    ((C[2]*.super$pars$mr)/
       (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
       C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend27 <- function(.,C,t,from,to) {
  .super$pars$pem * 
    ((C[2]*.super$pars$mr)/
       (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
       C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend52 <- function(.,C,t,from,to){
 # print(C)
  #print(from)
  #print(to)
  #print(t)
  #print(.super$pars$poolmax[[4]])
  #print(C[5])
  #print(C[4])
  (C[5]*((.super$pars$vmax[[from]]+.super$pars$mr)/.super$pars$cue[[from]])*(C[2]/(.super$pars$km[[from]]+C[5]))) / #Fu
  (C[5]*((.super$pars$vmax[[from]]+.super$pars$mr)/.super$pars$cue[[from]])*(C[2]/(.super$pars$km[[from]]+C[5])) +  #Fu
    C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]]))   #Fa
}

f_transfer_mend54 <- function(.,C,t,from,to){
  (C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]]))/ #Fu
    (C[5]*((.super$pars$vmax[[from]]+.super$pars$mr)/.super$pars$cue[[from]])*(C[2]/(.super$pars$km[[from]]+C[5])) +  #Fu
       C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]]))   #Fa
}





# APW: these have all been replaced by more generic functions in the new code

# doc uptake 
# APW: some kind of MM
#f_docuptake_mend <- function(., C, t, i )
#  (.super$pars$vmax[[i]]+.super$pars$mend[['Mr']])*C[i]*C[4] / (C[i]+.super$pars$km[[i]])

# micrbial growth respiration
# APW: some kind of MM
#f_growthresp_mend <- function(., C, t, i )
  #.super$pars$vmax[[5]]*C[i]*C[5] / (C[5]+.super$pars$km[[5]])
  # APW: 5 is the DOC pool, and DOC uptake is the 9th of
  # APW: i should be 4, the Microbial pool
  #.super$pars$vmax[[9]]*C[i]*C[5] / (C[5]+.super$pars$km[[9]])


# micrbial maintenance respiration
# APW: some kind of MM
#f_maintresp_mend <- function(., C, t, i )
  #.super$pars$k[[4]]*C[i]*C[5] / (C[5]+.super$pars$km[[5]])
  # APW: 5 is the DOC pool, and DOC uptake is the 9th olooks f
  # APW: i should be 4, the Microbial pool
  # APW: k[[4]] is total microbial turnover inc. enzyme production, which is of 6,7,8
  #.super$pars$k[[6]]*C[i]*C[5] / (C[5]+.super$pars$km[[9]])


# APW: WFT?? haha
#f_decomp_rmm_wieder <- function(.,C,t,i,cat_pool = 3){
#  if(cat_pool == 3){
#    if(i==7){
#      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
#      Km = .super$pars$km[[7]] * pscalar
#    } else if(i==6) {
#      #km par same as structural litter (km2)
#      Km = .super$pars$km[[2]] *(1/.super$pars$mimics[['ko_r']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
#    }  else {
#      Km = .super$pars$km[[i]]
#    }
#    #correcting Km for temperature
#    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
#    ### RMM equation
#    C[i] * .super$pars$vmax[[i]] * C[cat_pool] / (Km_cor + C[cat_pool])
#    ###
#
#  } else {
#    if(i==7){
#      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
#      Km = .super$pars$km2[[7]] * pscalar
#    } else if(i==6) {
#      Km = .super$pars$km2[[2]] *(1/.super$pars$mimics[['ko_k']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
#    }  else {
#      Km = .super$pars$km2[[i]]
#    }
#    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
#    ###RMM equation
#    C[i] * .super$pars$vmax2[[i]] * C[cat_pool] / (Km_cor + C[cat_pool])
#    ###
#  }
#}
#
#f_decomp_mm_wieder <- function(.,C,t,i,cat_pool = 3){
#  if(cat_pool == 3){
#    if(i==7){
#      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
#      Km = .super$pars$km[[7]] * pscalar
#    } else if(i==6) {
#      #km par same as structural litter (km2)
#      Km = .super$pars$km[[2]] *(1/.super$pars$mimics[['ko_r']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
#    }  else {
#      Km = .super$pars$km[[i]]
#    }
#    #correcting Km for temperature
#    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
#    ###MM equation
#    C[i] * .super$pars$vmax[[i]]*10 * C[cat_pool] / (Km_cor*10 + C[i])
#    ###
#  } else {
#    if(i==7){
#      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
#      Km = .super$pars$km2[[7]] * pscalar
#    } else if(i==6) {
#      Km = .super$pars$km2[[2]] *(1/.super$pars$mimics[['ko_k']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
#    }  else {
#      Km = .super$pars$km2[[i]]
#    }
#    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
#    ###MM equation
#    C[i] * .super$pars$vmax2[[i]]*10 * C[cat_pool] / (Km_cor*10 + C[i])
#    ###
#  }
#}


### END ###
