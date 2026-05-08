################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# Matt Craig, AWalker October 2019 
# Carlos Sierra, Markus Mueller (SoilR developers) 
#
################################

source('soil_decomp_temperature_functions.R')
source('soil_decomp_water_functions.R')


# FUNCTIONS
################################


# litter functions 
################################

# function for converting ANPP to litter inputs from MIMICS
f_input_mimics <- function(., t ) {
  # convert inputs
  EST_LIT_in <- .super$env$anpp / (365*24) # gC/m2/h (from gC/m2/y)
  EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h)
  input      <- EST_LIT/.super$env$depth 
 
  # partitioning inputs
  i_LITm <- input * .super$state_pars$fmet     * (1-.super$pars$mimics[['fi_LITm']])
  i_LITs <- input * (1-.super$state_pars$fmet) * (1-.super$pars$mimics[['fi_LITs']])
  i_SOMp <- input * .super$state_pars$fmet     * .super$pars$mimics[['fi_LITm']]
  i_SOMc <- input * (1-.super$state_pars$fmet) * .super$pars$mimics[['fi_LITs']]

  # input matrix
  matrix(c(i_LITm, i_LITs, 0, 0, i_SOMp, i_SOMc, 0 ), ncol=1 ) 
}



# decomp functions
# APW: strikes me that we might want a broader class of functions called flux functions given we have uptake, sorption, and other processes that cause a flux from one pool to another 
###################

# generic functions
f_decomp_none <- function(., C, t, i, ... ) 0

f_zero        <- function(., C, t, i ) 0

# APW: where is this used? 
f_identity    <- function(., C, t, i ) 1

# function where there are multiple fluxes from a single pool
f_decomp_fluxsum <- function(., i, ... )
  sum(unlist(.super$state$outflux[unlist(.super$pars[[paste0('decomp_outflux',i)]]) ]))
 
# function where there is a single flux from each pool, to work with outfluxes structure when all pools only have one flux
f_decomp_outflux <- function(., i, ... ) .super$state$outflux[[i]]
 
# linear decomp, Oleson 1963
# - cat and sat_pool are dummy arguments to allow switching from other functions, APW: can ... be used?
# - MEND O6,O7 (APW: not sure what this is referring to)
# APW: I think k_from_list can be removed and the pars k specified as default
#f_decomp_lin <- function(., C, t, i, k_from_list=TRUE, k=NULL, cat=NULL, sat_pool=NULL ) { 
#  if(k_from_list == TRUE){
#    C[i]*.super$pars$k[[i]]
#  } else {
#    C[i]*k
#  }
#}
f_decomp_lin <- function(., C, t, i, of=i, k=.super$pars$k[[of]], ... ) 
  C[i]*k
  #if(C[i]>0) C[i]*k else 0
f_decomp_lin <- function(., C, t, i, of=i, k=.super$state_pars$k[[of]], ... ) 
  C[i]*k


# density-dependent decomp, used for density-dependent turnover as in Georgiou et al. 2017  
f_decomp_dd_georgiou <- function(., C, t, i, of=i ) 
  (C[i]^.super$pars$beta) * .super$pars$k[[of]] 
  #if(C[i]>0) (C[i]^.super$pars$beta) * .super$pars$k[[of]] else 0      


# non-linear Michaelis-Menten decomp, a function of microbial biomass
# APW: this is the same as below
#f_decomp_MM_microbe <- function(., C, t, i ) 
#  (.super$pars$vmax[[i]]*C[2]*C[i]) / (.super$pars$km[[i]]+C[i])
f_decomp_MM_microbe <- function(., C, t, i, cat_pool=.super$pars$cat_pool ) 
  (.super$pars$vmax[[i]]*C[cat_pool]*C[i]) / (.super$pars$km[[i]]+C[i])


# Michaelis-Menten decomp
#f_decomp_mm <- function(., C, t, i, cat=4 )
#  (.super$pars$vmax[[i]]*C[cat]*C[i]) / (C[i] + .super$pars$km[[i]])
# - i refers to pool, of refers to outflux, where n_pools = n_outfluxes of = i 
f_decomp_mm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool )
  (.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[i] + .super$pars$km[[of]])
  #if(C[i]>0 & C[cat_pool]>0) (.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[i] + .super$pars$km[[of]]) else 0


# reverse Michaelis-Menten decomp
# - C[4] is microbial biomass in CORPSE
# - k is dummy argument to get MILLENNIAL to work for mic decay of maom
# - APW: should km be rkm here or does it not matter? 
#f_decomp_rmm <- function(., C, t, i, cat=4, k=NULL )
#  (.super$pars$vmax[[i]]*C[cat]*C[i]) / (C[cat] + .super$pars$km[[i]])
#f_decomp_rmm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool, ... ) { 
f_decomp_rmm <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool[[of]], ... ) { 
  #print("decomp_rmm")
  #print(i)
  #print(of)
  #print(C[i])
  #print(cat_pool) 
  #print(C[cat_pool]) 
  #print(.super$pars$vmax[[of]])
  #print(.super$state_pars$km[[of]]) 
  #(.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (.super$pars$km[[of]] + C[cat_pool])
  (.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (.super$state_pars$km[[of]] + C[cat_pool])
  #if(C[i]>0 & C[cat_pool]>0) (.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (.super$pars$km[[of]] + C[cat_pool]) else 0
}

# double Michaelis-Menten decomp
# APW: curvature of is a function of both substrate concentration and microbial pool size?
f_decomp_dmm <- function(., C, t, i )
  .super$pars$vmax[[i]] * (C[1]/(.super$pars$km[[1]] + C[1])) * (C[2]/(.super$pars$rkm[[1]]+C[2]))


# saturating sorption
#f_sorp_sat  <- function(., C, t, i, k_from_list = TRUE, k = NULL, sat_pool ) {
#  if(k_from_list == TRUE) {
#    C[i]*.super$pars$k[[i]]*(1-C[sat_pool]/.super$pars$poolmax[[sat_pool]])
#  } else {
#    C[i]*k*(1-C[sat_pool]/.super$pars$poolmax[[sat_pool]])
#  }
#}
#f_sorp_sat <- function(., C, t, i, k=.super$pars$k[[i]], sat_pool=.super$pars$sat_pool ) 
#  C[i]*k*(1-C[sat_pool]/.super$pars$poolmax[[sat_pool]])
# APW: this can be modified to access k from list with of/i indexing 
# APW: I think this might be a monod functino just written a different way
f_sorp_sat <- function(., C, t, i, of=i, k=.super$state_pars$kaff_lm, sat_pool=.super$pars$sat_pool )
  C[i]*k*(1-C[sat_pool]/.super$pars$poolmax[[sat_pool]])
  #if(C[i]>0 & C[sat_pool]>0) C[i]*k*(1-C[sat_pool]/.super$pars$poolmax[[sat_pool]]) else 0



# MEND-specific functions
#########

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



# CORPSE-specific functions
#########
f_decomp_rmm_sulman <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[4] + .super$pars$km[[i]]*(C[1]+C[2]+C[3]))
f_decomp_rmm_sulman_gen <- function(., C, t, i, of=i, cat_pool=.super$pars$cat_pool, ... ) 
  (.super$pars$vmax[[of]]*C[cat_pool]*C[i]) / (C[cat_pool] + .super$pars$km[[of]]*(C[1]+C[2]+C[3]))

f_micturn_sulman <- function(., C, t, i, of=i ) {
  (C[i] - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[of]] 
}

f_micturn_sulman_mdd <- function(.,C,t,i) {
  (C[i]^.super$pars$beta - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[i]] 
}



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

calc_state_pars_weider <- function(.) {

  # general parameters
  .super$state_pars$fmet     <- .$fmet()
  .super$state_pars$tau_mod1 <- .$tau_mod1()
  
  # k parameters
  # -- as MEC noted when there is a catalyst pool in mm or rmm decomp, Vmax units are per unit time only
  # -- i.e. they are a specific flux rate per unit mass of the catalyst, so mass units cancel
  # -- thus vmax is equivalent to k  
  for(i in c(5:7)) .super$state_pars$k[[i]] <- .[[paste0('k.k',i)]](i=i) 
  #.super$pars$k[[5]] <- .$k_tau_r()
  #.super$pars$k[[6]] <- .$k_tau_k()
  #.super$pars$k[[7]] <- .$k_tau_desorb()

  # km parameters
  #for(i in 1:.super$pars$n_outfluxes) .super$state_pars$km[[i]] <- .[[paste0('km.km',i)]](i)
  for(i in c(1:4,8:11)) .super$state_pars$km[[i]] <- .[[paste0('km.km',i)]](i=i) 

  # cue parameters
  for(i in c(1:6,10:11)) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
  for(i in 5:6)          .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](i=i) 
}


# dynamic parameters

# km functions
#Wieder et al. 2015 function for calculating Km using a base Km value and clay content
#cat_pool can be adjusted: 3 = r_selected microbes, 4 = k_selected microbes
#in MIMICS this function is applied to pool 7 (SOMa) for both types of microbe (r and k)
#cat_pool=.super$pars$cat_pool[[i]] ) {
f_km_wieder_temp_clay <- function(., i ) { 
  pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
  km1 = .super$pars$km[[i]] * pscalar
  exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] / km1
}

#Wieder et al. 2015 function for calculating Km using a base Km value and tuning coefficient (ko_r)
#This function currently specific to r_selected microbes (cat_pool = 3)
#in MIMICS this function is applied to pool 6 (SOMc) for both types of microbe (r and k)
f_km_wieder_temp_tuning1 <- function(., i ) {
  km1 = .super$pars$km[[i]] *(1/.super$pars$mimics[['ko_r']]) #
  exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] / km1
}

#Wieder et al. 2015 function for calculating Km using a base Km value and tuning coefficient (ko_r)
#This function currently specific to k_selected microbes (cat_pool = 4)
#in MIMICS this function is applied to pool 6 (SOMc) for both types of microbe (r and k)
f_km_wieder_temp_tuning2 <- function(., i ) {
  km1 = .super$pars$km[[i]] *(1/.super$pars$mimics[['ko_k']]) #
  exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] / km1
}

f_km_wieder_temp <- function(., i ) {
  exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /.super$pars$km[[i]]
}

# k functions
f_texture_wieder_fmet <- function(.) 
  .super$pars$mimics[['fmet_p1']] * (.super$pars$mimics[['fmet_p2']] - .super$pars$mimics[['fmet_p3']]*(.super$env$lignin/.super$env$N))

# ensures that tau_mod1 is between two values
# in Will's script, anpp is multipled by 0 in the manipulation scirpt... which would imply that 
# tau_mod1 might always be set to 0.6 in the simulations in Ben's paper
f_k_wieder_tau_mod1 <- function(.) 
  min(max(sqrt(.super$env$anpp/.super$pars$mimics[['tau_mod1_p1']]),.super$pars$mimics[['tau_mod1_p2']]),.super$pars$mimics[['tau_mod1_p3']])
f_k_wieder_tau_r <- function(., ... )
  .super$pars$mimics[['tau_r_p1']] * exp(.super$pars$mimics[['tau_r_p2']] * .super$state_pars$fmet) * .super$state_pars$tau_mod1 * .super$pars$mimics[['tau_mod2']]
f_k_wieder_tau_k <- function(., ... ) 
  .super$pars$mimics[['tau_k_p1']] * exp(.super$pars$mimics[['tau_k_p2']] * .super$state_pars$fmet) * .super$state_pars$tau_mod1 * .super$pars$mimics[['tau_mod2']]
f_k_wieder_desorb <- function(., ... ) 
  .super$pars$mimics[['desorb_p1']] * exp(.super$pars$mimics[['desorb_p2']] * .super$env$clay) * 0.1
 
f_cue_wieder_texture1 <- function(., ... ) 
  .super$pars$mimics[['fSOMp_r_p1']] * exp(.super$pars$mimics[['fSOMp_r_p2']]*.super$env$clay) * 0.1 #0.1 manual calibration from Will's script #fSOMp_r
f_cue_wieder_quality1 <- function(., ... ) 
  .super$pars$mimics[['fSOMc_r_p1']] * exp(.super$pars$mimics[['fSOMc_r_p2']]*.super$state_pars$fmet)*.super$pars$mimics[['fSOMc_r_p3']]  #fSOMc_r
#f_transfer_remainder_2cues1 <- function(., ... ) 
#  1 - .$transfer.t3_to_5() - .$transfer.t3_to_6() 
f_cue_wieder_texture2 <- function(., ... ) 
  .super$pars$mimics[['fSOMp_k_p1']] * exp(.super$pars$mimics[['fSOMp_k_p2']]*.super$env$clay) * 0.1 #fSOMp_k
f_cue_wieder_quality2 <- function(., ... ) 
  .super$pars$mimics[['fSOMc_k_p1']] * exp(.super$pars$mimics[['fSOMc_k_p2']]*.super$state_pars$fmet)*.super$pars$mimics[['fSOMc_k_p3']]  #fSOMc_k
#f_transfer_remainder_2cues2 <- function(., ... ) 
#  1 - .$transfer.t4_to_5() - .$transfer.t4_to_6() 

# APW: WFT?? haha
f_decomp_rmm_wieder <- function(.,C,t,i,cat_pool = 3){
  if(cat_pool == 3){
    if(i==7){
      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
      Km = .super$pars$km[[7]] * pscalar
    } else if(i==6) {
      #km par same as structural litter (km2)
      Km = .super$pars$km[[2]] *(1/.super$pars$mimics[['ko_r']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
    }  else {
      Km = .super$pars$km[[i]]
    }
    #correcting Km for temperature
    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
    ### RMM equation
    C[i] * .super$pars$vmax[[i]] * C[cat_pool] / (Km_cor + C[cat_pool])
    ###

  } else {
    if(i==7){
      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
      Km = .super$pars$km2[[7]] * pscalar
    } else if(i==6) {
      Km = .super$pars$km2[[2]] *(1/.super$pars$mimics[['ko_k']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
    }  else {
      Km = .super$pars$km2[[i]]
    }
    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
    ###RMM equation
    C[i] * .super$pars$vmax2[[i]] * C[cat_pool] / (Km_cor + C[cat_pool])
    ###
  }
}

f_decomp_mm_wieder <- function(.,C,t,i,cat_pool = 3){
  if(cat_pool == 3){
    if(i==7){
      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
      Km = .super$pars$km[[7]] * pscalar
    } else if(i==6) {
      #km par same as structural litter (km2)
      Km = .super$pars$km[[2]] *(1/.super$pars$mimics[['ko_r']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
    }  else {
      Km = .super$pars$km[[i]]
    }
    #correcting Km for temperature
    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
    ###MM equation
    C[i] * .super$pars$vmax[[i]]*10 * C[cat_pool] / (Km_cor*10 + C[i])
    ###
  } else {
    if(i==7){
      pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
      Km = .super$pars$km2[[7]] * pscalar
    } else if(i==6) {
      Km = .super$pars$km2[[2]] *(1/.super$pars$mimics[['ko_k']]) #this is 1/ko_r bc initial km value is in the denomitor of km_cor calc
    }  else {
      Km = .super$pars$km2[[i]]
    }
    Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
    ###MM equation
    C[i] * .super$pars$vmax2[[i]]*10 * C[cat_pool] / (Km_cor*10 + C[i])
    ###
  }
}

# alt mimics function
f_decomp_rmm_twocat <- function(.,C,t,i, cat1 = 3, cat2 = 4) (.super$pars$vmax[[i]]*(C[cat1]+C[cat2])*C[i]) / ((C[cat1]+C[cat2]) + .super$pars$km[[i]])



# MEND-specific functions 
#########

# desorption of C 
f_desorp_mend <- function(., C, t, i ) 
  .super$pars$k[[i]] * C[i] / .super$pars$poolmax[[i]]


# doc uptake 
# APW: some kind of MM
f_docuptake_mend <- function(., C, t, i )
  (.super$pars$vmax[[i]]+.super$pars$mend[['Mr']])*C[i]*C[4] / (C[i]+.super$pars$km[[i]])


# micrbial growth respiration
# APW: some kind of MM
f_growthresp_mend <- function(., C, t, i )
  .super$pars$vmax[[5]]*C[i]*C[5] / (C[5]+.super$pars$km[[5]])


# micrbial maintenance respiration
# APW: some kind of MM
f_maintresp_mend <- function(., C, t, i )
  .super$pars$k[[4]]*C[i]*C[5] / (C[5]+.super$pars$km[[5]])




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


# cat is dummy
# APW: can we use ... ?
#f_desorp_millennialv2_nosat <- function(., C, t, i, k=1, cat=NULL ) 
f_desorp_millennialv2_nosat <- function(., C, t, i, ... ) 
  .super$pars$millennialV2[['kld']] * C[i]



# transfer functions
###################

f_cue_constant <- function(., i ) 
  .super$pars$cue[[i]] 

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

# APW: this is an issue that needs solved
#      we don't need n fluxes functions for each transfer coefficient
#      can something like above be done so it's just one function with n fluxes wrappers
# APW: also this includes a temperature function to calculate CUE, could be a state_pars type calculation   
f_transfer_fluxsum_prop_three_cue <- function(., from, ... ) 
  (.super$pars$cue[[from]] - .super$pars$millennialV2[['cue_t']] * (.super$env$temp - .super$pars$millennialV2[['Taeref']])) * 
  .$transfer_fluxsum_prop(f=3, from=from, ... ) 

# CUE / transfer efficiency sets transfer subject to a maximum pool size 
# - can be used both for saturating MAOM pool and density dependent microbial growth efficiency
# APW: should this include a min 0 function in the final term?
f_transfer_cue_sat <- function(.,C,t,from,to) 
  .super$pars$cue[[from]] * (1-C[to]/.super$pars$poolmax[[to]])

# CENTURY specific CUE calculations
f_transfer_century_quality <- function(.,C,t,from,to) 
  #(1-.super$env$lignin) * .super$pars$century[['strlitter_to_active']] 
  (1-.super$env$lignin) * .super$pars$cue[[from]] 

f_transfer_century_quality2 <- function(.,C,t,from,to) 
  #.super$env$lignin * .super$pars$century[['strlitter_to_slow']] 
  .super$env$lignin * .super$pars$cue2[[from]] 

f_transfer_century_texture <- function(.,C,t,from,to) 
  #(1-f_TEX-.super$pars$century[['active_to_passive']]) 
  (1-.[[paste0('scor.s',from)]]()-.super$pars$cue2[[from]]) 


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



# scaling functions 
# APW: there are more in another file, needs alignment
#############################

f_tcor_none <- function(...) 1
f_wcor_none <- function(...) 1
f_scor_none <- function(...) 1


# correct soil protection rates
f_scor_sulman <- function(.,C,t,i) (.super$env$clay/.super$pars$clayref)^.super$pars$qslope_mayes

# this function not incorporated into the 'water_functions' script because it is not normalized (e.g. 0-1)
# - thus it is kind of integral to the corpse model under the current parameterization and not substitutable
f_wcor_sulman <- function(.,C,t,i) {
  theta <- .super$env$vwc/.super$env$porosity
  theta^3 * (1-theta)^2.5
}

f_tcor_wieder <- function(.,C,t,i) {
  exp(.super$env$temp * .super$pars$mimics[['V_slope']] + .super$pars$mimics[['V_int']]) * .super$pars$mimics[['aV']]
}

# f_tcor_abramoff <- function(.,C,t,i){
#   t1 = 15.4
#   t2 = 11.75
#   t3 = 29.7
#   t4 = 0.031
#   (t2 + (t3/pi)* atan(pi*t4*(.super$env$temp - t1))) / (t2 + (t3/pi)* atan(pi*t4*(.super$pars$reftemp - t1)))
# }

# f_wcor_abramoff <- function(.,C,t,i){
#   w1 = 30
#   w2 = 9
#   (1 / (1+w1*exp(-w2*.super$env$vwc/.35)))
#   #.35 is whc I think? this should be specified in env. Perhaps porosity as in sulman.
# }

# f_tcor_arrhenius <- function(.,C,t,i) {
#   # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
#   # Arrhenius equation
#   
#   # input parameters  
#   # Ea     -- rate of increase to optimum  (J mol-1)
#   # R      -- molar gas constant J mol-1 K-1
#   
#   # Tr     -- reference temperature (oC) 
#   # Trk    -- reference temperature (K) 
#   # Tsk    -- temperature to adjust parameter to (K) 
#   
#   #convert to Kelvin
#   Trk <- .super$pars$reftemp + 273.15
#   Tsk <- .super$env$temp + 273.15
#   
#   exp( .super$pars$ea[[i]]*(Tsk-Trk) / (.super$pars$R*Tsk*Trk) )
# }



### END ###
