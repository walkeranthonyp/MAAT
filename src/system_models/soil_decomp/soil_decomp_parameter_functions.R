################################
#
# MAAT soil_decomp process representation functions (PRFs)
# functions to calculate state parameters
#
# Matt Craig, AWalker 2026 
#
################################


# calculation of state parameters functions 
################################

f_calc_state_pars_generic <- function(.) {

  # unique parameters 
  #if(!is.na(.$calc_state_pars_unique)) .$calc_state_pars_unique()  
  if(!is.na(.super$fnames$calc_state_pars_unique)) .$calc_state_pars_unique()  

  # pool-indexed parameters
  #print(.super$pool_state_pars)
  for(p in .super$pool_state_pars)
    for(i in 1:.super$pars$n_pools) {
      #print(paste0(p,'.',p,i)) 
      #print(.[[paste0(p,'.',p,i)]](i=i)) 
      .super$state_pars[[p]][[i]]  <- .[[paste0(p,'.',p,i)]](i=i) 
    }

  # outflux-indexed parameters
  #print(.super$outflux_state_pars)
  for(p in .super$outflux_state_pars)
    for(i in 1:.super$pars$n_outfluxes) {
      #print(paste0(p,'.',p,i)) 
      #print(.[[paste0(p,'.',p,i)]](i=i))
      #print(.[[paste0(p,'.',p,i)]](.=., i=i ))
      # APW: .=. is only needed here for century cue function that calls scor (or any) function, for some reason . doesn't get automatically passed
      # APW: matters even when that scor function does not have ... as an argument 
      .super$state_pars[[p]][[i]]    <- .[[paste0(p,'.',p,i)]](i=i)
      #.super$state_pars[[p]][[i]]    <- .[[paste0(p,'.',p,i)]](.=., i=i )
    }
}


f_calc_state_pars_century_weider <- function(.) {

  # general parameters
  .super$state_pars$matpot_sat <- .$matpot_sat()
  .super$state_pars$fmet       <- .$fmet()
  .super$state_pars$tau_mod1   <- .$tau_mod1()
}
  

#f_calc_state_pars_weider <- function(.) {
#
#  # general parameters
#  .super$state_pars$fmet     <- .$fmet()
#  .super$state_pars$tau_mod1 <- .$tau_mod1()
#  
#  # k parameters
#  # -- as MEC noted when there is a catalyst pool in mm or rmm decomp, Vmax units are per unit time only
#  # -- i.e. they are a specific flux rate per unit mass of the catalyst, so mass units cancel
#  # -- thus vmax is equivalent to k  
#  for(i in c(5:7)) .super$state_pars$k[[i]] <- .[[paste0('k.k',i)]](i=i) 
#
#  # km parameters
#  #for(i in 1:.super$pars$n_outfluxes) .super$state_pars$km[[i]] <- .[[paste0('km.km',i)]](i)
#  for(i in c(1:4,8:11)) .super$state_pars$km[[i]] <- .[[paste0('km.km',i)]](i=i) 
#
#  # cue parameters
#  for(i in c(1:6,10:11)) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
#  for(i in 5:6)          .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](i=i) 
#}
#
#f_calc_state_pars_century <- function(.) {
#
#  # general parameters
#  .super$state_pars$matpot_sat <- .$matpot_sat()
#
#  # k parameters
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$k[[i]] <- .[[paste0('k.k',i)]](i=i) 
#
#  # cue parameters
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](.=., i=i ) 
#}
#
#f_calc_state_pars_corpse <- function(.) {
#
#  # pool-indexed parameters
#  for(i in 1:.super$pars$n_pools) .super$state_pars$poolmax[[i]]  <- .[[paste0('poolmax.poolmax',i)]](i=i) 
#
#  # outflux-indexed parameters
#  # k parameters
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$k[[i]]    <- .[[paste0('k.k',i)]](i=i) 
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$km[[i]]   <- .[[paste0('km.km',i)]](i=i) 
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$vmax[[i]] <- .[[paste0('vmax.vmax',i)]](i=i) 
#  #for(i in 1:.super$pars$n_outfluxes) .super$state_pars$ea[[i]]   <- .[[paste0('ea.ea',i)]](i=i) 
#
#  # cue parameters
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](.=., i=i ) 
#  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue3[[i]] <- .[[paste0('cue3.cue3',i)]](.=., i=i ) 
#}



# poolmax functions
################################

f_poolmax_constant <- function(., i ) 
  .super$pars$poolmax[[i]] 

f_poolmax_texture_abramoff <- function(., i )
  #.super$env$BD * .super$env$claysilt * .super$pars$millennialV2[['param_pc']] (in g m-2 version) 
  (100 - .super$env$sand) * .super$pars$millennialV2[['param_pc']] 



# k & Vmax functions
# - when a catalyst is involved in M-M dynamics, Vmax units are per unit time only and is therefore equivalent to k
################################

f_k_constant <- function(., i ) 
  .super$pars$k[[i]] 

f_vmax_constant <- function(., i ) 
  .super$pars$vmax[[i]] 

f_texture_wieder_fmet <- function(.) 
  .super$pars$mimics[['fmet_p1']] * (.super$pars$mimics[['fmet_p2']] - .super$pars$mimics[['fmet_p3']]*(.super$env$lignin/.super$env$N))

# ensures that tau_mod1 is between two values
# in Will's script, anpp is multipled by 0 in the manipulation script... which would imply that 
# tau_mod1 might always be set to 0.6 in the simulations in Ben's paper
f_k_wieder_tau_mod1 <- function(.) 
  min(max(sqrt(.super$env$anpp/.super$pars$mimics[['tau_mod1_p1']]),.super$pars$mimics[['tau_mod1_p2']]),.super$pars$mimics[['tau_mod1_p3']])
 
f_k_wieder_tau_r <- function(., ... )
  .super$pars$mimics[['tau_r_p1']] * exp(.super$pars$mimics[['tau_r_p2']] * .super$state_pars$fmet) * .super$state_pars$tau_mod1 * .super$pars$mimics[['tau_mod2']]
 
f_k_wieder_tau_k <- function(., ... ) 
  .super$pars$mimics[['tau_k_p1']] * exp(.super$pars$mimics[['tau_k_p2']] * .super$state_pars$fmet) * .super$state_pars$tau_mod1 * .super$pars$mimics[['tau_mod2']]
 
f_k_wieder_desorb <- function(., ... ) 
  .super$pars$mimics[['desorb_p1']] * exp(.super$pars$mimics[['desorb_p2']] * .super$env$clay) * 0.1

f_vmax_wang_cue_inflation <- function(., i )
  1/.super$pars$cue[[i]] * .super$pars$vmax[[i]]  
 
f_vmax_wang_k_from_cat_cue_inflation <- function(., i )
  1/.super$pars$cue[[i]] * .super$pars$k[[.super$pars$cat_pool[[i]] ]]  
 
f_k_abramoff_sorp_efficiency <- function(., i )
  # kaff_lm
  # APW: note kld is used elsewhere and is kld/1000 currently 
  exp(-.super$pars$millennialV2[['sorp_p1']] * .super$env$pH - .super$pars$millennialV2[['sorp_p2']]) * .super$pars$millennialV2[['kld']]



# km functions
################################

f_km_constant <- function(., i ) 
  .super$pars$km[[i]] 

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


 
# CUE functions
################################

f_cue_constant <- function(., i ) 
  .super$pars$cue[[i]] 

f_cue2_constant <- function(., i ) 
  .super$pars$cue2[[i]] 

f_cue3_constant <- function(., i ) 
  .super$pars$cue3[[i]] 


# MIMICS specific CUE calculations
f_cue_wieder_texture1 <- function(., ... ) 
  .super$pars$mimics[['fSOMp_r_p1']] * exp(.super$pars$mimics[['fSOMp_r_p2']]*.super$env$clay) * 0.1 #0.1 manual calibration from Will's script #fSOMp_r

f_cue_wieder_quality1 <- function(., ... ) 
  .super$pars$mimics[['fSOMc_r_p1']] * exp(.super$pars$mimics[['fSOMc_r_p2']]*.super$state_pars$fmet)*.super$pars$mimics[['fSOMc_r_p3']]  #fSOMc_r

f_cue_wieder_texture2 <- function(., ... ) 
  .super$pars$mimics[['fSOMp_k_p1']] * exp(.super$pars$mimics[['fSOMp_k_p2']]*.super$env$clay) * 0.1 #fSOMp_k

f_cue_wieder_quality2 <- function(., ... ) 
  .super$pars$mimics[['fSOMc_k_p1']] * exp(.super$pars$mimics[['fSOMc_k_p2']]*.super$state_pars$fmet)*.super$pars$mimics[['fSOMc_k_p3']]  #fSOMc_k


# CENTURY specific CUE calculations
f_cue_century_quality <- function(., i )
  (1-.super$env$lignin) * .super$pars$cue[[i]] 

f_cue_century_quality2 <- function(., i )
  .super$env$lignin * .super$pars$cue2[[i]] 

f_cue_century_texture <- function(., i ) 
  1 - .[[paste0('scor.scor',i)]]() - .super$pars$cue[[i]] 

# APW: this is a the same as above but calls the texture function differently
f_cue_century_texture_elm <- function(., i ) 
  #1 - .$fmet() - .super$pars$cue[[i]] 
  1 - .super$state_pars$fmet - .super$pars$cue[[i]] 

f_cue_temp_mod_abramoff  <- function(., i ) 
  .super$pars$cue[[i]] - .super$pars$millennialV2[['cue_t']] * (.super$env$temp - .super$pars$millennialV2[['reftemp_cue']]) 



### END ###
