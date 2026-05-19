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

calc_state_pars_weider <- function(.) {

  # general parameters
  .super$state_pars$fmet     <- .$fmet()
  .super$state_pars$tau_mod1 <- .$tau_mod1()
  
  # k parameters
  # -- as MEC noted when there is a catalyst pool in mm or rmm decomp, Vmax units are per unit time only
  # -- i.e. they are a specific flux rate per unit mass of the catalyst, so mass units cancel
  # -- thus vmax is equivalent to k  
  for(i in c(5:7)) .super$state_pars$k[[i]] <- .[[paste0('k.k',i)]](i=i) 

  # km parameters
  #for(i in 1:.super$pars$n_outfluxes) .super$state_pars$km[[i]] <- .[[paste0('km.km',i)]](i)
  for(i in c(1:4,8:11)) .super$state_pars$km[[i]] <- .[[paste0('km.km',i)]](i=i) 

  # cue parameters
  for(i in c(1:6,10:11)) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
  for(i in 5:6)          .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](i=i) 
}

calc_state_pars_century <- function(.) {

  # general parameters
  .super$state_pars$matpot_sat <- .$matpot_sat()

  # k parameters
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$k[[i]] <- .[[paste0('k.k',i)]](i=i) 

  # cue parameters
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](.=., i=i ) 
}

calc_state_pars_corpse <- function(.) {

  # k parameters
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$k[[i]]    <- .[[paste0('k.k',i)]](i=i) 
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$km[[i]]   <- .[[paste0('km.km',i)]](i=i) 
  #for(i in 1:.super$pars$n_outfluxes) .super$state_pars$vmax[[i]] <- .[[paste0('vmax.vmax',i)]](i=i) 
  #for(i in 1:.super$pars$n_outfluxes) .super$state_pars$ea[[i]]   <- .[[paste0('ea.ea',i)]](i=i) 

  # cue parameters
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue[[i]]  <- .[[paste0('cue.cue',i)]](i=i) 
  for(i in 1:.super$pars$n_outfluxes) .super$state_pars$cue2[[i]] <- .[[paste0('cue2.cue2',i)]](.=., i=i ) 
}



# k functions
################################

f_k_constant <- function(., i ) 
  .super$pars$k[[i]] 

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
  1 - .$fmet() - .super$pars$cue[[i]] 




### END ###
