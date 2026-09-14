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
  if(!is.na(.super$fnames$calc_state_pars_unique)) .$calc_state_pars_unique()  

  # pool-indexed parameters
  for(p in .super$pool_state_pars)
    for(i in 1:.super$pars$n_pools) {
      .super$state_pars[[p]][[i]]  <- .[[paste0(p,'.',p,i)]](i=i) 
    }

  # outflux-indexed parameters
  for(p in .super$outflux_state_pars)
    for(i in 1:.super$pars$n_outfluxes) {
      # APW: .=. is only needed here for century cue function that calls scor (or any) function, for some reason . doesn't get automatically passed
      # APW: matters even when that scor function does not have ... as an argument 
      #.super$state_pars[[p]][[i]]    <- .[[paste0(p,'.',p,i)]](.=., i=i )
      .super$state_pars[[p]][[i]]    <- .[[paste0(p,'.',p,i)]](i=i)
    }
}


f_calc_state_pars_weider <- function(.) {
  # general parameters
  .super$state_pars$matpot_sat  <- .$matpot_sat()
  .super$state_pars$quality_mod <- .$quality_mod()
  .super$state_pars$anpp_mod    <- .$anpp_mod()
}
  

f_calc_state_pars_century <- function(.) {
  # general parameters
  .super$state_pars$matpot_sat  <- .$matpot_sat()
  .super$state_pars$texture_mod <- .$texture_mod()
}
  


# poolmax functions
################################

f_poolmax_constant <- function(., i ) 
  .super$pars$poolmax[[i]] 

f_poolmax_texture_abramoff <- function(., i )
  #.super$env$BD * .super$env$claysilt * .super$pars$millennialV2[['param_pc']] (in g m-2 version) 
  #(100 - .super$env$sand) * .super$pars$millennialV2[['param_pc']] 
  (1.0 - .super$env$sand) * .super$pars$max_poolmax



# k & Vmax functions
# - when a catalyst is involved in M-M dynamics, Vmax units are per unit time only and is therefore equivalent to k
################################

f_k_constant <- function(., i ) 
  .super$pars$k[[i]] 

f_vmax_constant <- function(., i ) 
  .super$pars$vmax[[i]] 


f_quality_mod_wieder <- function(.) 
  .super$pars$quality_c * (.super$pars$quality_a - .super$pars$quality_b*(.super$env$lignin/.super$env$N))


# ensures that anpp_mod is between two values
# in Will's script, anpp is multipled by 0 in the manipulation script... which would imply that 
# anpp_mod might always be set to 0.6 in the simulations in Ben's paper
f_anpp_mod_wieder <- function(.) 
  .super$pars$anpp_mod_p4 * min(max(sqrt(.super$env$anpp/.super$pars$anpp_mod_p1),.super$pars$anpp_mod_p2),.super$pars$anpp_mod_p3)
 
 
#f_k_wieder <- function(., i )
#  .super$pars$k[[i]]*.super$state_pars$anpp_mod * exp(.super$pars$k_exp[[i]]*.super$state_pars$quality_mod) 
#f_k_wieder_tau_r <- function(., ... )
#  .super$pars$mimics[['tau_r_p1']] * exp(.super$pars$mimics[['tau_r_p2']] * .super$state_pars$fmet) * .super$state_pars$anpp_mod
# 
#f_k_wieder_tau_k <- function(., ... ) 
#  .super$pars$mimics[['tau_k_p1']] * exp(.super$pars$mimics[['tau_k_p2']] * .super$state_pars$fmet) * .super$state_pars$anpp_mod 
 
#f_k_wieder_desorb <- function(., i )
#  #.super$pars$mimics[['desorb_p1']] * exp(.super$pars$mimics[['desorb_p2']] * .super$env$clay) * 0.1
#  .super$pars$k_desorp_tuning*.super$pars$k[[i]] * exp(.super$pars$k_exp[[i]]*.super$env$clay) 


f_k_abramoff_sorp_efficiency <- function(., i ) 
  # APW: note kld is used elsewhere and is kld/1000 currently 
  #exp(-.super$pars$millennialV2[['sorp_p1']] * .super$env$pH - .super$pars$millennialV2[['sorp_p2']]) * .super$pars$millennialV2[['kld']]
  .super$pars$k[[i]] * exp(.super$pars$sorp_exp_a + .super$pars$sorp_exp_b*.super$env$pH) 


f_vmax_wang_cue_inflation <- function(., i )
  1/.super$pars$cue[[i]] * .super$pars$vmax[[i]]  
 

f_vmax_wang_k_from_cat_cue_inflation <- function(., i )
  1/.super$pars$cue[[i]] * .super$pars$k[[.super$pars$cat_pool[[i]] ]]  
 
 

# km functions
################################

f_km_constant <- function(., i ) 
  .super$pars$km[[i]] 

f_km_wieder_temp <- function(., i ) 
  exp(.super$env$temp * .super$pars$km_temp_slope + .super$pars$km_temp_int) * .super$pars$km_tuning /.super$pars$km[[i]]

# Wieder et al. 2015 function for calculating Km using a base Km value and clay content
# -- in MIMICS this function is applied to pool 7 (SOMa) for both types of microbe (r and k)
# APW: the bottom line of this function is the same as above, could combine
f_km_wieder_temp_clay <- function(., i ) { 
  km_texture_scalar <- .super$pars$km_texture_norm * exp(.super$pars$km_texture_exp*sqrt(.super$env$clay))
  km1               <- .super$pars$km[[i]] * km_texture_scalar
  exp(.super$env$temp * .super$pars$km_temp_slope + .super$pars$km_temp_int) * .super$pars$km_tuning / km1
}

##Wieder et al. 2015 function for calculating Km using a base Km value and tuning coefficient (ko_r)
##This function currently specific to r_selected microbes (cat_pool = 3)
##in MIMICS this function is applied to pool 6 (SOMc) for both types of microbe (r and k)
#f_km_wieder_temp_tuning1 <- function(., i ) {
#  km1 = .super$pars$km[[i]] *(1/.super$pars$mimics[['ko_r']]) #
#  exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] / km1
#}
#
##This function currently specific to k_selected microbes (cat_pool = 4)
##in MIMICS this function is applied to pool 6 (SOMc) for both types of microbe (r and k)
#f_km_wieder_temp_tuning2 <- function(., i ) {
#  km1 = .super$pars$km[[i]] *(1/.super$pars$mimics[['ko_k']]) #
#  exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] / km1
#}


 
# CUE functions
################################

f_cue_constant <- function(., i ) 
  .super$pars$cue[[i]] 

f_cue2_constant <- function(., i ) 
  .super$pars$cue2[[i]] 

f_cue3_constant <- function(., i ) 
  .super$pars$cue3[[i]] 


# MIMICS specific CUE calculations
f_cue_wieder_texture <- function(., i )
  .super$pars$cue[[i]]*.super$pars$cue_texture_tuning * exp(.super$pars$cue_texture_exp[[i]]*.super$env$clay)

f_cue_wieder_quality <- function(., i ) 
  .super$pars$cue2[[i]]*.super$pars$cue_quality_tuning * exp(.super$pars$cue_quality_exp[[i]]*.super$state_pars$quality_mod)

#f_cue_wieder_texture1 <- function(., i )
#  # 0.1 manual calibration from Will's script #fSOMp_r 
#  #.super$pars$mimics[['fSOMp_r_p1']] * exp(.super$pars$mimics[['fSOMp_r_p2']]*.super$env$clay) * 0.1 
#  .super$pars$cue[[i]]*.super$pars$cue_texture_tuning * exp(.super$pars$cue_texture_exp*.super$env$clay)
#
#f_cue_wieder_texture2 <- function(., ... ) 
#  # 0.1 manual calibration from Will's script #fSOMp_r 
#  .super$pars$mimics[['fSOMp_k_p1']] * exp(.super$pars$mimics[['fSOMp_k_p2']]*.super$env$clay) * 0.1 
#
#f_cue_wieder_quality1 <- function(., i ) 
#  #.super$pars$mimics[['fSOMc_r_p1']] * exp(.super$pars$mimics[['fSOMc_r_p2']]*.super$state_pars$quality_mod)*.super$pars$mimics[['fSOMc_r_p3']]  #fSOMc_r
#  #.super$pars$cue2[[i]] * exp(.super$pars$fSOMc_r_p2*.super$state_pars$quality_mod)*.super$pars$fSOMc_r_p3 
#  .super$pars$cue2[[i]]*.super$pars$cue_quality_tuning * exp(.super$pars$cue_quality_exp*.super$state_pars$quality_mod)
#
#f_cue_wieder_quality2 <- function(., ... ) 
#  .super$pars$mimics[['fSOMc_k_p1']] * exp(.super$pars$mimics[['fSOMc_k_p2']]*.super$state_pars$quality_mod)*.super$pars$mimics[['fSOMc_k_p3']]  #fSOMc_k


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
  1 - .super$state_pars$texture_mod - .super$pars$cue[[i]] 

f_cue_temp_mod_abramoff  <- function(., i ) 
  #.super$pars$cue[[i]] - .super$pars$millennialV2[['cue_t']] * (.super$env$temp - .super$pars$millennialV2[['reftemp_cue']]) 
  .super$pars$cue[[i]] - .super$pars$cue_temp_scalar*(.super$env$temp-.super$pars$reftemp_cue) 



### END ###
