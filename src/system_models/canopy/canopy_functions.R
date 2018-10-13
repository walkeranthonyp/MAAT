################################
#
# Canopy level functions
# 
# AWalker December 2015
#
################################



# LAI
################################

f_lai_constant <- function(.) .$pars$lai

f_lai_sphagnum <- function(.) {
  # calculate sphagnum 'lai' as a logistic function of water table depth
  b <- 0.1
    
  .$pars$lai_max*b / (b + exp(.$pars$lai_curve*(.$env$water_td-.$env$sphag_h)/10) )
}


# Light Scaling
################################

# partition direct and diffuse radiation
f_par_partition_spitters <- function(.) {
   #calculate diffuse irradiance (from Spitters et al 1986; and de Jong 1980)

#   if(subd_par) {
     # hourly data
     deJong_R <- 0.847 - 1.61*cos(.$env$zenith) + 1.04*cos(.$env$zenith)^2
     deJong_K <- (1.47-deJong_R)/1.66
     if(.$env$clearness<0.22) {
       diffprop <- 1.0
     } else if(.$env$clearness<0.35) {
       diffprop <- 1.0-6.4*(.$env$clearness-0.22)^2
     } else if(.$env$clearness<deJong_K) {
       diffprop <- 1.47-1.66*.$env$clearness
     } else
       diffprop <- deJong_R
#   } else {
#
#     # daily data
#     if(.$env$clearness<0.07) {
#       diffprop <- 1.0
#     } else if(.$env$clearness<0.35) {
#       diffprop <- 1.0-2.3*(.$env$clearness-0.07)^2
#     } else if(.$env$clearness<0.75) {
#       diffprop <- 1.33-1.46*.$env$clearness
#     } else
#       diffprop <- 0.23
#     }
#   }

   .$env$par_diff[] <- .$env$par*diffprop
   .$env$par_dir[]  <- .$env$par - .$env$par_diff

}


# calculate parameters for light interception
f_pars_init <- function(.) {
  # calculate parameters for canopy light scaling
  
  # extinction coefficents of direct diffuse radiation, assuming leaves are optically black   
  # these account for both canopy clumping and solar zenith angle, also informed by Bodin & Franklin 2012 GMD  
  .$state_pars$G_dir       <- .$pars$G * .$pars$can_clump
  .$state_pars$k_dir       <- .$state_pars$G_dir / cos(.$env$zenith)

  # this is the k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  .$state_pars$k_diff      <- .$state_pars$G_dir / (2/pi)

  # transmittance equals reflectance
  .$state_pars$lscattering <- 2 * .$pars$leaf_reflectance

  # adjustment of k to account for scattering, i.e. that leaves are not optically black
  .$state_pars$m           <- (1.0-.$state_pars$lscattering)^0.5
  .$state_pars$k_dirprime  <- .$state_pars$m * .$state_pars$k_dir
  .$state_pars$k_diffprime <- .$state_pars$m * .$state_pars$k_diff

  # calculate extinction coefficient for Vcmax
  .$state_pars$vcmax_k     <- get(.$fnames$vcmax_k)(.)
}


# Beer's Law
# - these need to be a number of different but similar functions that allow all the commonly used implementations of Beer's Law
# - e.g. accounting for diffuse light, adjusting k by m or not, etc.

# Beer's law direct light
f_rt_beerslaw_goudriaan <- function(.,l) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .$state$vert$sun$apar[]     <- .$state_pars$k_dir * .$env$par * exp(-.$state_pars$k_dir*(l-.$pars$k_layer))
  .$state$vert$sun$fraction[] <- 1.0 
}


# misconstrued Beer's Law - ignores conversion of incident light per unit ground area to incident light per unit leaf area
f_rt_beerslaw <- function(.,l) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .$state$vert$sun$apar[]     <- .$env$par * exp(-.$state_pars$k_dir*(l-.$pars$k_layer))
  .$state$vert$sun$fraction[] <- 1.0 
}


# Goudriaan's Law
f_rt_goudriaan <- function(.,l) {
  # as described in Walker et al. (2017) New Phyt, supplement; from Spitters 1986 and Wang 2003 Func.Plant Biol.     
  # all below values are for visible wavelengths 

  # set leaf absorptance to 1 as this scheme already accounts for leaf scattering
  # as currently implemented this will permanently alter a for a whole ensemble simulation
  .$leaf$pars$a <- 1.0
  
  # calculate albedos
  # after Wang 2003
  alb_h                     <- (1-.$state_pars$m)/(1+.$state_pars$m)
  .$state_pars$alb_dir_can  <- alb_h*2*.$state_pars$k_dir / (.$state_pars$k_dir + .$state_pars$k_diff)
  .$state_pars$alb_diff_can <- 4 * .$state_pars$G_dir * alb_h *
    ( .$state_pars$G_dir * (log(.$state_pars$G_dir)-log(.$state_pars$G_dir+.$state_pars$k_diff)) / .$state_pars$k_diff^2 + 1/.$state_pars$k_diff)

  .$state_pars$alb_dir      <- .$state_pars$alb_dir_can  + (.$pars$alb_soil-.$state_pars$alb_dir_can)  * exp(-2*.$state_pars$k_dirprime  * .$state$lai)
  .$state_pars$alb_diff     <- .$state_pars$alb_diff_can + (.$pars$alb_soil-.$state_pars$alb_diff_can) * exp(-2*.$state_pars$k_diffprime * .$state$lai)
  
  # calculate absorbed direct & diffuse radiation  
  qshade <- (1-.$state_pars$alb_diff)    * .$env$par_diff * .$state_pars$k_diffprime * exp(-.$state_pars$k_diffprime*l) + 
            (1-.$state_pars$alb_dir)     * .$env$par_dir  * .$state_pars$k_dirprime  * exp(-.$state_pars$k_dirprime*l) - 
            (1-.$state_pars$lscattering) * .$env$par_dir  * .$state_pars$k_dir       * exp(-.$state_pars$k_dir*l)
  qshade <- ifelse(qshade<0, 0, qshade)
  
  .$state$vert$shade$apar[] <- qshade
  .$state$vert$sun$apar[]   <- (1-.$state_pars$lscattering)*.$state_pars$k_dir*.$env$par_dir + .$state$vert$shade$apar
  
  # calculate fraction sunlit vs shaded leaves
  .$state$vert$sun$fraction[]   <- exp(-.$state_pars$k_dir*l)
  .$state$vert$shade$fraction[] <- 1 - .$state$vert$sun$fraction

}



# Nitrogen Scaling
################################

f_scale_n_CLMuniform <- function(.,l) {
  rep( (.$state$mass_a/ceiling(.$state$lai)) / .$state$C_to_N, length(l) ) 
}

# Use Beer's Law to scale leaf N through the canopy
f_scale_n_beerslaw <- function(.,l) {
  # for use with a multilayer phototsynthesis scheme
  
  .$state$totalN * exp(-.$state_pars$k_dir*l) /  sum(exp(-.$state_pars$k_dir*1:.$state$lai))  
}

# Use Beer's Law to scale leaf vcmax through the canopy 
f_scale_vcmax_beerslaw <- function(.,l) {
  # for use with a multilayer phototsynthesis scheme
  
  .$state$vcmax0 * exp(-.$state_pars$vcmax_k*l) 
}

f_scale_vcmax_uniform <- function(., layers ) {
  rep(.$state$vcmax0, length(layers) ) 
}

f_vcmax_k_constant  <- function(.) .$pars$vcmax_k

f_vcmax_k_lloyd2012 <- function(.) {
   exp( .$pars$vcmax_k_expa + .$pars$vcmax_k_expb*.$state$vcmax0 )
}


# Canopy Environment 
################################

# Ca
f_scale_ca_uniform <- function(., layers ) {
  rep(.$env$ca_conc, length(layers) ) 
}

# VPD
f_scale_vpd_uniform <- function(., layers ) {
  rep(.$env$vpd, length(layers) ) 
}


# Canopy/Leaf water status  
################################

f_water_status_none <- function(.) .$leaf$state$fwdw_ratio <- NA

# set sphagnum water status
f_water_status_sphagnum <- function(.) {
  .$leaf$state$fwdw_ratio <- get(.$fnames$fwdw)(.) 
}


# Sphagnum functions
################################

# Calculates Sphagnum fresh weight : dry weight ratio as a function of water level (mm) - linear
f_fwdw_wtd_lin <- function(.) {
  
  if((.$env$water_td - .$env$sphag_h) > 0) .$pars$fwdw_wl_sat 
  else .$pars$fwdw_wl_sat + .$pars$fwdw_wl_slope * -(.$env$water_td - .$env$sphag_h) 
}

# Calculates Sphagnum fresh weight : dry weight ratio as a function of water level (mm) - exponential
f_fwdw_wtd_exp <- function(.) {
  # Strack & Price 2009
  
  if((.$env$water_td - .$env$sphag_h) > 0) exp( .$pars$fwdw_wl_exp_b + .$pars$fwdw_wl_exp_a*0 )
  else exp( .$pars$fwdw_wl_exp_b + .$pars$fwdw_wl_exp_a * (-(.$env$water_td - .$env$sphag_h)/10) ) 
}



### END ###
