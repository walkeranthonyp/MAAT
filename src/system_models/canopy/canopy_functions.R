################################
#
# Canopy level functions
# 
# AWalker December 2015
#
################################


# Canopy Environment 
################################

# LAI
f_constant <- function(.) .$pars$lai

f_lai_sphagnum <- function(.) {
  # calculate sphagnum 'lai' as a logistic function of water table height
  b <- 0.1
    
  .$pars$lai_max*b / (b + exp(.$pars$lai_curve*(.$leaf$env$water_l-.$leaf$env$sphag_l)/10) )
}

# Ca
f_Ca_uniform <- function(.,layers,...) {
  rep(.$env$ca_conc,layers) 
}

# VPD
f_vpd_uniform <- function(.,layers,...) {
  rep(.$env$vpd,layers) 
}



# Light Scaling
################################

# calculate parameters for light interception
f_canlight_pars <- function(.){
  # calculate parameters for canopy light scaling
  
  # extinction coefficents of direct  diffuse radiation 
  # these account for both canopy clumping and solar zenith angle, also informed by Bodin & Franklin 2012 GMD  
  .$state_pars$G_dir       <- .$pars$G * .$pars$can_clump
  .$state_pars$k_dir       <- .$state_pars$G_dir / cos(.$env$zenith)
  # this is the k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  .$state_pars$k_diff      <- .$state_pars$G_dir / (2/pi)
  # transmittance equals reflectance
  .$state_pars$lscattering <- 2 * .$pars$leaf_reflectance
  # adjustment of k to account for scattering
  .$state_pars$m           <- (1.0-.$state_pars$lscattering)^0.5
  .$state_pars$k_dirprime  <- .$state_pars$m * .$state_pars$k_dir
  .$state_pars$k_diffprime <- .$state_pars$m * .$state_pars$k_diff
}


# Beer's Law
# - these need to be a number of different but similar functions that allow all the commonly used implementations of Beer's Law
# - e.g. accounting for diffuse light, adjusting k by m or not, e.t.c.
f_canlight_beerslaw_goudriaan <- function(.,l,...){
  # calculates direct beam light attenuation through the canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .$state_pars$k_dir * .$env$par_dir * exp(-.$state_pars$k_dir*(l-.$pars$k_layer))
}

f_canlight_beerslaw <- function(.,l,...){
  # misconstrued Beer's Law - ignores conversion of incident light per unit ground area to incident light per unit leaf area
  
  # calculates direct beam light attenuation through the canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .$env$par_dir * exp(-.$state_pars$k_dir*(l-.$pars$k_layer))
}



# Goudriaan's Law
f_canlight_goudriaan <- function(.,l,...){
  # as described in Walker ... from Spitters 1986 and Wang 2003 Func.Plant Biol.     
  # all below values are for visible wavelengths 
  
  # calculate albedos
  # after Wang 2003
  alb_h                     <- (1-.$state_pars$m)/(1+.$state_pars$m)
  .$state_pars$alb_dir_can  <- alb_h*2*.$state_pars$k_dir / (.$state_pars$k_dir + .$state_pars$k_diff)
  .$state_pars$alb_diff_can <- 4 * .$state_pars$G_dir * alb_h *
    ( .$state_pars$G_dir*(log(.$state_pars$G_dir)-log(.$state_pars$G_dir+.$state_pars$k_diff)) / .$state_pars$k_diff^2 + 1/.$state_pars$k_diff)

  .$state_pars$alb_dir      <- .$state_pars$alb_dir_can  + (.$pars$alb_soil - .$state_pars$alb_dir_can)  * exp(-2*.$state_pars$k_dirprime  * .$state$lai)
  .$state_pars$alb_diff     <- .$state_pars$alb_diff_can + (.$pars$alb_soil - .$state_pars$alb_diff_can) * exp(-2*.$state_pars$k_diffprime * .$state$lai)
  
  
  # calculate absorped direct & diffuse radiation  
  qshade <- (1-.$state_pars$alb_diff)    * .$env$par_diff * .$state_pars$k_diffprime  * exp(-.$state_pars$k_diffprime*l) + 
            (1-.$state_pars$alb_dir)     * .$env$par_dir  * .$state_pars$k_dirprime   * exp(-.$state_pars$k_dirprime*l) - 
            (1-.$state_pars$lscattering) * .$env$par_dir  * .$state_pars$k_dir        * exp(-.$state_pars$k_dir*l)
  if(qshade<0) qshade <- 0
  
  .$state$vert$apar_shade <- qshade
  .$state$vert$apar_sun   <- (1-.$state_pars$lscattering)*.$state_pars$k_dir*.$env$par_dir + .$state$vert$apar_shade
  
  # calculate fraction sunlit vs shaded leaves
  .$state$vert$f_sun      <- exp(-.$state_pars$k_dir*l)
  .$state$vert$f_shade    <- 1 - .state$vert$f_sun

}


# Nitrogen Scaling
################################

f_leafN_CLMuniform <- function(.,layers,...){
  rep(( .$state$mass_a / ceiling(.$state$lai) ) / .$state$C_to_N , layers) 
}


f_leafN_beerslaw <- function(.,l,...){
  # this function uses Beer's Law to scale leaf N through the canopy
  # for use with a multilayer phototsynthesis scheme
  
  .$state$totalN * exp(-.$state_pars$k_dir*l) /  sum(exp(-.$state_pars$k_dir*1:.$state$lai))  
}



### END ###
