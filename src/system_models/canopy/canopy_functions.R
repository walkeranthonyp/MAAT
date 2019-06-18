################################
#
# Canopy level functions
# 
# AWalker December 2015
#
################################



# LAI
################################

f_lai_constant <- function(.) .super$pars$lai

f_lai_sphagnum <- function(.) {
  # calculate sphagnum 'lai' as a logistic function of water table depth
  b <- 0.1
    
  .super$pars$lai_max*b / (b + exp(.super$pars$lai_curve*(.super$env$water_td-.super$env$sphag_h)/10) )
}


# Light Scaling
################################

# partitioned par passed as an environmental variable
f_par_partition_env <- function(.) print('/nExpects direct and diffuse PAR passed as environment variables/n')

# partition direct and diffuse radiation
f_par_partition_spitters <- function(.) {
   #calculate diffuse irradiance (from Spitters et al 1986; and de Jong 1980)

#   if(subd_par) {
     # hourly data
     deJong_R <- 0.847 - 1.61*cos(.super$env$zenith) + 1.04*cos(.super$env$zenith)^2
     deJong_K <- (1.47-deJong_R)/1.66
     if(.super$env$clearness<0.22) {
       diffprop <- 1.0
     } else if(.super$env$clearness<0.35) {
       diffprop <- 1.0-6.4*(.super$env$clearness-0.22)^2
     } else if(.super$env$clearness<deJong_K) {
       diffprop <- 1.47-1.66*.super$env$clearness
     } else
       diffprop <- deJong_R
#   } else {
#
#     # daily data
#     if(.super$env$clearness<0.07) {
#       diffprop <- 1.0
#     } else if(.super$env$clearness<0.35) {
#       diffprop <- 1.0-2.3*(.super$env$clearness-0.07)^2
#     } else if(.super$env$clearness<0.75) {
#       diffprop <- 1.33-1.46*.super$env$clearness
#     } else
#       diffprop <- 0.23
#     }
#   }

   .super$env$par_diff[] <- .super$env$par*diffprop
   .super$env$par_dir[]  <- .super$env$par - .super$env$par_diff

}


# calculate parameters for light interception
f_pars_init <- function(.) {
  # calculate parameters for canopy light scaling
  
  # extinction coefficents of direct diffuse radiation, assuming leaves are optically black   
  # these account for both canopy clumping and solar zenith angle, also informed by Bodin & Franklin 2012 GMD  
  .super$state_pars$G_dir[]       <- .super$pars$G * .super$pars$can_clump
  .super$state_pars$k_dir[]       <- .super$state_pars$G_dir / cos(.super$env$zenith)

  # this is the k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  .super$state_pars$k_diff[]      <- .super$state_pars$G_dir / (2/pi)

  # transmittance equals reflectance
  .super$state_pars$lscattering[] <- 2 * .super$pars$leaf_reflectance

  # adjustment of k to account for scattering, i.e. that leaves are not optically black
  .super$state_pars$m[]           <- (1.0-.super$state_pars$lscattering)^0.5
  .super$state_pars$k_dirprime[]  <- .super$state_pars$m * .super$state_pars$k_dir
  .super$state_pars$k_diffprime[] <- .super$state_pars$m * .super$state_pars$k_diff

  # calculate Vcmax0 and extinction coefficient for Vcmax
  .super$state$vcmax0[]           <- get(.$fnames$vcmax0)(.)
  .super$state_pars$k_vcmax[]     <- get(.$fnames$k_vcmax)(.)
}


# Beer's Law
# - these need to be a number of different but similar functions that allow all the commonly used implementations of Beer's Law
# - e.g. accounting for diffuse light, adjusting k by m or not, etc.

# Beer's law - from Sellers (1992) and citing Goudriaan (1977) 
f_rt_beerslaw_goudriaan <- function(.,l) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .super$state$vert$sun$apar[]     <- (1 - .super$state_pars$lscattering) * .super$state_pars$k_dir * .super$env$par * exp(-.super$state_pars$k_dirprime*l)
  .super$state$vert$sun$fraction[] <- 1.0 
}


# misconstrued Beer's Law - ignores conversion of incident light per unit ground area to incident light per unit leaf area
f_rt_beerslaw <- function(.,l) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .super$state$vert$sun$apar[]     <- (1 - .super$state_pars$lscattering) * .super$env$par * exp(-.super$state_pars$k_dir*l)
  .super$state$vert$sun$fraction[] <- 1.0 
}


# Goudriaan's Law
f_rt_goudriaan <- function(.,l) {
  # as described in Walker et al. (2017) New Phyt, supplement; from Spitters 1986 and Wang 2003 Func.Plant Biol.     
  # all below values are for visible wavelengths 

  # calculate albedos
  # after Wang 2003
  alb_h                     <- (1-.super$state_pars$m)/(1+.super$state_pars$m)
  .super$state_pars$alb_dir_can  <- alb_h*2*.super$state_pars$k_dir / (.super$state_pars$k_dir + .super$state_pars$k_diff)
  .super$state_pars$alb_diff_can <- 4 * .super$state_pars$G_dir * alb_h *
    ( .super$state_pars$G_dir * (log(.super$state_pars$G_dir)-log(.super$state_pars$G_dir+.super$state_pars$k_diff)) / .super$state_pars$k_diff^2 + 1/.super$state_pars$k_diff)

  .super$state_pars$alb_dir      <- .super$state_pars$alb_dir_can  + (.super$pars$alb_soil-.super$state_pars$alb_dir_can)  * exp(-2*.super$state_pars$k_dirprime  * .super$state$lai)
  .super$state_pars$alb_diff     <- .super$state_pars$alb_diff_can + (.super$pars$alb_soil-.super$state_pars$alb_diff_can) * exp(-2*.super$state_pars$k_diffprime * .super$state$lai)
  
  # calculate absorbed direct & diffuse radiation  
  qshade <- (1-.super$state_pars$alb_diff)    * .super$env$par_diff * .super$state_pars$k_diffprime * exp(-.super$state_pars$k_diffprime*l) + 
            (1-.super$state_pars$alb_dir)     * .super$env$par_dir  * .super$state_pars$k_dirprime  * exp(-.super$state_pars$k_dirprime*l) - 
            (1-.super$state_pars$lscattering) * .super$env$par_dir  * .super$state_pars$k_dir       * exp(-.super$state_pars$k_dir*l)
  qshade <- ifelse(qshade<0, 0, qshade)
  
  .super$state$vert$shade$apar[] <- qshade
  .super$state$vert$sun$apar[]   <- (1-.super$state_pars$lscattering)*.super$state_pars$k_dir*.super$env$par_dir + .super$state$vert$shade$apar
  
  # calculate fraction sunlit vs shaded leaves
  .super$state$vert$sun$fraction[]   <- exp(-.super$state_pars$k_dir*l)
  .super$state$vert$shade$fraction[] <- 1 - .super$state$vert$sun$fraction

}



# Nitrogen Scaling
################################

f_scale_n_CLMuniform <- function(.,l) {
  rep( (.super$state$mass_a/ceiling(.super$state$lai)) / .super$state$C_to_N, length(l) ) 
}

# Use Beer's Law to scale leaf N through the canopy
f_scale_n_beerslaw <- function(.,l) {
  # for use with a multilayer phototsynthesis scheme
  
  .super$state$totalN * exp(-.super$state_pars$k_dir*l) /  sum(exp(-.super$state_pars$k_dir*1:.super$state$lai))  
}


# Vcmax Scaling
################################

f_vcmax0_constant <- function(.) {
  .super$leaf$fnames$vcmax <- 'f_vcmax_constant'
  .super$pars$vcmax0
}

# Use Beer's Law to scale leaf vcmax through the canopy 
f_scale_vcmax_beerslaw <- function(.,l) {
  # for use with a multilayer phototsynthesis scheme
  
  .super$state$vcmax0 * exp(-.super$state_pars$k_vcmax*l) 
}

f_scale_vcmax_uniform <- function(., layers ) {
  rep(.super$state$vcmax0, length(layers) ) 
}

f_k_vcmax_constant  <- function(.) .super$pars$k_vcmax

f_k_vcmax_lloyd2012 <- function(.) {
   exp( .super$pars$k_vcmax_expa + .super$pars$k_vcmax_expb*.super$state$vcmax0 )
}


# Canopy Environment 
################################

# Ca
f_scale_ca_uniform <- function(., layers ) {
  rep(.super$env$ca_conc, length(layers) ) 
}

# VPD
f_scale_vpd_uniform <- function(., layers ) {
  rep(.super$env$vpd, length(layers) ) 
}


# Canopy/Leaf water status  
################################

f_water_status_none <- function(.) .super$leaf$state$fwdw_ratio <- NA

# set sphagnum water status
f_water_status_sphagnum <- function(.) {
  .super$leaf$state$fwdw_ratio <- get(.$fnames$fwdw)(.) 
}


# Sphagnum functions
################################

# Calculates Sphagnum fresh weight : dry weight ratio as a function of water level (mm) - linear
f_fwdw_wtd_lin <- function(.) {
  
  if((.super$env$water_td - .super$env$sphag_h) > 0) .super$pars$fwdw_wl_sat 
  else .super$pars$fwdw_wl_sat + .super$pars$fwdw_wl_slope * -(.super$env$water_td - .super$env$sphag_h) 
}

# Calculates Sphagnum fresh weight : dry weight ratio as a function of water level (mm) - exponential
f_fwdw_wtd_exp <- function(.) {
  # Strack & Price 2009
  
  if((.super$env$water_td - .super$env$sphag_h) > 0) exp( .super$pars$fwdw_wl_exp_b + .super$pars$fwdw_wl_exp_a*0 )
  else exp( .super$pars$fwdw_wl_exp_b + .super$pars$fwdw_wl_exp_a * (-(.super$env$water_td - .super$env$sphag_h)/10) ) 
}



### END ###
