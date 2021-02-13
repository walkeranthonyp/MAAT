################################
#
# Canopy level functions
# 
# AWalker December 2015
#
################################



# Light Scaling
################################

# partitioned par passed as an environmental variable
f_par_partition_env <- function(.) print('/nExpects direct and diffuse PAR passed as environment variables/n')

# partition direct and diffuse radiation
f_par_partition_spitters_hourly <- function(.) {
   #calculate diffuse irradiance (from Spitters et al 1986; and de Jong 1980)

   # hourly data
   deJong_R <- 0.847 - 1.61*cos(.super$env$zenith) + 1.04*cos(.super$env$zenith)^2
   deJong_K <- (1.47-deJong_R)/1.66
   if(.super$env$clearness<0.22) {
     diffprop <- 1.0
   } else if(.super$env$clearness<0.35) {
     diffprop <- 1.0-6.4*(.super$env$clearness-0.22)^2
   } else if(.super$env$clearness<deJong_K) {
     diffprop <- 1.47-1.66*.super$env$clearness
   } else {
     diffprop <- deJong_R
   }

   .super$env$par_diff[] <- .super$env$par*diffprop
   .super$env$par_dir[]  <- .super$env$par - .super$env$par_diff
}

f_par_partition_spitters_daily <- function(.) {
   #calculate diffuse irradiance (from Spitters et al 1986; and de Jong 1980)

   # daily data
   if(.super$env$clearness<0.07) {
     diffprop <- 1.0
   } else if(.super$env$clearness<0.35) {
     diffprop <- 1.0-2.3*(.super$env$clearness-0.07)^2
   } else if(.super$env$clearness<0.75) {
     diffprop <- 1.33-1.46*.super$env$clearness
   } else {
     diffprop <- 0.23
   }

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

  # albedo
  .$albedo()

  # calculate Vcmax0 and extinction coefficient for Vcmax
  .super$state$vcmax0[]           <- .$fns$vcmax0()
  .super$state_pars$k_vcmax[]     <- .$fns$k_vcmax()
}

f_pars_init_full <- function(.) {
  # calculate parameters for canopy light scaling
 
  # zenith angles (radians) for numerical approximation of diffuse light extinction & albedo
  dz <- 1/.super$pars$nz
  .super$state_pars$zi       <- zi <- seq(dz/2,1,dz) * pi/2
  .super$state_pars$delta_zi <- delta_zi <- sum(.super$state_pars$zi[1:2] * c(-1,1))  
  
  # extinction coefficents of direct & diffuse radiation, assuming leaves are optically black   
  .super$state_pars$G_dir[]       <- .$gz(z=.super$env$zenith) 
  .super$state_pars$k_dir[]       <- .super$state_pars$G_dir / cos(.super$env$zenith)
  .super$state_pars$k_dir_zi      <- .$gz(z=zi) / cos(zi)
  # k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  numint                          <- sin(zi)*cos(zi)*delta_zi * exp(-.super$state_pars$k_dir_zi * .super$env$lai) 
  .super$state_pars$k_diff[]      <- -log(2*sum(numint)) / .super$env$lai 

  # scattering transmittance + reflectance
  .super$state_pars$lscattering[] <- .super$pars$leaf_reflectance + .super$pars$leaf_transmitance
  .super$state_pars$m[]           <- (1.0-.super$state_pars$lscattering)^0.5
  .super$state_pars$k_dirprime[]  <- .super$state_pars$m * .super$state_pars$k_dir
  .super$state_pars$k_diffprime[] <- .super$state_pars$m * .super$state_pars$k_diff

  # adjustment of k_Xprime 
  # - to account for clumping,  
  .super$state_pars$k_dirprime[]  <- .super$state_pars$k_dirprime[]  * .super$pars$can_clump
  .super$state_pars$k_diffprime[] <- .super$state_pars$k_diffprime[] * .super$pars$can_clump

  # albedo
  .$albedo()

  # adjustment of k_X 
  # - to account for clumping,  
  .super$state_pars$k_dir[]       <- .super$state_pars$k_dir[]  * .super$pars$can_clump
  .super$state_pars$k_diff[]      <- .super$state_pars$k_diff[] * .super$pars$can_clump

  # calculate Vcmax0 and extinction coefficient for Vcmax
  .super$state$vcmax0[]           <- .$fns$vcmax0()
  .super$state_pars$k_vcmax[]     <- .$fns$k_vcmax()
}

# calculate G(z), Bonan 2019
f_gz_rossgoudriaan <- function(., chi_l=.super$pars$chi_l, z=.super$env$zenith) {
  # Ross-Goudriaan function to approxiamte G(z) from 
  # Ross index (chi_l) which indicates deviation from a spherical leaf angle distribution   
  
  # this semi-empirical function not applicable outside range 
  if(chi_l>0.6|chi_l<(-0.4)){
    print('chi_l outside interval: -0.4 to 0.6, changed to boundary')
    chi_l = min(max(chi_l, -0.4), 0.6)
  }
  
  phi1 <- 0.5 - 0.633*chi_l - 0.33*chi_l^2
  phi2 <- 0.877*(1 - 2*phi1)
  phi1 + phi2*cos(z)
}


# Beer's Law
# - these need to be a number of different but similar functions that allow all the commonly used implementations of Beer's Law
# - e.g. accounting for diffuse light, adjusting k by m or not, etc.

# Beer's law - from Sellers (1992) and citing Goudriaan (1977) 
f_rt_beerslaw_goudriaan <- function(.,l) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .super$state$vert$sun$apar[]     <- (1-.super$state_pars$lscattering)*.super$state_pars$k_dir*.super$env$par * exp(-.super$state_pars$k_dirprime*l)
  .super$state$vert$sun$fraction[] <- 1.0 
}


# misconstrued Beer's Law
# - ignores conversion of incident light per unit ground area to incident light per unit leaf area
# - considers scattering in absorption because the leaf models does not when nested in canopy
# - but scattering is not considered in RT, k_dir not k_dirprime is used
f_rt_beerslaw <- function(.,l) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .super$state$vert$sun$apar[]     <- (1-.super$state_pars$lscattering)*.super$env$par * exp(-.super$state_pars$k_dir*l)
  .super$state$vert$sun$fraction[] <- 1.0 
}


# Goudriaan's Law
f_rt_goudriaan <- function(.,l) {
  # as described in Walker et al. (2017) New Phyt, supplement; from Spitters 1986 and Wang 2003 Func.Plant Biol.     
  # all below values are for visible wavelengths 

  # calculate albedos
  #.$albedo()
  
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


# albedo functions for Goudriaan's model
f_albedo_goudriaan <- function(.) {
  # after Wang 2003, but after rereading looks different
  # informed by Bonan 2019, Climate Change and Terrestrial Ecosystem Modeling
  
  # albedo for canopies of infinite LAI
  # with horizontal leaves
  .super$state_pars$alb_h         <- (1-.super$state_pars$m) / (1+.super$state_pars$m)
  # adjust direct beam albedo for leaf angle distribution
  .super$state_pars$alb_dir_can  <- .super$state_pars$alb_h*2*.super$state_pars$k_dir / (.super$state_pars$k_dir + .super$state_pars$k_diff)
  # adjust diffuse beam albedo for leaf angle distribution
  .$diffalbedo()
  print(.super$state_pars$alb_dir_can)
  print(.super$state_pars$alb_diff_can)
#  .super$state_pars$alb_diff_can <- 4 * .super$state_pars$G_dir * alb_h * ( .super$state_pars$G_dir *
#                                      (log(.super$state_pars$G_dir)-log(.super$state_pars$G_dir+.super$state_pars$k_diff)) / .super$state_pars$k_diff^2 + 
#                                      1/.super$state_pars$k_diff )

  # adjust for finite LAI
  .super$state_pars$alb_dir      <- .super$state_pars$alb_dir_can  + (.super$pars$alb_soil-.super$state_pars$alb_dir_can)  * 
                                       exp(-2*.super$state_pars$k_dirprime  * .super$env$lai)
  .super$state_pars$alb_diff     <- .super$state_pars$alb_diff_can + (.super$pars$alb_soil-.super$state_pars$alb_diff_can) * 
                                       exp(-2*.super$state_pars$k_diffprime * .super$env$lai)
  print(.super$state_pars$alb_dir)
  print(.super$state_pars$alb_diff)
}


f_diffalbedo_goudriaan <- function(., zi=.$state_pars$zi, delta_zi=.$state_pars$delta_zi ) {
  
  numint <- sin(zi)*cos(zi)*delta_zi * .super$state_pars$alb_h*.super$state_pars$k_dir_zi / (.super$state_pars$k_dir_zi + .super$state_pars$k_diff) 
  .super$state_pars$alb_diff_can  <- 2 * sum(numint) 
} 

f_diffalbedo_approx <- function(.) {
  .super$state_pars$alb_diff_can <- 4 * .super$state_pars$G_dir * .super$state_pars$alb_h * ( .super$state_pars$G_dir *
                                      (log(.super$state_pars$G_dir)-log(.super$state_pars$G_dir+.super$state_pars$k_diff)) / .super$state_pars$k_diff^2 + 
                                      1/.super$state_pars$k_diff )
} 


# Nitrogen Scaling
################################

f_scale_n_CLMuniform <- function(.,l) {
  rep( (.super$state$mass_a/ceiling(.super$env$lai)) / .super$state$C_to_N, length(l) ) 
}

# Use Beer's Law to scale leaf N through the canopy
f_scale_n_beerslaw <- function(.,l) {
  # for use with a multilayer phototsynthesis scheme
  
  .super$state$totalN * exp(-.super$state_pars$k_dir*l) /  sum(exp(-.super$state_pars$k_dir*1:.super$env$lai))  
}

## Use Beer's Law to scale leaf vcmax through the canopy 
#f_scale_vcmax_beerslaw <- function(.,l) {
#  # for use with a multilayer phototsynthesis scheme
#  
#  .$state$vcmax0 * exp(-.$state_pars$k_vcmax*l) 
#}
#
#f_scale_vcmax_uniform <- function(., layers ) {
#  rep(.$state$vcmax0, length(layers) ) 
#}
#
#f_k_vcmax_constant  <- function(.) .$pars$k_vcmax
#
#f_k_vcmax_lloyd2012 <- function(.) {
#   exp( .$pars$k_vcmax_expa + .$pars$k_vcmax_expb*.$state$vcmax0 )
#}



# Vcmax Scaling
################################

f_vcmax0_constant <- function(.) {
  #.super$leaf$fnames$vcmax <- 'f_vcmax_constant'
  .super$configure(.=.super, vlist='fnames', df=c(leaf.vcmax='f_vcmax_constant') )
  .super$pars$vcmax0
}


f_k_vcmax_constant  <- function(.) .super$pars$k_vcmax


f_k_vcmax_lloyd2012 <- function(.) {
   exp( .super$pars$k_vcmax_expa + .super$pars$k_vcmax_expb*.super$state$vcmax0 )
}


# Use Beer's Law to scale leaf vcmax through the canopy 
f_scale_vcmax_beerslaw <- function(., l, var ) {
  # for use with a multilayer phototsynthesis scheme
  # currently var is a dummy argument for compatibility with two_layer
  
  .super$state$vcmax0 * exp(-.super$state_pars$k_vcmax*l) 
}


f_scale_vcmax_uniform <- function(., layers, var ) {
  # currently var is a dummy argument for compatibility with two_layer
  rep(.super$state$vcmax0, length(layers) ) 
}


f_scale_two_layer <- function(., l, var ) {
  can_vector <- l
  can_vector[] <- .super$pars[[var]]$layer2
  can_vector[l<.super$pars$can_lai_upper] <- .super$pars[[var]]$layer1
  can_vector
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



### END ###
