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
f_canlight_beerslaw <- function(.,l,...){
  # calculates direct beam light attenuation through the canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .$state_pars$k_dir * .$env$par_dir * exp(-.$state_pars$k_dir*(l-.$pars$k_layer))
}

f_canlight_beerslaw_wrong <- function(.,l,...){
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





# run functions
###############################

# Big Leaf canopy scaling 
# Sellers et al 1992
f_bigleaf_sellers1992 <- function(.,k=.$state_pars$k_dirprime,...) {
  # this function was described in Sellers to deal with time intergrated values of fpar and k 
  # could write wrapper function or if to pass different k coefficients
  
  # calculate fPAR, after Sellers 1992, see also Friend 2001
  fpar <- 1 - exp(-k*.$state$lai)
  
  # set leaf environment
  # I0 is incident light
  .$leaf$env$par     <- .$k_dir * (1-.$leafscattering) * .$env$par
  # assume no variation in CO2 concentration
  .$leaf$env$ca_conc <- .$env$ca_conc
  # VPD
  .$leaf$env$vpd     <- .$env$vpd
  
  # set leaf N0 
  .$leaf$state$leafN_area <- .$state$totalN * k / fpar
  
  # calculate A0
  # leaf model will calculate Vcmax0 and Jmax0 according to N0, temp, etc 
  .$leaf$run()
  
  # scale
  .$state$integrated$A             <- .$leaf$state$A * fpar/k
  .$state$integrated$respiration   <- .$leaf$state$respiration * fpar/k
  .$state$integrated$wc_lim        <- if(.$leaf$state$lim=='wc') .$state$integrated$A else 0
  .$state$integrated$wj_lim        <- if(.$leaf$state$lim=='wj') .$state$integrated$A else 0 
  .$state$integrated$wp_lim        <- if(.$leaf$state$lim=='wp') .$state$integrated$A else 0
  .$state$integrated$layers_wc_lim <- if(.$leaf$state$lim=='wc') .$state$lai
  .$state$integrated$layers_wj_lim <- if(.$leaf$state$lim=='wj') .$state$lai
  .$state$integrated$layers_wp_lim <- if(.$leaf$state$lim=='wp') .$state$lai
  # resistance - convert to conductance, minus minimum conductance, scale, add min conductance multiplied by LAI, convert back to resistance
  # - a somewhat complicated version of eq 37f & 35 from Sellers 1992 and a conductance to resistance conversion
  .$state$integrated$rs            <- 1 / ( (1/.$leaf$state$rs - .$leaf$pars$g0) * fpar/k + .$leaf$pars$g0*.$state$lai )
#   .$state$integrated$gi            <- sum(.$state$gi)
  # canopy mean values of Ci and Cc
  .$state$integrated$ci            <- .$leaf$state$ci * fpar/k / .$state$lai
  .$state$integrated$cc            <- .$leaf$state$cc * fpar/k / .$state$lai
  
}


# Two Big Leaf canopy scaling 
# - accounts for direct and diffuse light separately but only calculates A once for each stream
f_2bigleaf <- function(.) {
  # Thornton calculates the mean of all these canopy values, what does Dai do? 
  
  # calculate LAIsun and LAIshade - Dai 2004
  Lsun   <- (1 - exp(-.$state_pars$k_dir*.$state$lai)) / .$state_pars$k_dir
  Lshade <- .$state$lai - Lsun 
  
  # calculate APARsun and APARshade  
  # APARshade is the scattered part of the direct beam plus diffuse radiation
  
  # APARsun is the direct beam plus the scattered part of the direct beam plus diffuse radiation
  
  
  # calculate Nsun and Nshade 
  
  # calculate Asun and Ashade
  
  # scale
}


# Multilayer canopy scaling
f_multilayer <- function(.) {
  
  # initialise number of layers
  # explicitely scale canopy N
  # explicitely scale PARsun and PARshade
  # explicitely calculate Asun and Ashade in each layer
  # scale
  
  # initialise layers
  layers      <- ceiling(.$state$lai) # this could be a function specifying either no. of layers or the below lai assignment   
  
  # canopy layer 
  # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
  # populate layer dataframe      
  .$state$vert$leaf.ca_conc    <- get(.$fnames$can_scale_Ca)(.,1:layers)
  .$state$vert$leaf.vpd        <- get(.$fnames$can_scale_vpd)(.,1:layers)
  .$state$vert$leaf.par        <- get(.$fnames$can_scale_light)(.,1:layers)
  .$state$vert$leaf.leafN_area <- get(.$fnames$can_scale_N)(.,l=1:layers,layers=layers)
  
  # run leaf
  lout   <- as.data.frame(do.call(rbind,lapply(1:layers,.$run_leaf)))
  #       lout   <- mclapply(1:layers,.$leaf$wrapper,df=leafdf,mc.cores=.$pars$leaf_cores) # this can be used to run each canopy layer in parrallel, in practise this was found to be slower
  
  # assign data
  .$state$A           <- unlist(lout$A)
  .$state$cc          <- unlist(lout$cc)
  .$state$ci          <- unlist(lout$ci)
  .$state$ri          <- unlist(lout$ri)
  .$state$rs          <- unlist(lout$rs)
  .$state$lim         <- unlist(lout$lim)
  .$state$respiration <- unlist(lout$respiration)      
  
  #scale bottom layer fluxes/pars to leaf area in that layer
  # - need to generalise this to a multilayer canopy
  if(ceiling(.$state$lai)-floor(.$state$lai)==1){
    scale <- .$state$lai %% 1
    .$state$A[layers]           <- .$state$A[layers]           * scale        
    .$state$ri[layers]          <- .$state$ri[layers]          * scale        
    .$state$rs[layers]          <- .$state$rs[layers]          * scale        
    .$state$respiration[layers] <- .$state$respiration[layers] * scale      
    .$state$cc[layers]          <- .$state$cc[layers]          * scale        
    .$state$ci[layers]          <- .$state$ci[layers]          * scale        
  }      
  
  #integrate canopy layers
  # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
  # canopy sum values
  .$state$integrated$A             <- sum(.$state$A)
  .$state$integrated$respiration   <- sum(.$state$respiration)
  .$state$integrated$wc_lim        <- sum(.$state$A * (.$state$lim=='wc')) 
  .$state$integrated$wj_lim        <- sum(.$state$A * (.$state$lim=='wj'))
  .$state$integrated$wp_lim        <- sum(.$state$A * (.$state$lim=='wp'))
  .$state$integrated$layers_wc_lim <- sum(.$state$lim=='wc')
  .$state$integrated$layers_wj_lim <- sum(.$state$lim=='wj')
  .$state$integrated$layers_wp_lim <- sum(.$state$lim=='wp')
  .$state$integrated$ri            <- 1 / sum(1/.$state$ri)
  .$state$integrated$rs            <- 1 / sum(1/.$state$rs)
  # canopy mean values
  .$state$integrated$cc            <- sum(.$state$cc) / .$state$lai
  .$state$integrated$ci            <- sum(.$state$ci) / .$state$lai
  
}






















