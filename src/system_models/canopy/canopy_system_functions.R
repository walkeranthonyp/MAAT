################################
#
# Canopy system functions 
# 
# AWalker February 2018
#
################################



# Big Leaf canopy scaling 
###############################
# Sellers et al 1992
f_cansys_bigleaf_s1992 <- function(.,k=.$state_pars$k_dirprime,...) {
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
  # canopy mean values of Ci and Cc
  .$state$integrated$ci            <- .$leaf$state$ci * fpar/k / .$state$lai
  .$state$integrated$cc            <- .$leaf$state$cc * fpar/k / .$state$lai
  
}



# Two Big Leaf canopy scaling 
###############################
# - accounts for direct and diffuse light separately but only calculates A once for each stream
f_cansys_2bigleaf <- function(.) {
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
###############################
f_cansys_multilayer <- function(.) {
  
  # initialise number of layers
  # explicitely scale canopy N
  # explicitely scale PARsun and PARshade
  # explicitely calculate Asun and Ashade in each layer
  # scale
  
  # initialise layers
  layers      <- ceiling(.$state$lai) # this could be a function specifying either no. of layers or the below lai assignment   
  #?.$init_vert(l=layers)
  
  # canopy layer 
  # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
  # populate layer dataframe      
  .$state$vert$leaf.ca_conc    <- get(.$fnames$can_scale_Ca)(.,1:layers)
  .$state$vert$leaf.vpd        <- get(.$fnames$can_scale_vpd)(.,1:layers)
  .$state$vert$leaf.par        <- get(.$fnames$can_scale_light)(.,1:layers)
  .$state$vert$leaf.leafN_area <- get(.$fnames$can_scale_N)(.,l=1:layers,layers=layers)
  
  # run leaf
  lout   <- as.data.frame(do.call(rbind,lapply(1:layers,.$run_leaf)))
  
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
  # - need to generalise this to allow direct light only, both direct and diffuse (and sunlit and shade leaves)
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



### END ###
