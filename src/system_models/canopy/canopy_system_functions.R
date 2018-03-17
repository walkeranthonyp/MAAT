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
  # leaf model will calculate Vcmax0 and Jmax0 according to leaf process specifications, e.g. from N0, temp, etc 
  .$leaf$run()
  
  # scale
  .$state$integrated$A             <- .$leaf$state$A * fpar/k
  .$state$integrated$respiration   <- .$leaf$state$respiration * fpar/k
  .$state$integrated$Acg_lim        <- if(.$leaf$state$lim==2) .$state$integrated$A else 0
  .$state$integrated$Ajg_lim        <- if(.$leaf$state$lim==3) .$state$integrated$A else 0 
  .$state$integrated$Apg_lim        <- if(.$leaf$state$lim==7) .$state$integrated$A else 0
  .$state$integrated$layers_Acg_lim <- if(.$leaf$state$lim==2) .$state$lai
  .$state$integrated$layers_Ajg_lim <- if(.$leaf$state$lim==3) .$state$lai
  .$state$integrated$layers_Apg_lim <- if(.$leaf$state$lim==7) .$state$lai
  # resistance - convert to conductance, minus minimum conductance, scale, add min conductance multiplied by LAI, convert back to resistance
  # - a somewhat complicated version of eq 37f & 35 from Sellers 1992 and a conductance to resistance conversion
  .$state$integrated$rs            <- 1 / ( (1/.$leaf$state$rs - .$leaf$pars$g0) * fpar/k + .$leaf$pars$g0*.$state$lai )
  # canopy mean values of Ci and Cc
  .$state$integrated$ci            <- .$leaf$state$ci * fpar/k / .$state$lai
  .$state$integrated$cc            <- .$leaf$state$cc * fpar/k / .$state$lai
  
}



# Two Big Leaf canopy scaling 
###############################
# - accounts for direct and diffuse light separately but only calculates A once for each radiation type 
f_cansys_2bigleaf <- function(.) {
  # Thornton calculates the mean of all these canopy values, what does Dai do? 
  
  # calculate LAIsun and LAIshade - Dai 2004
  Lsun   <- (1 - exp(-.$state_pars$k_dir*.$state$lai)) / .$state_pars$k_dir
  Lshade <- .$state$lai - Lsun 
  
  # calculate APARsun and APARshade  
  # APARshade is the scattered part of the direct beam plus diffuse radiation
  
  # APARsun is the direct beam plus the scattered part of the direct beam plus diffuse radiation
  
  # calculate Nsun and Nshade - this doesn't really make sense as leaf N is not able to vary on the time scales that on which sun and shade leaves vary  
  
  # calculate Asun and Ashade
  
  # scale
}



# Multilayer canopy scaling
###############################
f_cansys_multilayer <- function(.) {
  
  # initialise layers
  linc           <- .$state$lai / .$pars$layers
  ca_calc_points <- seq(linc, .$state$lai, linc ) 
  layers         <- .$pars$layers
  #layers         <- ceiling(.$state$lai) # this could be a function specifying either no. of layers or the below lai assignment   
  .$init_vert(l=layers) # reallocating this memory is unnecessary in cases where layers is a fixed parameter. 
  
  # canopy leaf layer properties 
  .$state$vert$leaf.leafN_area[] <- get(.$fnames$can_scale_N)(., ca_calc_points )
  .$state$vert$leaf.ca_conc[]    <- get(.$fnames$can_scale_Ca)(., ca_calc_points )
  .$state$vert$leaf.vpd[]        <- get(.$fnames$can_scale_vpd)(., ca_calc_points )

  # Light scaling  
  # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
  get(.$fnames$can_scale_light)(.,ca_calc_points )

  # direct light
  #.$state$vert$leaf.par        <- get(.$fnames$can_scale_light)(.,ca_calc_points )
  .$state$vert$leaf.par[] <- .$state$vert$apar_sun 

  # create leaf environment  matrix
  lmatrix <- vapply(.$state$vert[c('leaf.leafN_area','leaf.ca_conc','leaf.vpd','leaf.par')], function(v) v, numeric(layers) )  
  print(lmatrix)

  # run leaf
  leaf_out <- vapply(1:layers, .$run_leaf, .$leaf$output(), df=lmatrix )
  print(leaf_out)
 
  # assign data
  .$state$vert$A[]           <- leaf_out['A',]
  .$state$vert$cc[]          <- leaf_out['cc',]
  .$state$vert$ci[]          <- leaf_out['ci',]
  .$state$vert$lim[]         <- leaf_out['lim',]
  .$state$vert$respiration[] <- leaf_out['respiration',]      
  .$state$vert$rb[]          <- leaf_out['rb',]
  .$state$vert$rs[]          <- leaf_out['rs',]
  .$state$vert$ri[]          <- leaf_out['ri',]
  
  #integrate canopy layers
  # - need to generalise this to allow direct light only, both direct and diffuse (and sunlit and shade leaves)
  # canopy sum values
  .$state$integrated$A[]              <- sum(.$state$vert$A) * linc
  .$state$integrated$respiration[]    <- sum(.$state$vert$respiration) * linc
  .$state$integrated$Acg_lim[]        <- sum(.$state$vert$A * (.$state$lim==2)) * linc 
  .$state$integrated$Ajg_lim[]        <- sum(.$state$vert$A * (.$state$lim==3)) * linc
  .$state$integrated$Apg_lim[]        <- sum(.$state$vert$A * (.$state$lim==7)) * linc
  .$state$integrated$layers_Acg_lim[] <- sum(.$state$vert$lim==2)
  .$state$integrated$layers_Ajg_lim[] <- sum(.$state$vert$lim==3)
  .$state$integrated$layers_Apg_lim[] <- sum(.$state$vert$lim==7)
  .$state$integrated$rb[]             <- 1 / sum(1/.$state$vert$rb * linc )
  .$state$integrated$rs[]             <- 1 / sum(1/.$state$vert$rs * linc )
  .$state$integrated$ri[]             <- 1 / sum(1/.$state$vert$ri * linc )
  # canopy mean values
  .$state$integrated$cc[]             <- sum(.$state$vert$cc) / .$state$lai * linc
  .$state$integrated$ci[]             <- sum(.$state$vert$ci) / .$state$lai * linc
  
}



### END ###
