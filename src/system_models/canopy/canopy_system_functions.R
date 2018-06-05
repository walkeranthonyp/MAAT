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
  .$state$integrated$A              <- .$leaf$state$A * fpar/k
  .$state$integrated$respiration    <- .$leaf$state$respiration * fpar/k
  .$state$integrated$Acg_lim        <- if(.$leaf$state$lim==2) .$state$integrated$A else 0
  .$state$integrated$Ajg_lim        <- if(.$leaf$state$lim==3) .$state$integrated$A else 0 
  .$state$integrated$Apg_lim        <- if(.$leaf$state$lim==7) .$state$integrated$A else 0
  .$state$integrated$layers_Acg_lim <- if(.$leaf$state$lim==2) .$state$lai
  .$state$integrated$layers_Ajg_lim <- if(.$leaf$state$lim==3) .$state$lai
  .$state$integrated$layers_Apg_lim <- if(.$leaf$state$lim==7) .$state$lai
  # resistance - convert to conductance, minus minimum conductance, scale, add min conductance multiplied by LAI, convert back to resistance
  # - a somewhat complicated version of eq 37f & 35 from Sellers 1992 and a conductance to resistance conversion
  .$state$integrated$rs             <- 1 / ( (1/.$leaf$state$rs - .$leaf$pars$g0) * fpar/k + .$leaf$pars$g0*.$state$lai )
  # canopy mean values of Ci and Cc
  .$state$integrated$ci             <- .$leaf$state$ci * fpar/k / .$state$lai
  .$state$integrated$cc             <- .$leaf$state$cc * fpar/k / .$state$lai
  
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
  get(.$fnames$rt)(., 1:.$state$lai )
  #take the mean of these? - need to check, is there an analytical method?
  
  # Leaf environment 
  .$leaf$env$ca_conc <- .$env$ca_conc
  .$leaf$env$vpd     <- .$env$vpd
  
  # calculate Nsun and Nshade - this doesn't really make sense as leaf N is not able to vary on the time scales that on which sun and shade leaves vary  
  .$leaf$state$leafN_area <- .$state$totalN * k / fpar
  
  # calculate Asun and Ashade
  
  # scale & combine
}



# Multilayer canopy scaling
###############################
f_cansys_multilayer <- function(.) {
 
  # initialise layers
  linc           <- .$state$lai / .$pars$layers
  ca_calc_points <- seq(linc, .$state$lai, linc ) 
  layers         <- .$pars$layers # this could be a function to dynamically specify the no. of layers 
  .$init_vert(l=layers) # reallocating this memory is unnecessary in cases where layers is a fixed parameter. 
  
  # canopy leaf layer properties 
  .$state$vert$leaf$leaf.leafN_area[] <- get(.$fnames$scale_n)(.,   ca_calc_points )
  .$state$vert$leaf$leaf.ca_conc[]    <- get(.$fnames$scale_ca)(.,  ca_calc_points )
  .$state$vert$leaf$leaf.vpd[]        <- get(.$fnames$scale_vpd)(., ca_calc_points )

  # Light scaling  
  get(.$fnames$rt)(., ca_calc_points )

  # sunlit leaves / direct light
  .$state$vert$leaf$leaf.par[] <- .$state$vert$sun$apar 
  # create leaf environment  matrix
  lmatrix  <- vapply(.$state$vert$leaf[c('leaf.leafN_area','leaf.ca_conc','leaf.vpd','leaf.par')], function(v) v, numeric(layers) )
  # run leaf
  leaf_out <- vapply(1:layers, .$run_leaf, .$leaf$output(), df=lmatrix )
  # assign data to canopy object data structure
  for(vname in row.names(leaf_out)) .$state$vert$sun[[vname]][] <- leaf_out[vname,]
  if(.$cpars$verbose) {
    print('Sun leaves:', quote=F )
    print(lmatrix)
    print(leaf_out)
  }

  # shade leaves
  if(any(.$state$vert$sun$fraction < 1) ) { 
    .$state$vert$leaf$leaf.par[] <- .$state$vert$shade$apar 
    # create leaf environment matrix
    lmatrix  <- vapply(.$state$vert$leaf[c('leaf.leafN_area','leaf.ca_conc','leaf.vpd','leaf.par')], function(v) v, numeric(layers) )
    # run leaf
    leaf_out <- vapply(1:layers, .$run_leaf, .$leaf$output(), df=lmatrix )
    # assign data to canopy object data structure
    for(vname in row.names(leaf_out)) .$state$vert$shade[[vname]][] <- leaf_out[vname,]
    if(.$cpars$verbose) {
      print('Shade leaves:', quote=F )
      print(lmatrix)
      print(leaf_out)
    }
  }

  # combine sun and shade leaves  
  for(vname in names(.$state$vert$layer)) .$state$vert$layer[[vname]][] <- .$state$vert$sun[[vname]] * .$state$vert$sun$fraction + .$state$vert$shade[[vname]] * .$state$vert$shade$fraction

  # integrate canopy layers
  # canopy sum values
  .$state$integrated$A[]              <- sum(.$state$vert$layer$A) * linc
  .$state$integrated$respiration[]    <- sum(.$state$vert$layer$respiration) * linc
  .$state$integrated$Acg_lim[]        <- sum(.$state$vert$layer$A * (.$state$lim==2)) * linc 
  .$state$integrated$Ajg_lim[]        <- sum(.$state$vert$layer$A * (.$state$lim==3)) * linc
  .$state$integrated$Apg_lim[]        <- sum(.$state$vert$layer$A * (.$state$lim==7)) * linc
  .$state$integrated$layers_Acg_lim[] <- sum(.$state$vert$layer$lim==2)
  .$state$integrated$layers_Ajg_lim[] <- sum(.$state$vert$layer$lim==3)
  .$state$integrated$layers_Apg_lim[] <- sum(.$state$vert$layer$lim==7)
  .$state$integrated$rb[]             <- 1 / sum(1/.$state$vert$layer$rb * linc )
  .$state$integrated$rs[]             <- 1 / sum(1/.$state$vert$layer$rs * linc )
  .$state$integrated$ri[]             <- 1 / sum(1/.$state$vert$layer$ri * linc )
  # canopy mean values
  .$state$integrated$cc[]             <- sum(.$state$vert$layer$cc) / .$state$lai * linc
  .$state$integrated$ci[]             <- sum(.$state$vert$ci) / .$state$lai * linc
  
}



### END ###
