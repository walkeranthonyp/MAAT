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
f_sys_bigleaf_s1992 <- function(., k=.super$state_pars$k_dirprime, ... ) {
  # this function was described in Sellers to deal with time intergrated values of fpar and k
  # could write wrapper function or if to pass different k coefficients

  # calculate fPAR, after Sellers 1992, see also Friend 2001
  fpar <- 1 - exp(-k*.super$state$lai)

  # set leaf environment
  # incident light - F_0 * first half of B_2 in Eq 37b (Sellers 1992)
  .super$leaf$env$par     <- .super$state_pars$k_dir * (1-.super$state_pars$lscattering) * .super$env$par
  # assume no variation in CO2 concentration, VPD, and T
  #.super$leaf$env$ca_conc <- .super$env$ca_conc
  #.super$leaf$env$vpd     <- .super$env$vpd
  #.super$leaf$env$temp    <- .super$env$temp

  # set leaf N0 or Vcmax0 - as with multi-layer model need to choose one and initialise leaf fnames correctly
  .super$leaf$state$leafN_area <- .super$state$totalN * k / fpar
  .super$leaf$pars$atref$vcmax <- .super$pars$vcmax0 

  # calculate A0
  # leaf model will calculate Vcmax0 and Jmax0 according to leaf process specifications, e.g. from N0, temp, etc
  .super$leaf$run()

  # scale
  .super$state$integrated$apar <- fpar * .super$env$par
  # Eq 34 Sellers (1992)
  .super$state$integrated$A              <- .super$leaf$state$A * fpar/k
  .super$state$integrated$rd             <- .super$leaf$state$rd * fpar/k
  .super$state$integrated$Acg_lim        <- if(.super$leaf$state$lim==2) .super$state$integrated$A else 0
  .super$state$integrated$Ajg_lim        <- if(.super$leaf$state$lim==3) .super$state$integrated$A else 0
  .super$state$integrated$Apg_lim        <- if(.super$leaf$state$lim==7) .super$state$integrated$A else 0
  .super$state$integrated$layers_Acg_lim <- if(.super$leaf$state$lim==2) .super$state$lai
  .super$state$integrated$layers_Ajg_lim <- if(.super$leaf$state$lim==3) .super$state$lai
  .super$state$integrated$layers_Apg_lim <- if(.super$leaf$state$lim==7) .super$state$lai
  # convert reistance to conductance, minus minimum conductance, scale, add min conductance multiplied by LAI
  # A combination of Eq 35 and Eq 37f in Sellers (1992)
  .super$state$integrated$gs             <- (1/.super$leaf$state$rs - .super$leaf$pars$g0) * fpar/k + .super$leaf$pars$g0*.super$state$lai
  # canopy mean values of Ci and Cc - not sure this is correct
  .super$state$integrated$ci             <- .super$leaf$state$ci * fpar/k / .super$state$lai
  .super$state$integrated$cc             <- .super$leaf$state$cc * fpar/k / .super$state$lai
  .super$state$integrated$cb             <- .super$leaf$state$cb * fpar/k / .super$state$lai

}



# Two Big Leaf canopy scaling
###############################
# - accounts for direct and diffuse light separately but only calculates A once for each radiation type
f_sys_2bigleaf <- function(.) {
  # Thornton calculates the mean of all these canopy values, what does Dai do?

  # calculate LAIsun and LAIshade - Dai 2004
  Lsun   <- (1 - exp(-.super$state_pars$k_dir*.super$state$lai)) / .super$state_pars$k_dir
  Lshade <- .super$state$lai - Lsun

  # calculate APARsun and APARshade
  # APARshade is the scattered part of the direct beam plus diffuse radiation
  .$rt(1:.super$state$lai)
  #take the mean of these? - need to check, is there an analytical method?

  # Leaf environment
  #.super$leaf$env$ca_conc <- .super$env$ca_conc
  #.super$leaf$env$vpd     <- .super$env$vpd

  # calculate Nsun and Nshade - this doesn't really make sense as leaf N is not able to vary on the time scales that on which sun and shade leaves vary
  .super$leaf$state$leafN_area <- .super$state$totalN * k / fpar

  # calculate Asun and Ashade

  # scale & combine
}



# Multilayer canopy scaling
###############################
f_sys_multilayer <- function(.) {

  # initialise layers
  # k_layer determines where in the layer photosynthesis etc is calculated, a value of 0.5 calculates at the center of the layer
  linc           <- .super$state$lai / .super$pars$layers
  ca_calc_points <- seq((linc-linc*.super$pars$k_layer), (.super$state$lai-linc*.super$pars$k_layer), linc )
  layers         <- .super$pars$layers # this could be a function to dynamically specify the no. of layers
  # consider putting this function in fns to avoid needing to use .super and .=.super
  .super$init_vert(.=.super, l=layers ) # reallocating this memory is unnecessary in cases where layers is a fixed parameter.
  #print(ca_calc_points)

  # canopy leaf layer properties
  .super$state$vert$leaf$leaf.leafN_area[]  <- .$scale_n(ca_calc_points)
  .super$state$vert$leaf$leaf.atref.vcmax[] <- .$scale_vcmax(ca_calc_points)
  .super$state$vert$leaf$leaf.ca_conc[]     <- .$scale_ca(ca_calc_points)
  .super$state$vert$leaf$leaf.vpd[]         <- .$scale_vpd(ca_calc_points)

  # Light scaling
  .$rt(ca_calc_points)

  # sunlit leaves / direct light
  .super$state$vert$leaf$leaf.par[] <- .super$state$vert$sun$apar
  # create leaf environment  matrix
  lmatrix  <- vapply(.super$state$vert$leaf[c('leaf.leafN_area','leaf.atref.vcmax','leaf.ca_conc','leaf.vpd','leaf.par')], function(v) v, numeric(layers) )
  # run leaf
  leaf_out <- vapply(1:layers, .$run_leaf, .super$leaf$output(), df=lmatrix )
  # assign data to canopy object data structure
  for(vname in row.names(leaf_out)) .super$state$vert$sun[[vname]][] <- leaf_out[vname,]
  if(.$cpars$verbose) {
    print('Sun leaves:', quote=F )
    print(lmatrix)
    print(leaf_out)
  }

  # shade leaves
  if(any(.super$state$vert$sun$fraction < 1) ) {
    .super$state$vert$leaf$leaf.par[] <- .super$state$vert$shade$apar
    # create leaf environment matrix
    lmatrix  <- vapply(.super$state$vert$leaf[c('leaf.leafN_area','leaf.atref.vcmax','leaf.ca_conc','leaf.vpd','leaf.par')], function(v) v, numeric(layers) )
    # run leaf
    leaf_out <- vapply(1:layers, .$run_leaf, .super$leaf$output(), df=lmatrix )
    # assign data to canopy object data structure
    for(vname in row.names(leaf_out)) .super$state$vert$shade[[vname]][] <- leaf_out[vname,]
    if(.$cpars$verbose) {
      print('Shade leaves:', quote=F )
      print(lmatrix)
      print(leaf_out)
    }
  }

  # combine sun and shade leaves
  for(vname in c('apar','A','gb','gs','gi','g','rd','cb','ci','cc') ) {
    .super$state$vert$layer[[vname]][] <-
      .super$state$vert$sun[[vname]] * .super$state$vert$sun$fraction + .super$state$vert$shade[[vname]] * .super$state$vert$shade$fraction
  }

  # partition A among limiting rates
  .super$state$vert$layer$Acg_lim[] <- .super$state$vert$sun$A * (.super$state$vert$sun$lim==2) * .super$state$vert$sun$fraction + .super$state$vert$shade$A * (.super$state$vert$shade$lim==2) * .super$state$vert$shade$fraction
  .super$state$vert$layer$Ajg_lim[] <- .super$state$vert$sun$A * (.super$state$vert$sun$lim==3) * .super$state$vert$sun$fraction + .super$state$vert$shade$A * (.super$state$vert$shade$lim==3) * .super$state$vert$shade$fraction
  .super$state$vert$layer$Acg_lim[] <- .super$state$vert$sun$A * (.super$state$vert$sun$lim==7) * .super$state$vert$sun$fraction + .super$state$vert$shade$A * (.super$state$vert$shade$lim==7) * .super$state$vert$shade$fraction

  # integrate canopy layers
  # canopy sum values
  for(vname in c('apar','A','gb','gs','gi','g','rd','Acg_lim','Ajg_lim','Apg_lim') ) {
    .super$state$integrated[[vname]][] <- sum(.super$state$vert$layer[[vname]]) * linc
  }
  # canopy mean values - not sure that this is correct
  .super$state$integrated$cc[] <- sum(.super$state$vert$layer$cc) * linc / .super$state$lai
  .super$state$integrated$ci[] <- sum(.super$state$vert$layer$ci) * linc / .super$state$lai
  .super$state$integrated$cb[] <- sum(.super$state$vert$layer$cb) * linc / .super$state$lai

}



### END ###
