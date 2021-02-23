################################
#
# Canopy system functions
#
# AWalker February 2018
#
################################

source('../../functions/solver_tridiag.R') 



# Big Leaf canopy scaling
###############################

# Sellers et al 1992
f_sys_bigleaf_s1992 <- function(., ... ) {
  # this function was described in Sellers to deal with time intergrated values of fpar and k
  # could write wrapper function or if to pass different k coefficients

  # calculate fPAR, after Sellers 1992, see also Friend 2001
  k    <- .super$state_pars$k_dirprime
  fpar <- 1 - exp(-k*.super$env$lai)

  # set leaf APAR 
  # absorbed light in leaf layer 0 - F_0 * first half of B_2 in Eq 37b (Sellers 1992)
  # APW: same as beerslaw RT function when l = 0
  .super$leaf$env$par[] <- (1-.super$state_pars$lscattering) * .super$state_pars$k_dir * .super$env$par

  # set leaf N0 or Vcmax0 - as with multi-layer model need to choose one and initialise leaf fnames correctly
  .super$leaf$pars$atref$vcmax[] <- .super$pars$vcmax$layer0 
  # APW: is this really necessary? Sets values for leaf, but using layer1, not layer0
  if(.super$fnames$scale$vcmax=='f_scale_two_layer') {
    .super$leaf$pars$atref$vcmax[] <- .$scale.vcmax(1, var='vcmax' )
    .super$leaf$pars$atref$jmax[]  <- .$scale.jmax(1, var='jmax' )
    .super$leaf$pars$f[]           <- .$scale.f(1, var='f' )
    .super$leaf$pars$g1_medlyn[]   <- .$scale.g1(1, var='g1' )
  }  

  # calculate A0
  # leaf model will calculate Vcmax0 and Jmax0 according to leaf process specifications, e.g. from N0, temp, etc
  .super$leaf$run()

  # print leaf state
  if(.$cpars$verbose) {
    print('', quote=F )
    print('Bigleaf leaf pars:', quote=F )
    print(.super$leaf$fnames$etrans)
    print(.super$leaf$pars$theta_j)
    print(.super$leaf$pars$g1_medlyn)
    print(.super$leaf$pars$f)
    print(.super$leaf$pars$atref)
    print(.super$leaf$state_pars)
  }

  # scale
  .super$state$integrated$apar <- fpar * .super$env$par
  # Eq 34 Sellers (1992)
  .super$state$integrated$A              <- .super$leaf$state$A * fpar/k
  .super$state$integrated$rd             <- .super$leaf$state$rd * fpar/k
  .super$state$integrated$Acg_lim        <- if(.super$leaf$state$lim==2) .super$state$integrated$A else 0
  .super$state$integrated$Ajg_lim        <- if(.super$leaf$state$lim==3) .super$state$integrated$A else 0
  .super$state$integrated$Apg_lim        <- if(.super$leaf$state$lim==7) .super$state$integrated$A else 0
  .super$state$integrated$layers_Acg_lim <- if(.super$leaf$state$lim==2) .super$env$lai
  .super$state$integrated$layers_Ajg_lim <- if(.super$leaf$state$lim==3) .super$env$lai
  .super$state$integrated$layers_Apg_lim <- if(.super$leaf$state$lim==7) .super$env$lai
  # convert reistance to conductance, minus minimum conductance, scale, add min conductance multiplied by LAI
  # A combination of Eq 35 and Eq 37f in Sellers (1992)
  .super$state$integrated$gs             <- (1/.super$leaf$state_pars$rs - .super$leaf$pars$g0) * fpar/k + .super$leaf$pars$g0*.super$env$lai 
  # canopy mean values of Ci and Cc - not sure this is correct
  .super$state$integrated$ci             <- .super$leaf$state$ci * fpar/k / .super$env$lai
  .super$state$integrated$cc             <- .super$leaf$state$cc * fpar/k / .super$env$lai
  .super$state$integrated$cb             <- .super$leaf$state$cb * fpar/k / .super$env$lai
}



# Two Big Leaf canopy scaling
###############################
# - accounts for direct and diffuse light & sun and shaded leaves but only calculates A once for each leaf class 
f_sys_2bigleaf <- function(.) {

  # initialise layers
  # sun/shade
  .super$state$vert$sun$apar       <- numeric(.super$pars$layers_2bigleaf)
  .super$state$vert$shade$apar     <- numeric(.super$pars$layers_2bigleaf)
  .super$state$vert$sun$fraction   <- numeric(.super$pars$layers_2bigleaf)  
  .super$state$vert$shade$fraction <- numeric(.super$pars$layers_2bigleaf)
  # leaf traits
  .super$state$vert$leaf$leaf.atref.vcmax <- numeric(.super$pars$layers_2bigleaf)

  # calculate Vcmax sun/shade, other traits yet to add
  .super$state$vert$leaf$leaf.atref.vcmax[] <- .$scale.vcmax(ca_calc_points, var='vcmax' )
  
  # calculate LAI sun and LAI shade - Dai 2004
  lai_sun   <- (1 - exp(-.super$state_pars$k_dir*.super$env$lai)) / .super$state_pars$k_dir
  lai_shade <- .super$env$lai - lai_sun
#  print('lai:')
#  print(lai_sun)
#  print(lai_shade)

  # calculate APARsun and APARshade
  linc           <- .super$env$lai / .super$pars$layers_2bigleaf
  ca_calc_points <- seq((linc-linc*.super$pars$k_layer), (.super$env$lai-linc*.super$pars$k_layer), linc )
  .$rt(ca_calc_points)
#  print(linc)
#  print(ca_calc_points)
#  print('')
#  print('APAR:')
#  print(.super$state$vert$sun$apar) 
#  print(.super$state$vert$shade$apar) 
#  print(.super$state$vert$sun$fraction) 
#  print(.super$state$vert$shade$fraction) 
  
  # scale APAR and traits
  .super$state$vert$sun$apar   <- sum(.super$state$vert$sun$apar*.super$state$vert$sun$fraction*linc) / lai_sun 
  .super$state$vert$shade$apar <- sum(.super$state$vert$shade$apar*.super$state$vert$shade$fraction*linc) / lai_shade
#  print(.super$state$vert$sun$apar) 
#  print(.super$state$vert$shade$apar) 
  vcmax_sun                    <- sum(.super$state$vert$leaf$leaf.atref.vcmax*.super$state$vert$sun$fraction*linc) / lai_sun
  vcmax_shade                  <- sum(.super$state$vert$leaf$leaf.atref.vcmax*.super$state$vert$shade$fraction*linc) / lai_shade

  # calculate Asun and Ashade
  # sun 
  .super$leaf$env$par          <- .super$state$vert$sun$apar
  .super$leaf$pars$atref$vcmax <- vcmax_sun
  leaf_out                     <- .super$leaf$run()
  # assign data to canopy object data structure
  for(vname in names(leaf_out)) .super$state$vert$sun[[vname]][] <- leaf_out[vname]
#  print('')
#  print('Vcmax:')
#  print(vcmax_sun)
#  print(vcmax_shade)
#  print('')
#  print('A etc:')
#  print(leaf_out)
#  print(.super$state$vert$sun) 
  
  # shade  
  .super$leaf$env$par          <- .super$state$vert$shade$apar
  .super$leaf$pars$atref$vcmax <- vcmax_shade
  leaf_out                     <- .super$leaf$run()
  # assign data to canopy object data structure
  for(vname in names(leaf_out)) .super$state$vert$shade[[vname]][] <- leaf_out[vname]
#  print(leaf_out)
#  print(.super$state$vert$shade) 

  # scale & combine
  #for(vname in c('apar','A','gb','gs','gi','g','rd','cb','ci','cc') ) {
  for(vname in c('apar','A','gb','gs','gi','g','rd','ci','cc') ) {
    #print(vname)
    .super$state$integrated[[vname]][] <-
      .super$state$vert$sun[[vname]]*lai_sun + .super$state$vert$shade[[vname]]*lai_shade 
  }
  .super$state$integrated$Acg_lim[] <-
    .super$state$vert$sun$A * (.super$state$vert$sun$lim==2) * lai_sun +
    .super$state$vert$shade$A * (.super$state$vert$shade$lim==2) * lai_shade 
  .super$state$integrated$Ajg_lim[] <-
    .super$state$vert$sun$A * (.super$state$vert$sun$lim==3) * lai_sun +
    .super$state$vert$shade$A * (.super$state$vert$shade$lim==3) * lai_shade 
  .super$state$integrated$Apg_lim[] <-
    .super$state$vert$sun$A * (.super$state$vert$sun$lim==7) * lai_sun +
    .super$state$vert$shade$A * (.super$state$vert$shade$lim==7) * lai_shade 
}



# Multilayer canopy scaling
###############################
f_sys_multilayer <- function(.) {

  # initialise layers
  linc           <- .super$env$lai / .super$pars$layers
  ca_calc_points <- seq((linc-linc*.super$pars$k_layer), (.super$env$lai-linc*.super$pars$k_layer), linc )
  layers         <- .super$pars$layers # this could be a function to dynamically specify the no. of layers
  # APW: consider putting this function in fns to avoid needing to use .super and .=.super
  # APW: could be associated with sys function with a suffix and used for this and 2bigleaf 
  .super$init_vert(.=.super, l=layers ) # reallocating this memory is unnecessary in cases where layers is a fixed parameter.
  #print(ca_calc_points)

  # canopy leaf layer parameters
  # APW: could be a separate function shared by multilayer and 2bigleaf
  .super$state$vert$leaf$leaf.leafN_area[]   <- .$scale.n(ca_calc_points, var='n' )
  .super$state$vert$leaf$leaf.atref.vcmax[]  <- .$scale.vcmax(ca_calc_points, var='vcmax' )
  leaf_vars <- c('leaf.leafN_area', 'leaf.atref.vcmax' )
  if(.super$fnames$scale$vcmax=='f_scale_two_layer') {
    .super$state$vert$leaf$leaf.atref.jmax[] <- .$scale.jmax(ca_calc_points, var='jmax' )
    .super$state$vert$leaf$leaf.f[]          <- .$scale.f(ca_calc_points, var='f' )
    .super$state$vert$leaf$leaf.g1_medlyn[]  <- .$scale.g1(ca_calc_points, var='g1' )
    leaf_vars <- c(leaf_vars, 'leaf.atref.jmax', 'leaf.f', 'leaf.g1_medlyn' )
  }  

  # canopy leaf layer environment
  # APW: when called as a character string for some reason the proto . is not fully recognised within the function
  #      parameters can be found, but aliased function calls fail, maybe named function calls work
  .super$state$vert$leaf$leaf.ca_conc[]      <- .$scale.ca_conc(ca_calc_points, var='ca_conc', vlist='env' )
  .super$state$vert$leaf$leaf.vpd[]          <- .$scale.vpd(ca_calc_points, var='vpd', vlist='env' )
  leaf_vars <- c(leaf_vars, 'leaf.ca_conc', 'leaf.vpd', 'leaf.par' )
  #print(.super$state$vert$leaf)

  # Light scaling
  .$rt(ca_calc_points)

  # sunlit leaves / direct light
  .super$state$vert$leaf$leaf.par[] <- .super$state$vert$sun$apar
  # create leaf environment  matrix
  lmatrix  <- vapply(.super$state$vert$leaf[leaf_vars], function(v) v, numeric(layers) )
  # run leaf
  leaf_out <- vapply(1:layers, .$run_leaf, .super$leaf$output(), df=lmatrix )
  # assign data to canopy object data structure
  for(vname in row.names(leaf_out)) .super$state$vert$sun[[vname]][] <- leaf_out[vname,]
  if(.$cpars$verbose) {
    print('', quote=F )
    print('Sun leaves:', quote=F )
    print(lmatrix)
    print(leaf_out)
    print(.super$state$vert$sun)
  }

  # shade leaves
  if(any(.super$state$vert$sun$fraction < 1) ) {
    .super$state$vert$leaf$leaf.par[] <- .super$state$vert$shade$apar
    # create leaf environment matrix
    lmatrix  <- vapply(.super$state$vert$leaf[leaf_vars], function(v) v, numeric(layers) )
    # run leaf
    leaf_out <- vapply(1:layers, .$run_leaf, .super$leaf$output(), df=lmatrix )
    # assign data to canopy object data structure
    for(vname in row.names(leaf_out)) .super$state$vert$shade[[vname]][] <- leaf_out[vname,]
    if(.$cpars$verbose) {
      print('', quote=F )
      print('Shade leaves:', quote=F )
      print(lmatrix)
      print(leaf_out)
      print(.super$state$vert$shade)
    }
  }

  # combine sun and shade leaves
  for(vname in c('apar','A','gb','gs','gi','g','rd','cb','ci','cc') ) {
    .super$state$vert$layer[[vname]][] <-
      .super$state$vert$sun[[vname]] * .super$state$vert$sun$fraction + .super$state$vert$shade[[vname]] * .super$state$vert$shade$fraction
  }
  # partition A among limiting rates
  .super$state$vert$layer$Acg_lim[] <-
    .super$state$vert$sun$A * (.super$state$vert$sun$lim==2) * .super$state$vert$sun$fraction +
    .super$state$vert$shade$A * (.super$state$vert$shade$lim==2) * .super$state$vert$shade$fraction
  .super$state$vert$layer$Ajg_lim[] <- 
    .super$state$vert$sun$A * (.super$state$vert$sun$lim==3) * .super$state$vert$sun$fraction + 
    .super$state$vert$shade$A * (.super$state$vert$shade$lim==3) * .super$state$vert$shade$fraction
  .super$state$vert$layer$Apg_lim[] <- 
    .super$state$vert$sun$A * (.super$state$vert$sun$lim==7) * .super$state$vert$sun$fraction + 
    .super$state$vert$shade$A * (.super$state$vert$shade$lim==7) * .super$state$vert$shade$fraction
  if(.$cpars$verbose) {
    print('', quote=F )
    print('Layer:', quote=F )
    print(.super$state$vert$layer)
    print(.super$state$vert$sun['lim'])
    print(.super$state$vert$shade['lim'])
  }

  # integrate canopy layers
  # canopy sum values
  for(vname in c('apar','A','gb','gs','gi','g','rd','Acg_lim','Ajg_lim','Apg_lim') ) {
    .super$state$integrated[[vname]][] <- sum(.super$state$vert$layer[[vname]]) * linc
  }
  # canopy mean values - not sure that this is correct
  .super$state$integrated$cc[] <- sum(.super$state$vert$layer$cc) * linc / .super$env$lai
  .super$state$integrated$ci[] <- sum(.super$state$vert$layer$ci) * linc / .super$env$lai
  .super$state$integrated$cb[] <- sum(.super$state$vert$layer$cb) * linc / .super$env$lai
}



### END ###
