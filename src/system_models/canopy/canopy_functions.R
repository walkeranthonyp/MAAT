################################
#
# Canopy level functions
# 
# AWalker December 2015
#
################################



# Canopy discretisation 
################################

# break canopy into <layers> equally sized bins
f_canopy_discretisation_fixednumber_const <- function(.) {

  .super$state_pars$linc <- linc <- .super$env$lai/.super$pars$layers
  .super$state_pars$ca_calc_points <- seq((0+linc*.super$pars$k_layer), (.super$env$lai-linc+linc*.super$pars$k_layer), linc )
  .super$state_pars$linc <- rep(linc, length(.super$state_pars$ca_calc_points) ) 
}

# break canopy into equally sized bins, up to and including <maxlai>
f_canopy_discretisation_fixedbins_const <- function(.) {

  .super$state_pars$linc <- linc <- .super$pars$linc
  cum_lai <- seq(1, .super$env$lai, linc )
  .super$state_pars$ca_calc_points <- c(0,cum_lai) + linc*.super$pars$k_layer
  .super$state_pars$cum_lai <- c(cum_lai, cum_lai[length(cum_lai)] + linc )

  # calculate linc etc for final layer
  lai_diff <- .super$env$lai - .super$state_pars$cum_lai[length(subset)-1]  
  .super$state_pars$cum_lai[length(subset)]        <- .super$env$lai
  .super$state_pars$linc[length(subset)]           <- lai_diff 
  #.super$state_pars$ca_calc_points[length(subset)] <- .super$state_pars$ca_calc_points[subset]
}

# break canopy into <layers> exponentially sized bins
f_canopy_discretisation_fixednumber_exp <- function(.) {
  # fixed number of bins, exponential discretisation
  # - fixed exponent, variable initial layer width

  layer1 <-  .super$env$lai / sum(exp(.super$pars$k$can_disc * 1:.super$pars$layers))  
  .super$state_pars$linc <- layer1 * exp(.super$pars$k$can_disc * 1:.super$pars$layers)   
  .super$state_pars$cum_lai <- cum_lai <- cumsum(.super$state_pars$linc)
   if(abs(cum_lai[length(cum_lai)]-.super$env$lai)>1e-5)  
     stop(paste('Canopy layer sum (',cum_lai[length(cum_lai)],') != lai (',.super$env$lai,')')) 
  .super$state_pars$ca_calc_points <- c(0,cum_lai[-length(cum_lai)]) + .super$state_pars$linc*.super$pars$k_layer

}

# break canopy into <layers> exponentially sized bins variable exponent with LAI
f_canopy_discretisation_fixednumber_expvarlai <- function(.) {
  # fixed number of bins, exponential discretisation
  # - exponent varies as a power law of LAI, variable initial layer width

  k_can_disc <- .super$env$lai^0.095 - 1  
  layer1     <- .super$env$lai / sum(exp(k_can_disc * 1:.super$pars$layers))  
  .super$state_pars$linc <- layer1 * exp(k_can_disc * 1:.super$pars$layers)   
  .super$state_pars$cum_lai <- cum_lai <- cumsum(.super$state_pars$linc)
   if(abs(cum_lai[length(cum_lai)]-.super$env$lai)>1e-5)  
     stop(paste('Canopy layer sum (',cum_lai[length(cum_lai)],') != lai (',.super$env$lai,')')) 
  .super$state_pars$ca_calc_points <- c(0,cum_lai[-length(cum_lai)]) + .super$state_pars$linc*.super$pars$k_layer

}

# break canopy into exponentially increasing sized bins, up to and including <maxlai>
f_canopy_discretisation_fixedbins_exp <- function(.) {

  .super$state_pars$linc <- .super$pars$can_disc$layer1 * exp(.super$pars$k$can_disc * 0:1000 )  
  cum_lai <- cumsum(.super$state_pars$linc)
   if(cum_lai[length(cum_lai)]<.super$env$lai) stop('Canopy discretisation parameters result in canopy layers > max layers of 1000.') 
  .super$state_pars$ca_calc_points <- c(0,cum_lai[-length(cum_lai)]) + .super$state_pars$linc*.super$pars$k_layer
  
  # reduce vectors to include maxlai and no more 
  subset  <- which(cum_lai<.super$env$lai)
  subset  <- c(subset,length(subset)+1)
  .super$state_pars$cum_lai        <- cum_lai[subset]
  .super$state_pars$linc           <- .super$state_pars$linc[subset]
  .super$state_pars$ca_calc_points <- .super$state_pars$ca_calc_points[subset]

  # calculate linc etc for final layer
  lai_diff <- .super$env$lai - .super$state_pars$cum_lai[length(subset)-1]  
  .super$state_pars$cum_lai[length(subset)]        <- .super$env$lai
  .super$state_pars$linc[length(subset)]           <- lai_diff 
  #.super$state_pars$ca_calc_points[length(subset)] <- .super$state_pars$ca_calc_points[subset]
}



# RT schemes 
################################

# Beer's Law
# - these need to be a number of different but similar functions that allow all the commonly used implementations of Beer's Law
# - e.g. accounting for diffuse light, adjusting k by m or not, etc.

# Beer's law - from Sellers (1992) and citing Goudriaan (1977) 
f_rt_beerslaw_goudriaan <- function(., l=.super$state_pars$ca_calc_points ) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .super$state$vert$sun$apar[]     <- (1-.super$state_pars$lscattering)*.super$state_pars$k_dir*.super$env$par * exp(-.super$state_pars$k_dirprime*l)
  .super$state$vert$sun$fraction[] <- 1.0 
}


# misconstrued Beer's Law
# - ignores conversion of incident light per unit ground area to incident light per unit leaf area
# - considers scattering in absorption, leaf object does not when nested in canopy object
# - but scattering is not considered in RT, k_dir not k_dirprime is used
f_rt_beerslaw <- function(., l=.super$state_pars$ca_calc_points ) {
  # calculates direct beam light attenuation through the canopy
  # for use with a multilayer canopy
  # returns incident radiation in canopy layer 'l', can take 'l' as a vector
  
  .super$state$vert$sun$apar[]     <- (1-.super$state_pars$lscattering)*.super$env$par * exp(-.super$state_pars$k_dir*l)
  .super$state$vert$sun$fraction[] <- 1.0 
}


# Goudriaan
f_rt_goudriaan <- function(., l=.super$state_pars$ca_calc_points ) {
  # as described in Walker et al. (2017) New Phyt, supplement; from Spitters 1986 and Wang 2003 Func.Plant Biol.     
  # all below values are for visible wavelengths 

  # calculate absorbed direct & diffuse radiation  
  qshade <- (1-.super$state_pars$alb_diff)    * .super$env$par_diff * .super$state_pars$k_diffprime * exp(-.super$state_pars$k_diffprime*l) + 
            (1-.super$state_pars$alb_dir)     * .super$env$par_dir  * .super$state_pars$k_dirprime  * exp(-.super$state_pars$k_dirprime*l) - 
            (1-.super$state_pars$lscattering) * .super$env$par_dir  * .super$state_pars$k_dir       * exp(-.super$state_pars$k_dir*l)
  qshade <- ifelse(qshade<0, 0, qshade)
  
  .super$state$vert$shade$apar[] <- qshade
  .super$state$vert$sun$apar[]   <- (1-.super$state_pars$lscattering)*.super$state_pars$k_dir*.super$env$par_dir + .super$state$vert$shade$apar
  
  # calculate fraction sunlit vs shaded leaves
  .super$state$vert$sun$fraction[]   <- .super$pars$can_clump * exp(-.super$state_pars$k_dir*l) 
  #.super$state$vert$sun$fraction[]   <- exp(-.super$state_pars$k_dir*l)
  .super$state$vert$shade$fraction[] <- 1 - .super$state$vert$sun$fraction
}


# calculate transmission of diffuse radiation through leaf layer dLAI thick 
# - by calculating transmissio from all sky angles
# - zi, dLAI, can_clump, and chi_l (used in gz) can all be vectors
# - but only either zi or dLAI, can_clump, chi_l 
f_td <- function(., zi, dLAI ) { 
  gz  <- .$gz(z=zi) 
  exp(-gz/cos(zi)*dLAI*.super$pars$can_clump) * sin(zi)*cos(zi)
}


# Norman
# - using Bonan's tridiagonal solution
f_rt_norman <- function(., l=.super$state_pars$ca_calc_points, 
                        leaf_reflectance=.super$pars$leaf_reflectance,
                        leaf_transmitance=.super$pars$leaf_transmitance,
                        soil_reflectance_dir=.super$pars$soil_reflectance_dir,
                        soil_reflectance_diff=.super$pars$soil_reflectance_diff
                        ) {
  
  # Norman scheme calculates layers + 1 to include soil 
  # code by Bonan has the first layer as ground and last being the top canopy layer
  # - so l reversed
  nlayers <- length(l)
  sumlai  <- c(NA, rev(l) )
  dLAI    <- if(length(.super$state_pars$linc)==1) {
               c(NA, rep(.super$state_pars$linc, nlayers )) 
             } else {
               c(NA, rev(.super$state_pars$linc) ) 
             }

  # APW: these initial calculations assume that leaf angle distribution is constant through the canopy    
  # APW: k_dir assumed fixed here   
  .super$state$vert$sun$fraction[]   <- .super$pars$can_clump * exp(-.super$state_pars$k_dir*rev(l)) 
  .super$state$vert$shade$fraction[] <- 1 - .super$state$vert$sun$fraction
 
  # tb (frac direct beam not intercepted in each layer)
  # APW: k_dir assumed fixed here  
  tb <- exp(-.super$state_pars$k_dir * dLAI )
 
  # td (frac diffuse beam not intercepted in each layer)
  # APW: k_dir assumed fixed here in gz() call  
  # APW: if variable leaf angle with canopy layer needed this will want to loop over dLAI and vectorise zi  
  # APW: - but maybe it can work, gz should return a vector if chi_l is a vector 
#  td <- 0
#  for(zi in .super$state_pars$zi) {
#    gz  <- .$gz(z=zi) 
#    td  <- td + exp(-gz/cos(zi)*dLAI*.super$pars$can_clump) * sin(zi)*cos(zi)
#  }                     
#  td <- td*2*.super$state_pars$delta_zi
  td_mat <- sapply(.super$state_pars$zi, .$td, dLAI=dLAI )
  td     <- 2*apply(td_mat, 1, sum ) * .super$state_pars$delta_zi

  # tbcum (frac toc direct beam incident at top of each layer)
  tbcum     <- rep(NA,nlayers+1)
  iv        <- nlayers+1
  tbcum[iv] <- 1
  cumlai    <- 0
  for(iv in (nlayers+1):2) {
    cumlai      <- cumlai + dLAI[iv]
    # APW: k_dir assumed fixed here, can this be developed as f(tb) - should be cummulative product of tb   
    tbcum[iv-1] <- exp(-.super$state_pars$k_dir * cumlai )
  }
  
  # initiate arrays
  swup <- swdn <- a <- b <- c <- d <- rep(0, nlayers+1 )
  direct <- diffuse <- sun <- shade <- rep(NA, nlayers+1 )
  
  # Soil: upward flux
  m  <- 1
  iv <- 1
  a[m] <- 0
  b[m] <- 1
  c[m] <- -soil_reflectance_diff
  d[m] <- .super$env$par_dir*tbcum[m] * soil_reflectance_dir
  
  # Soil: downward flux
  # APW: leaf_refelctance and leaf_transmitance assumed constant here and below 
  refld <- (1-td[iv+1]) * leaf_reflectance 
  trand <- (1-td[iv+1]) * leaf_transmitance + td[iv+1]
  aiv   <- refld - trand*trand/refld
  biv   <- trand/refld
  
  m <- 2
  a[m] <- -aiv
  b[m] <- 1
  c[m] <- -biv
  d[m] <- .super$env$par_dir*tbcum[iv+1] * (1-tb[iv+1]) * (leaf_transmitance - leaf_reflectance*biv)
  
  # Leaf layers, excluding top layer
  for (iv in 2:nlayers) {
    
    # Upward flux
    refld <- (1-td[iv]) * leaf_reflectance
    trand <- (1-td[iv]) * leaf_transmitance + td[iv]
    fiv  <- refld - trand*trand/refld
    eiv  <- trand / refld
    
    m <- m + 1
    a[m] <- -eiv
    b[m] <- 1
    c[m] <- -fiv
    d[m] <- .super$env$par_dir*tbcum[iv] * (1-tb[iv]) * (leaf_reflectance - leaf_transmitance*eiv)
    
    # Downward flux
    refld <- (1-td[iv+1]) * leaf_reflectance
    trand <- (1-td[iv+1]) * leaf_transmitance + td[iv+1]
    aiv  <- refld - trand*trand/refld
    biv  <- trand/refld
    
    m <- m + 1
    a[m] <- -aiv
    b[m] <- 1
    c[m] <- -biv
    d[m] <- .super$env$par_dir*tbcum[iv+1] * (1-tb[iv+1]) * (leaf_transmitance - leaf_reflectance*biv)
  }
  
  # Top canopy layer: upward flux
  iv <- nlayers+1
  refld <- (1-td[iv]) * leaf_reflectance
  trand <- (1-td[iv]) * leaf_transmitance + td[iv]
  fiv  <- refld - trand*trand/refld
  eiv  <- trand/refld
  
  m <- m + 1
  a[m] <- -eiv
  b[m] <- 1
  c[m] <- -fiv
  d[m] <- .super$env$par_dir*tbcum[iv] * (1-tb[iv]) * (leaf_reflectance - leaf_transmitance*eiv)
  
  # Top canopy layer: downward flux
  m <- m + 1
  a[m] <- 0
  b[m] <- 1
  c[m] <- 0
  d[m] <- .super$env$par_diff 
  
  # Solve tridiagonal equations for fluxes
  u <- .$solver_tridiagonal(a, b, c, d, m )
  
  
  # Now copy the solution for diffuse fluxes (u) to the upward (swup) and downward (swdn) fluxes for each layer
  # swup - Upward diffuse solar flux above layer
  # swdn - Downward diffuse solar flux onto layer
  
  # Soil diffuse fluxes
  iv <- 1
  m  <- 1
  swup[iv] <- u[m]
  m  <- m + 1
  swdn[iv] <- u[m]
  
  # Leaf layer diffuse fluxes
  for(iv in 2:(nlayers+1)) {
    m <- m + 1
    swup[iv] <- u[m]
    m <- m + 1
    swdn[iv] <- u[m] 
  }
  
  # Absorbed direct beam and diffuse for ground (soil)
  iv <- 1
  direct[iv]  <- .super$env$par_dir*tbcum[iv] * (1-soil_reflectance_dir)
  diffuse[iv] <- swdn[iv] * (1-soil_reflectance_diff)
  swsoi       <- direct[iv] + diffuse[iv]
  
  # Absorbed direct beam and diffuse for each leaf layer and sum
  swveg    <- 0
  swvegsun <- 0
  swvegsha <- 0
  swleafsun <- swleafsha <- rep(NA, nlayers+1 )
  for(iv in 2:(nlayers+1)) {
    
    # Per unit ground area (W/m2 ground)
    # APW: leaf_scattering assumed constant here 
    direct[iv]  <- .super$env$par_dir*tbcum[iv] * (1-tb[iv]) * (1-.super$state_pars$lscattering)
    diffuse[iv] <- (swdn[iv] + swup[iv-1])      * (1-td[iv]) * (1-.super$state_pars$lscattering)
    
    # Absorbed solar radiation for shaded and sunlit portions of leaf layer
    # per unit ground area
    sun[iv]   <- diffuse[iv] * .super$state$vert$sun$fraction[iv-1] + direct[iv]
    shade[iv] <- diffuse[iv] * .super$state$vert$shade$fraction[iv-1] 
   
    # Convert to per unit sunlit and shaded leaf area 
    swleafsun[iv] <- sun[iv]   / (.super$state$vert$sun$fraction[iv-1] * dLAI[iv])
    swleafsha[iv] <- shade[iv] / (.super$state$vert$shade$fraction[iv-1] * dLAI[iv])
    
    # Sum fluxes over all leaf layers
    swveg <- swveg + (direct[iv] + diffuse[iv])
    swvegsun <- swvegsun + sun[iv]
    swvegsha <- swvegsha + shade[iv]
  }
  
  # Albedo
  incoming  <- .super$env$par_dir + .super$env$par_diff
  reflected <- swup[nlayers+1]
  if(incoming>0) {
    albcan <- reflected/incoming
  } else {
    albcan <- 0
  }
 
  # return calculated values 
 .super$state$vert$sun$apar       <- rev(swleafsun[2:(nlayers+1)])
 .super$state$vert$shade$apar     <- rev(swleafsha[2:(nlayers+1)])
 .super$state$vert$sun$fraction   <- rev(.super$state$vert$sun$fraction)
 .super$state$vert$shade$fraction <- rev(.super$state$vert$shade$fraction)
  
  # Conservation check
  # Total radiation balance: absorbed <- incoming - outgoing
  suminc <- incoming 
  sumref <- albcan * incoming 
  sumabs <- suminc - sumref
  
  err <- sumabs - (swveg + swsoi)
  if(abs(err)>1e-03) {
    warning('NormanRadiation: Total solar conservation error of, ', err )
  }
  
  # Sunlit and shaded absorption
  err <- (swvegsun + swvegsha) - swveg
  if(abs(err)>1e-03) {
    warning('NormanRadiation: Sunlit/shade solar conservation error of, ', err )
  }

  if(.super$cpars$verbose) {
    print(''); print('')
    print(paste('Radiation model for a total LAI of:', cumlai[length(cumlai)]))
    print(''); print('')
    print('dlai:')
    print(dLAI)
    print('sumlai:')
    print(sumlai)
    print(''); print('')
    print('tb (frac direct beam not intercepted in each layer):')
    print(tb)
    print('td (frac diffuse beam not intercepted in each layer):')
    print(td)
    print('')
    print('tbcum (frac toc direct beam incident at top of each layer):')
    print(tbcum)
    print('')
    print('cumlai:')
    print(cumlai)
    print(''); print('')
    print('per unit ground area fluxes')
    print('Up diffuse incident at upper layer boundary (Iup, u[odd]:')
    print(swup)
    print('Down diffuse incident at upper layer boundary (Idown, u[even]:')
    print(swdn)
    print('Abs direct:')
    print(direct)
    print(sum(direct))
    print('Abs diffuse:')
    print(diffuse)
    print(sum(diffuse))
    print('Abs sun:')
    print(sun)
    print(paste('Sum:',sum(sun,na.rm=T)))
    print('Abs shade:')
    print(shade)
    print(paste('Shade:',sum(shade,na.rm=T)))
    
    print(''); print('') 
    print('Leaf PARsun:') 
    print(swleafsun[1:(nlayers+1)])
    print('Leaf PARsha:') 
    print(swleafsha[1:(nlayers+1)])
    print('Leaf fsun:') 
    print(rev(.super$state$vert$sun$fraction[1:(nlayers)]))
    print('Leaf fsha:') 
    print(rev(.super$state$vert$shade$fraction[1:(nlayers)]))
    
    print('')
    print('APARveg, APARsoil, reflected:') 
    print(paste(swveg, swsoi, sumref, sep=', ')) 
  }
}



# calculate parameters for RT / light interception
################################

# Beer
f_pars_init <- function(.) {
  # calculate parameters for Beer's Law canopy light scaling
  
  # extinction coefficents of direct diffuse radiation, assuming leaves are optically black   
  # these account for both canopy clumping and solar zenith angle, also informed by Bodin & Franklin 2012 GMD  
  # APW - a consistent approach would consider clumping where k_dir & k_diff are calculated 
  .super$state_pars$G_dir[]       <- .super$pars$G * .super$pars$can_clump
  .super$state_pars$k_dir[]       <- .super$state_pars$G_dir / cos(.super$env$zenith)

  # this is the k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  .super$state_pars$k_diff[]      <- .super$state_pars$G_dir / (2/pi)

  # transmitance equals reflectance
  # APW - cahange this to the sum of trans and refl 
  .super$state_pars$lscattering[] <- 2 * .super$pars$leaf_reflectance

  # adjustment of k to account for scattering, i.e. that leaves are not optically black
  .super$state_pars$m[]           <- (1.0-.super$state_pars$lscattering)^0.5
  .super$state_pars$k_dirprime[]  <- .super$state_pars$m * .super$state_pars$k_dir
  .super$state_pars$k_diffprime[] <- .super$state_pars$m * .super$state_pars$k_diff

  # albedo
  # APW - this is calculated with clumping applied to k, but below clumping is aplied afterwards
  # APW - which is correct? Generally the combination of G(z) and clumping is inconsistent, needs unified 
  # APW - Bonan includes clumping in the calculatiion of albedo, but not in every use of kb/kd: sp_14_03.m
  # APW - maybe best to be explicit about clumping everywhere that it is used
  .$albedo()

  # calculate Vcmax0 and extinction coefficient for Vcmax
  .super$pars$vcmax$layer0[]      <- .$fns$layer0.vcmax(var='vcmax')
  .super$state_pars$k_vcmax[]     <- .$fns$k.vcmax(var='vcmax')
}


# Goudriaan & Norman
f_pars_init_full <- function(.) {
  # calculate parameters for Goudriaan's canopy light scaling
 
  # zenith angles (radians) for numerical approximation of diffuse light extinction & albedo
  dz <- 1/.super$pars$nz
  .super$state_pars$zi       <- zi <- seq(dz/2,1,dz) * pi/2
  .super$state_pars$delta_zi <- delta_zi <- sum(.super$state_pars$zi[1:2] * c(-1,1))  
  
  # extinction coefficents of direct & diffuse radiation, assuming leaves are optically black   
  .super$state_pars$G_dir[]       <- .$gz(z=.super$env$zenith) 
  .super$state_pars$k_dir[]       <- .super$state_pars$G_dir / cos(.super$env$zenith)
  # Prevent large k_dir at low sun angle
  .super$state_pars$k_dir[]       <- min(.super$state_pars$k_dir, 20)
  #.super$state_pars$k_dir_zi      <- .$gz(z=zi) / cos(zi)
  # k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  #numint                          <- sin(zi)*cos(zi)*delta_zi * exp(-.super$state_pars$k_dir_zi * .super$pars$can_clump * .super$env$lai) 
  #.super$state_pars$k_diff[]      <- -log(2*sum(numint)) / (.super$env$lai * .super$pars$can_clump)
  #.super$state_pars$k_dir_zi      <- .$gz(z=zi) / cos(zi)
  # k_dir assuming light is coming from all points of the hemisphere (i.e. solar zenith angle between 0 and pi/2)
  # APW: this will fail if chi_l is a vector -- needs a catch
  # APW: does the simpler calculation in pars_init provide the same result?
  .super$state_pars$k_dir_zi      <- .$gz(z=zi) / cos(zi)
  numint                          <- .$td(.super$state_pars$zi, dLAI=.super$env$lai) 
  .super$state_pars$k_diff[]      <- -log(2*sum(numint)*.super$state_pars$delta_zi) / (.super$env$lai*.super$pars$can_clump)

  # scattering transmitance + reflectance
  .super$state_pars$lscattering[] <- .super$pars$leaf_reflectance + .super$pars$leaf_transmitance
  .super$state_pars$m[]           <- (1.0-.super$state_pars$lscattering)^0.5
  .super$state_pars$k_dirprime[]  <- .super$state_pars$m * .super$state_pars$k_dir
  .super$state_pars$k_diffprime[] <- .super$state_pars$m * .super$state_pars$k_diff

  # adjustment of k_Xprime 
  # - to account for clumping 
  .super$state_pars$k_dirprime[]  <- .super$state_pars$k_dirprime[]  * .super$pars$can_clump
  # APW - seems like k_diff is already calculated with clumping considered 
  .super$state_pars$k_diffprime[] <- .super$state_pars$k_diffprime[] * .super$pars$can_clump

  # albedo
  .$albedo()

  # adjustment of k_X 
  # - to account for clumping,  
  .super$state_pars$k_dir[]       <- .super$state_pars$k_dir[]  * .super$pars$can_clump
  # APW - seems like k_diff is already calculated with clumping considered 
  .super$state_pars$k_diff[]      <- .super$state_pars$k_diff[] * .super$pars$can_clump

  # calculate Vcmax0 and extinction coefficient for Vcmax
  .super$pars$vcmax$layer0[]      <- .$fns$layer0.vcmax(var='vcmax')
  .super$state_pars$k_vcmax[]     <- .$fns$k.vcmax(var='vcmax')
}

# calculate G(z), Bonan 2019
f_gz_rossgoudriaan <- function(., chi_l=.super$pars$chi_l, z=.super$env$zenith ) {
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

f_gz_constant <- function(., ... ) {
  .super$pars$G
}

# albedo functions for Goudriaan's model
f_albedo_goudriaan <- function(.) {
  # after Wang 2003
  # informed by Bonan 2019, Climate Change and Terrestrial Ecosystem Modeling
  
  # albedo for canopies of infinite LAI
  # - with horizontal leaves
  .super$state_pars$alb_h        <- (1-.super$state_pars$m) / (1+.super$state_pars$m)
  # - adjust direct beam albedo for leaf angle distribution
  .super$state_pars$alb_dir_can  <- .super$state_pars$alb_h*2*.super$state_pars$k_dir / (.super$state_pars$k_dir + .super$state_pars$k_diff)
  # - adjust diffuse beam albedo for leaf angle distribution
  .$diffalbedo()

  # adjust for finite LAI
  .super$state_pars$alb_dir  <- .super$state_pars$alb_dir_can  + (.super$pars$soil_reflectance_dir-.super$state_pars$alb_dir_can) * 
                                  exp(-2*.super$state_pars$k_dirprime  * .super$env$lai)
  .super$state_pars$alb_diff <- .super$state_pars$alb_diff_can + (.super$pars$soil_reflectance_diff-.super$state_pars$alb_diff_can) * 
                                  exp(-2*.super$state_pars$k_diffprime * .super$env$lai)
}

# albedo for diffuse radiation
f_diffalbedo_goudriaan <- function(., zi=.$state_pars$zi, delta_zi=.$state_pars$delta_zi ) {
  # proper solution - assumes random azimuth distribution

  numint <- sin(zi)*cos(zi)*delta_zi * .super$state_pars$alb_h*.super$state_pars$k_dir_zi / (.super$state_pars$k_dir_zi + .super$state_pars$k_diff) 
  .super$state_pars$alb_diff_can  <- 2 * sum(numint) 
} 

f_diffalbedo_approx <- function(.) {
  # assumes G(z) is 0.5 at all zenith angles - incorrect
  
  .super$state_pars$alb_diff_can <- 4 * .super$state_pars$G_dir * .super$state_pars$alb_h * ( .super$state_pars$G_dir *
                                      (log(.super$state_pars$G_dir)-log(.super$state_pars$G_dir+.super$state_pars$k_diff)) / .super$state_pars$k_diff^2 + 
                                      1/.super$state_pars$k_diff )
} 



# partition par
###############################

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



# Scaling of traits/env through the canopy
# - for use with a multilayer phototsynthesis scheme
################################

# function to apply scaling functions to leaf traits
f_traits_scale <- function(., l=.super$state_pars$ca_calc_points ) {

  .super$state$vert$leaf$leaf.leafN_area[]   <- .$scale.n(l, var='n' )
  .super$state$vert$leaf$leaf.atref.vcmax[]  <- .$scale.vcmax(l, var='vcmax' )
  leaf_vars <- c('leaf.leafN_area', 'leaf.atref.vcmax' )

  # APW: if scale-jmax etc are not two-layer this will be inconsistent
  if(.super$fnames$scale$vcmax=='f_scale_two_layer') {
    .super$state$vert$leaf$leaf.atref.jmax[] <- .$scale.jmax(l, var='jmax' )
    .super$state$vert$leaf$leaf.f[]          <- .$scale.f(l, var='f' )
    .super$state$vert$leaf$leaf.g1_medlyn[]  <- .$scale.g1(l, var='g1' )
    leaf_vars <- c(leaf_vars, 'leaf.atref.jmax', 'leaf.f', 'leaf.g1_medlyn' )
  }
  leaf_vars
}


# Uniform scaling
f_scale_uniform <- function(., l, var, vlist='env' ) {
  .super[[vlist]][[var]]
}

f_scale_uniform_layer0 <- function(., l, var, vlist='pars' ) {
  .super[[vlist]][[var]]$layer0
}

# Use Beer's Law to scale leaf traits/env through the canopy 
# - based on canopy layer 0 value 
f_scale_beerslaw <- function(., l, var, vlist='pars' ) {
  # for use with a multilayer phototsynthesis scheme
  # vcmax should probably be an env variable
  
  .super$state_pars$k[[var]] <- .[[paste0('k.',var)]](., var=var )
  .super[[vlist]][[var]]$layer0 * exp(-.super$state_pars$k[[var]] * l ) 
}

# Use Beer's Law to to scale leaf traits/env through the canopy 
# - based on total canopy value 
f_scale_beerslaw_total <- function(., l, var, vlist='pars' ) {
  .super$state_pars$k[var] <- .[[paste0('k.',var)]]()
  .super[[vlist]][[var]]$total * exp(-.super$state_pars$k[[var]]*l) /  sum(exp(-.super$state_pars$k[[var]]*1:.super$env$lai))  
}


# two values, high & low LAI
f_scale_two_layer <- function(., l, var, vlist='pars' ) {
  can_vector <- l
  can_vector[] <- .super[[vlist]][[var]]$layer2
  can_vector[l<.super$pars$can_lai_upper] <- .super[[vlist]][[var]]$layer1
  can_vector
}


# set zero layer value
# - tell leaf object to read var value from pars list, i.e. as a constant that is assigned by the canopy object
# APW: how to handle the leaf functions when either Vcmax or N is passed to the leaf model? 
f_layer0_constant <- function(., var ) {
  fnames_update        <- paste('f', var, 'constant', sep='_' )
  names(fnames_update) <- c(paste0('leaf.',var)) 
  .super$configure(.=.super, vlist='fnames', df=fnames_update )

  .super$pars[[var]]$layer0
}


# set k for exponential/Beer's law scaling
f_k_constant  <- function(., var ) .super$pars$k[[var]]

f_k_kdirect   <- function(., var ) .super$state_pars$k_dir

f_k_exponential_exponent_flayer0 <- function(., var ) {
  exp( .super$pars$k_expa[[var]] + .super$pars$k_expb[[var]] * .super$state[[var]]$layer0 )
}



## Nitrogen Scaling
#################################
#
#f_scale_n_CLMuniform <- function(.,l) {
#  rep( (.super$state$mass_a/ceiling(.super$env$lai)) / .super$state$C_to_N, length(l) ) 
#}
#
## Use Beer's Law to scale leaf N through the canopy
#f_scale_n_beerslaw <- function(.,l) {
#  # for use with a multilayer phototsynthesis scheme
#  
#  .super$state$totalN * exp(-.super$state_pars$k_dir*l) /  sum(exp(-.super$state_pars$k_dir*1:.super$env$lai))  
#}
#
### Use Beer's Law to scale leaf vcmax through the canopy 
##f_scale_vcmax_beerslaw <- function(.,l) {
##  # for use with a multilayer phototsynthesis scheme
##  
##  .$state$vcmax0 * exp(-.$state_pars$k_vcmax*l) 
##}
##
##f_scale_vcmax_uniform <- function(., layers ) {
##  rep(.$state$vcmax0, length(layers) ) 
##}
##
##f_k_vcmax_constant  <- function(.) .$pars$k_vcmax
##
##f_k_vcmax_lloyd2012 <- function(.) {
##   exp( .$pars$k_vcmax_expa + .$pars$k_vcmax_expb*.$state$vcmax0 )
##}
#
#
#
## Vcmax Scaling
#################################
#
#f_vcmax0_constant <- function(.) {
#  #.super$leaf$fnames$vcmax <- 'f_vcmax_constant'
#  .super$configure(.=.super, vlist='fnames', df=c(leaf.vcmax='f_vcmax_constant') )
#  .super$pars$vcmax0
#}
#
#
#f_k_vcmax_constant  <- function(.) .super$pars$k_vcmax
#
#
#f_k_vcmax_lloyd2012 <- function(.) {
#   exp( .super$pars$k_vcmax_expa + .super$pars$k_vcmax_expb*.super$state$vcmax0 )
#}
#
#
## Use Beer's Law to scale leaf vcmax through the canopy 
#f_scale_vcmax_beerslaw <- function(., l, var ) {
#  # for use with a multilayer phototsynthesis scheme
#  # currently var is a dummy argument for compatibility with two_layer
#  
#  .super$state$vcmax0 * exp(-.super$state_pars$k_vcmax*l) 
#}
#
#
#f_scale_vcmax_uniform <- function(., layers, var ) {
#  # currently var is a dummy argument for compatibility with two_layer
#  rep(.super$state$vcmax0, length(layers) ) 
#}
#
#
#
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
