################################
#
# Canopy object
#
# AWalker Jun 2019
#
################################

source('canopy_system_functions.R')
source('canopy_functions.R')



# CANOPY OBJECT
###############################################################################

# use generic template
setwd('..')
source('generic_model_object.R')
canopy_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('canopy')



# assign object functions
###########################################################################
canopy_object$name              <- 'canopy'
canopy_object$child_list        <- list('leaf')
canopy_object$build_child       <- build_child
canopy_object$configure_child   <- configure_child



# function to configure unique elements of the object
# - adds functions to fns that are not in fnames
# - or functions that are derivations of other functions:
####################################
canopy_object$configure_unique <- function(., init=F, flist=NULL ) {

#  if(any(names(flist)=='rs')) {
#   .$fns$rs_fe <- get(paste0(.$fnames$rs,'_fe'), pos=1 )
#   .$fns$rs_r0 <- get(paste0(.$fnames$rs,'_r0'), pos=1 )
#  }

  if(any(names(flist)=='rt'))
    if(grepl('norman',flist[['rt']])) .$fns$solver_tridiagonal <- f_solver_tridiagonal
}



# assign unique run function
###########################################################################
canopy_object$run <- function(.) {
  if(.$cpars$verbose) print('canopy run()')

  # initialise canopy
  .$fns$pars_init()

  # assign canopy environment to leaf environment
  # any canopy scaling of these variables will overwrite the values written here
  envss     <- which(names(.$env) %in% names(.$leaf$env) )
  df        <- as.data.frame(.$env[envss])
  names(df) <- paste0('leaf.',names(df))
  .$leaf$configure(vlist='env', df=df )

  # set leaf absorptance to 1 as all cansys functions (should) account for leaf scattering
  # APW: shift this to a leaf init XML that configures a leaf for the canopy
  .$leaf$pars$a <- 1.0

  # calculate diffuse and direct radiation
  .$fns$par_partition()

  # run canopy model
  .$fns$sys()

  # output
  .$state$integrated$fapar <- .$state$integrated$apar / .$env$par 
  .$output()
}



# functions unique to object that do not live in fnames/fns, i.e. do not vary ever
# - will these cause trouble if called from within the fns object? probably need to be called as .super
####################################

# initialise the number of layers in the canopy
canopy_object$init_vert <- function(.,l) {
  .$state$vert$leaf  <- lapply(.$state$vert$leaf,  function(v,leng) numeric(leng), leng=l )
  .$state$vert$sun   <- lapply(.$state$vert$sun,   function(v,leng) numeric(leng), leng=l )
  .$state$vert$shade <- lapply(.$state$vert$shade, function(v,leng) numeric(leng), leng=l )
  .$state$vert$layer <- lapply(.$state$vert$layer, function(v,leng) numeric(leng), leng=l )
}


# function to run the leaves within the canopy
canopy_object$run_leaf <- function(., ii, df ) {
  # This wrapper function is called from an (v/l)apply function to run over each leaf in the canopy
  # assumes that each row of the dataframe are independent and non-sequential

  .$leaf$configure(vlist='env',   df=df[ii,] )
  .$leaf$configure(vlist='state', df=df[ii,] )
  .$leaf$configure(vlist='pars',  df=df[ii,] )

  # run leaf
  .$leaf$run()
}



# assign object variables
###########################################################################

# function names
####################################
canopy_object$fnames <- list(
  sys           = 'f_sys_multilayer',
  pars_init     = 'f_pars_init',
  traits_scale  = 'f_traits_scale',
  par_partition = 'f_par_partition_spitters_hourly',
  rt            = 'f_rt_beerslaw_goudriaan',
  gz            = 'f_gz_rossgoudriaan',
  albedo        = 'f_albedo_goudriaan',
  diffalbedo    = 'f_diffalbedo_goudriaan',
  layer0 = list(
    vcmax = 'f_layer0_constant'
  ),
  k = list(
    vcmax = 'f_k_constant',
    n     = 'f_k_kdirect'
  ),
  scale = list(
    n       = 'f_scale_uniform_layer0',
    vcmax   = 'f_scale_beerslaw',
    jmax    = 'f_scale_two_layer',
    f       = 'f_scale_two_layer',
    g1      = 'f_scale_two_layer',
    ca_conc = 'f_scale_uniform',
    vpd     = 'f_scale_uniform'
  )
)

# environment
####################################
canopy_object$env <- list(
  lai       = 6, 
  temp      = 25,
  par       = 1000,
  par_dir   = numeric(1),
  par_diff  = numeric(1),
  ca_conc   = 400,
  vpd       = 1,
  clearness = 1,
  zenith    = 0
)


# state
####################################
canopy_object$state <- list(

  # canopy layer vectors
  vert    = list(
    # variable canopy environment etc
    leaf = list(
      leaf.ca_conc     = numeric(1),
      leaf.vpd         = numeric(1),
      leaf.par         = numeric(1),
      leaf.atref.vcmax = numeric(1), 
      leaf.atref.jmax  = numeric(1), 
      leaf.f           = numeric(1), 
      leaf.g1_medlyn   = numeric(1), 
      leaf.leafN_area  = numeric(1)
    ),
    # variable canopy light & physiology by sun and shade leaves
    sun = list(
      apar     = numeric(1),
      fraction = numeric(1),
      A        = numeric(1),
      rd       = numeric(1),
      ci       = numeric(1),
      cc       = numeric(1),
      gb       = numeric(1),
      gs       = numeric(1),
      gi       = numeric(1),
      g        = numeric(1),
      lim      = numeric(1),
      Acg      = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg      = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
      Apg      = numeric(1)         # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
    ),
    shade = list(
      apar     = numeric(1),
      fraction = numeric(1),
      A        = numeric(1),
      rd       = numeric(1),
      ci       = numeric(1),
      cc       = numeric(1),
      gb       = numeric(1),
      gs       = numeric(1),
      gi       = numeric(1),
      g        = numeric(1),
      lim      = numeric(1),
      Acg      = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg      = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
      Apg      = numeric(1)         # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
    ),
    layer = list(
      apar     = numeric(1),
      A        = numeric(1),
      rd       = numeric(1),
      ci       = numeric(1),
      cc       = numeric(1),
      gb       = numeric(1),
      gs       = numeric(1),
      gi       = numeric(1),
      g        = numeric(1),
      Acg_lim  = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg_lim  = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)
      Apg_lim  = numeric(1)         # assimilation rate of canopy layers Ap limited    (umol m-2s-1)
    )
  ),

  # integrated canopy values
  integrated = list(
    apar       = numeric(1),        # canopy absorbed PAR
    fapar      = numeric(1),        # canopy fraction incoming PAR absorbed 
    A          = numeric(1),        # canopy assimilation rate                         (umol m-2s-1)
    Acg_lim    = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
    Ajg_lim    = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)
    Apg_lim    = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)
    cb         = numeric(1),        # canopy mean boundary layer CO2                   (Pa)
    ci         = numeric(1),        # canopy mean leaf internal CO2                    (Pa)
    cc         = numeric(1),        # canopy mean chloroplast CO2                      (Pa)
    gb         = numeric(1),        # canopy leaf boundary conductance                 (mol H2O m-2 s-1)
    gs         = numeric(1),        # canopy stomatal conductance                      (mol H2O m-2 s-1)
    gi         = numeric(1),        # canopy leaf internal conductance                 (mol CO2 m-2 s-1)
    g          = numeric(1),        # canopy total conductance                         (mol H2O m-2 s-1)
    rd         = numeric(1)         # canopy respiration rate                          (umol m-2s-1)
  )
)


# state parameters (i.e. calculated parameters)
####################################
canopy_object$state_pars <- list(
  lscattering  = numeric(1),        # leaf reflectance + transmitance    
  m            = numeric(1),        # goudriaan's scattering adjustment for non-optically black leaves
  zi           = numeric(1),        # nz zenith angles (radians) for numerical approximation 
  delta_zi     = numeric(1),        # difference between each zi value
  G_dir        = numeric(1),        # mean leaf projected area in direction of beam
  k_dir        = numeric(1),        # extinction coefficient for direct beam
  k_dir_zi     = numeric(1),        # extinction coefficient for direct beam at all zi
  k_diff       = numeric(1),        # extinction coefficient for diffuse beam
  k_dirprime   = numeric(1),        # extinction coefficient for direct beam adjusted for scattering
  k_diffprime  = numeric(1),        # extinction coefficient for diffuse beam adjusted for scattering
  k_vcmax      = numeric(1),        # extinction coefficient for vcmax  
  alb_h        = numeric(1),        # albedo of infinite lai canopy with horizontal leaves  
  alb_dir_can  = numeric(1),        # albedo of infinite lai canopy for direct beam and actual leaf angle distribution  
  alb_diff_can = numeric(1),        # albedo of infinite lai canopy for diffuse beam and actual leaf angle distribution  
  alb_dir      = numeric(1),        # albedo of finite lai canopy for direct beam and actual leaf angle distribution  
  alb_diff     = numeric(1)         # albedo of finite lai canopy for diffuse beam and actual leaf angle distribution  
)


# parameters
####################################
canopy_object$pars   <- list(
  layers                = 10,       # number of layers in multi-layer canopy calculation
  layers_2bigleaf       = 10,       # number of layers in two-bigleaf RT integration 
  nz                    = 9,        # number of discrete solar zenith angles in numerical approximation of diffuse light parameters
  can_lai_upper         = 2.6,      # LAI of upper canopy for two layer scaling scheme 
  leaf_cores            = 1,
  G                     = 0.5,      # light extinction coefficient assuming leaves are black bodies and randomly distributed horizontally, 0.5 assumes random or spherical leaf orientation, 1.5 for Sphagnum Williams & Flannagan, 1998
  chi_l                 = 0.0,      # Ross index indicating departure from spherical leaf angle distribution 
  can_clump             = 1,        # canopy clumping coefficient, 1 - random horizontal distribution, leaves become more clumped as coefficient goes towards zero.
  k_layer               = 0.5,      # for multilayer/ numerical 2bigleaf canopy, where in the layer to calculate APAR & physiology, 0 - bottom, 0.5 - midway, 1 - top; not the correct solution to the simplifying assumption of Beer's law (Wang 2003)
  alb_soil              = 0.15,     # soil albedo
  soil_reflectance_dir  = 0.1,      # soil reflectance for direct radiation (visible) 
  soil_reflectance_diff = 0.1,      # soil reflectance for diffuse radiation (visible) 
  leaf_reflectance      = 0.10,     # leaf reflectance
  leaf_transmitance     = 0.05,     # leaf reflectance
 
  # could be env variables 
  mass_a  = 10,
  C_to_N  = 40,
  
  k = list(                         # scaling exponent for vcmax through canopy
    vcmax  = 0.2
  ),    
  k_expa = list(                    # intercept parameter in exponnent to calculate scaling exponent for vcmax through canopy
    vcmax = -2.43    
  ),
  k_expb = list(                    # slope parameter in exponent to calculate scaling exponent for vcmax through canopy
    vcmax = 9.63e-3   
  ),
  n = list(                         # N area:
    layer0 = 3,                     #   at extreme top of canopy
    layer1 = 2,                     #   in first canopy layer
    layer2 = 1,                     #   in second canopy layer
    total  = 7                      #   canopy total sum
  ),
  vcmax = list(                     # vcmax:
    layer0 = 35,                    #   at extreme top of canopy
    layer1 = 108.7,                 #   in first canopy layer
    layer2 = 67.9                   #   in second canopy layer
  ),
  jmax = list(                      # jmax:
    layer0 = 65,                    #   at extreme top of canopy
    layer1 = 170,                   #   in first canopy layer
    layer2 = 76.3                   #   in second canopy layer
  ),
  f = list(                         # f:
    layer0 = 0.41,                  #   at extreme top of canopy
    layer1 = 0.41,                  #   in first canopy layer
    layer2 = 0.34                   #   in second canopy layer
  ),
  g1 = list(                        # g1 (medlyn):
    layer0 = 3,                     #   at extreme top of canopy
    layer1 = 3.41,                  #   in first canopy layer
    layer2 = 8.82                   #   in second canopy layer
  )
)



# run control parameters
####################################
canopy_object$cpars <- list(
  verbose       = F,          # write diagnostic output during runtime
  cverbose      = F,          # write configuration output during runtime
  output        = 'run'       # type of output from run function
)



# output functions
#######################################################################

f_output_canopy_eval <- f_output_eval

f_output_canopy_state <- function(.) {
  unlist(.$state)
}

f_output_canopy_full <- function(.) {
  unlist(c(.$state, .$state_pars ))
}

f_output_canopy_run <- function(.) {
  c(A=.$state$integrated$A, gs=.$state$integrated$gs, rd=.$state$integrated$rd)
}

f_output_canopy_mcmc <- function(.) {
  unlist(.$state$integrated$A)
}

f_output_canopy_leaf <- function(.) {
  c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci,
    gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd, lim=NA)
}

f_output_canopy_all_lim <- function(.) {
  c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci,
    gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd, lim=NA,
    Acg_lim=.$state$integrated$Acg_lim,
    Ajg_lim=.$state$integrated$Ajg_lim,
    Apg_lim=.$state$integrated$Apg_lim
    )
}

f_output_canopy_canopy_structure <- function(.) {
  vapply(.$state$vert$layer, function(v) v, .$state$vert$layer[[1]] )
}

f_output_canopy_wtc <- function(.) {
  c(A=.$state$integrated$A, 
    gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd,
    Acg_lim=.$state$integrated$Acg_lim, 
    Ajg_lim=.$state$integrated$Ajg_lim, 
    Apg_lim=.$state$integrated$Apg_lim#,
    #Alayer=.$state$vert$layer$A 
    )
}



# test functions
#######################################################################

canopy_object$.test <- function(., verbose=T,
                                canopy.par=2000, canopy.ca_conc=400, 
                                canopy.lai=6, canopy.layers=10,
                                canopy.diffalbedo='f_diffalbedo_approx', 
                                canopy.rt='f_rt_beerslaw_goudriaan'
                                ) {

  # Build, assign fnames, configure
  .$build(switches=c(F,verbose,F))
  .$leaf$cpars$verbose  <- F

  .$env$par        <- canopy.par
  .$env$ca_conc    <- canopy.ca_conc
  .$env$lai        <- canopy.lai
  .$pars$layers    <- canopy.layers
  .$state$mass_a   <- 175
  .$state$C_to_N   <- 40
  
  if(grepl('norman',canopy.rt) | grepl('goudriaan',canopy.rt)) {
    .$fnames$pars_init <- 'f_pars_init_full'
  }
  .$fnames$diffalbedo <- canopy.diffalbedo
  .$fnames$rt         <- canopy.rt
  
  .$configure_test() 
  .$leaf$configure_test() 

  .$run()
}


canopy_object$.test_2bigleaf <- function(., verbose=F,
                                         canopy.par=2000, canopy.ca_conc=400, 
                                         canopy.lai=6, canopy.layers_2bigleaf=100,
                                         canopy.diffalbedo='f_diffalbedo_goudriaan', 
                                         canopy.rt='f_rt_goudriaan'
                                         ) {
         
  # Build, assign fnames, configure
  .$build(switches=c(F,verbose,F))
  .$leaf$cpars$verbose  <- F

  .$env$par        <- canopy.par
  .$env$ca_conc    <- canopy.ca_conc
  .$env$lai        <- canopy.lai
  .$pars$layers_2bigleaf <- canopy.layers_2bigleaf
  
  .$fnames$sys        <- 'f_sys_2bigleaf' 
  .$fnames$pars_init  <- 'f_pars_init_full'
  .$fnames$diffalbedo <- canopy.diffalbedo
  .$fnames$rt         <- canopy.rt
  
  .$configure_test() 
  .$leaf$configure_test() 

  .$run()
}

canopy_object$.test_aca <- function(., verbose=F, cverbose=F,
                                    canopy.par=c(100,1000), canopy.ca_conc=seq(50,1200,50),
                                    canopy.lai=6, canopy.layers=10,
                                    leaf.rs='f_r_zero', canopy.rt='f_rt_beerslaw_goudriaan' ) {

  # Build, assign fnames, configure
  .$build(switches=c(F,verbose,cverbose))
  if(verbose) str.proto(canopy_object)
  .$leaf$cpars$verbose  <- F

  .$env$lai        <- canopy.lai
  .$pars$layers    <- canopy.layers
  .$state$mass_a   <- 175
  .$state$C_to_N   <- 40

  if(grepl('norman',canopy.rt) | grepl('goudriaan',canopy.rt)) {
    .$fnames$pars_init <- 'f_pars_init_full'
  }
  .$fnames$rt      <- canopy.rt
  .$leaf$fnames$rs <- leaf.rs

  .$configure_test() 
  .$leaf$configure_test() 
  
  .$dataf          <- list()
  .$dataf$met      <- t(as.matrix(expand.grid(mget(c('canopy.ca_conc','canopy.par')))))
  .$dataf$lm       <- dim(.$dataf$met)[2]
  .$dataf$mout     <- .$output()
  .$dataf$out      <- .$run_met() 
  .$dataf$out_full <- as.data.frame(cbind(t(.$dataf$met), .$dataf$out ))
  print(.$dataf$out_full)
  
  p1 <- xyplot(A ~ canopy.ca_conc | as.factor(canopy.par), .$dataf$out_full, abline=0,
               ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'), xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'))
  print(p1)
}



### END ###
