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



# assign unique run function
###########################################################################
canopy_object$run <- function(.) {
  if(.$cpars$verbose) print('canopy run()')
 
  # initialise canopy
  .$state$lai <- .$fns$lai() # this also could point to a higher level plant object  
  .$fns$pars_init()

  # assign canopy environment to leaf environment 
  # any canopy scaling of these variables will overwrite the values written here 
  envss     <- which(names(.$env) %in% names(.$leaf$env) )
  df        <- as.data.frame(.$env[envss])
  names(df) <- paste0('leaf.',names(df))
  .$leaf$configure(vlist='env', df=df ) 
  
  # set leaf absorptance to 1 as all cansys functions (should) account for leaf scattering
  .$leaf$pars$a <- 1.0

  # calculate water status
  .$fns$water_status()

  # calculate diffuse and direct radiation
  .$fns$par_partition()      
  
  # run canopy model
  .$fns$sys()
  
  # output
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
  rt            = 'f_rt_beerslaw_goudriaan',
  scale_n       = 'f_scale_n_CLMuniform',
  scale_vcmax   = 'f_scale_vcmax_beerslaw',
  vcmax0        = 'f_vcmax0_constant',
  k_vcmax       = 'f_k_vcmax_constant',
  scale_ca      = 'f_scale_ca_uniform',
  scale_vpd     = 'f_scale_vpd_uniform',
  lai           = 'f_lai_constant',
  par_partition = 'f_par_partition_spitters',
  water_status  = 'f_water_status_none',
  fwdw          = 'f_fwdw_wtd_lin'
)


# environment
####################################
canopy_object$env <- list(
  temp      = 25,      
  par       = 1000,      
  par_dir   = numeric(1),      
  par_diff  = numeric(1),      
  ca_conc   = 400,
  vpd       = 1,
  clearness = 1,
  zenith    = 0,
  water_td  = numeric(1),
  sphag_h   = numeric(1)
)


# state
####################################
canopy_object$state <- list(
  # External
  lai     = numeric(1), 
  mass_a  = 10,
  C_to_N  = 40,
  totalN  = 7,
  vcmax0  = numeric(1),      
 
  # Calculated state
  # canopy layer vectors
  vert    = list(
    # variable canopy environment etc
    leaf = list( 
      leaf.ca_conc     = numeric(1),
      leaf.vpd         = numeric(1),
      leaf.par         = numeric(1),
      leaf.atref.vcmax = numeric(1), 
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
      lim      = numeric(1)
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
      lim      = numeric(1)
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
  m            = numeric(1),    
  G_dir        = numeric(1),
  k_dir        = numeric(1),
  k_diff       = numeric(1),
  k_dirprime   = numeric(1),
  k_diffprime  = numeric(1),
  k_vcmax      = numeric(1),
  lscattering  = numeric(1),
  alb_dir      = numeric(1),
  alb_diff     = numeric(1),
  alb_dir_can  = numeric(1),
  alb_diff_can = numeric(1)
)


# parameters
####################################
canopy_object$pars   <- list(
  layers           = 10,
  lai              = 10,
  lai_max          = 4,
  lai_curve        = 0.5,
  leaf_cores       = 1,
  G                = 0.5,     # light extinction coefficient assuming leaves are black bodies and randomly distributed horizontally, 0.5 assumes random or spherical leaf orientation, 1.5 for Sphagnum Williams & Flannagan, 1998
  can_clump        = 1,       # canopy clumping coefficient, 1 - random horizontal distribution, leaves become more clumped as coefficient goes towards zero.
  k_layer          = 0.5,     # for multilayer canopy, where in the layer to calculate physiology, 0 - bottom, 0.5 - midway, 1 - top; not the correct solution to the simplifying assumption of Beer's law (Wang 2003) 
  alb_soil         = 0.15,    # soil albedo
  leaf_reflectance = 0.075,   # leaf reflectance
  vcmax0           = 35,      # vcmax at extreme top of canopy
  k_vcmax          = 0.2,     # scaling exponent for vcmax through canopy
  k_vcmax_expa     = -2.43,   # intercept parameter in exponnent to calculate scaling exponent for vcmax through canopy
  k_vcmax_expb     = 9.63e-3, # slope parameter in exponnent to calculate scaling exponent for vcmax through canopy
  fwdw_wl_slope    = -0.022,  # delta sphagnum fwdw ratio per mm of decrease in water level      (mm-1), currently from Adkinson & Humpfries 2010, Rydin 1985 has similar intercept but slope seems closer to -0.6 
  fwdw_wl_sat      = 16,      # sphagnum fwdw ratio at 0 water level, currently from Adkinson & Humpfries 2010     
  fwdw_wl_exp_a    = -0.037,  # decrease in sphagnum fwdw ratio as an exponential f of water level (cm), currently from Strack & Price 2009
  fwdw_wl_exp_b    = 3.254    # decrease in sphagnum fwdw ratio as an exponential f of water level (cm) 
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
  c(A=.$state$integrated$A)
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



# test functions
#######################################################################        

canopy_object$.test <- function(., verbose=T, 
                                canopy.par=2000, canopy.ca_conc=400, canopy.lai=6 ) {
  
  # Build, assign fnames, configure 
  .$build(switches=c(F,verbose,F))
  .$leaf$cpars$verbose  <- F
  
  .$env$par        <- canopy.par
  .$env$ca_conc    <- canopy.ca_conc
  .$pars$lai       <- canopy.lai
  .$state$mass_a   <- 175
  .$state$C_to_N   <- 40
  
  .$run()
}


canopy_object$.test_aca <- function(., verbose=F, cverbose=F, 
                                    canopy.par=c(100,1000), canopy.ca_conc=seq(50,1200,50),
                                    leaf.rs = 'f_r_zero', canopy.rt = 'f_rt_beerslaw_goudriaan' ) {
  
  # Build, assign fnames, configure 
  .$build(switches=c(F,verbose,cverbose))
  if(verbose) str.proto(canopy_object)
  .$leaf$cpars$verbose  <- F
  
  .$pars$lai       <- 10
  .$state$mass_a   <- 175
  .$state$C_to_N   <- 40
  
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
