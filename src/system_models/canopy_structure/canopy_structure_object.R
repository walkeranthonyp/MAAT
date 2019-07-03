################################
#
# Canopy structure object 
# 
# AWalker Jun 2019 
#
################################

library(proto)

source('canopy_structure_functions.R')
source('canopy_structure_system_functions.R')



# CANOPY STRUCTURE OBJECT
###############################################################################

# use generic template
setwd('..')
source('generic_model_object.R')
canopy_structure_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('canopy_structure')



# assign object functions
###########################################################################
canopy_structure_object$name              <- 'canopy_structure'
canopy_structure_object$child_list        <- list('canopy') 
canopy_structure_object$build             <- build
canopy_structure_object$build_child       <- build_child  
canopy_structure_object$configure         <- configure
canopy_structure_object$configure_sublist <- configure_sublist  
canopy_structure_object$configure_child   <- configure_child  
canopy_structure_object$run_met           <- run_met  



# assign unique run function
###########################################################################
canopy_structure_object$run <- function(.) {
  if(.$cpars$verbose) print('canopy_structure run()')

  # assign canopy_structure environment to canopy environment 
  # any canopy_structure scaling of these variables will overwrite the values written here 
  envss     <- which(names(.$env) %in% names(.$canopy$env) )
  df        <- as.data.frame(.$env[envss])
  names(df) <- paste0('canopy.',names(df))
  .$canopy$configure(vlist='env', df=df ) 

  # run canopy_structure model
  .$fns$sys()
  
  # output
  .$output()      
}


# initialise the number of layers in the canopy_structure
canopy_structure_object$init_vert <- function(.,l) {
  .$state$vert$young  <- lapply(.$state$vert$young,  function(v,leng) numeric(leng), leng=l )
  .$state$vert$mature <- lapply(.$state$vert$mature, function(v,leng) numeric(leng), leng=l )
  .$state$vert$old    <- lapply(.$state$vert$old,    function(v,leng) numeric(leng), leng=l )
  .$state$vert$layer  <- lapply(.$state$vert$layer,  function(v,leng) numeric(leng), leng=l )
}


# function to run the canopy within the canopy_structure
canopy_structure_object$run_canopy <- function(.,ii,df){
  # This wrapper function is called from an (v/l)apply function to run over each canopy in the canopy_structure
  # assumes that each row of the dataframe are independent and non-sequential
  
  #.$canopy$configure(vlist='env',   df=df[ii,] )
  #.$canopy$configure(vlist='state', df=df[ii,] )
  .$canopy$configure(vlist='pars', df=df[ii,] )
  
  # run canopy
  .$canopy$run()        
}


# function to configure unique elements of the object 
####################################
canopy_structure_object$configure_unique <- NULL 



# assign object variables 
###########################################################################

# function names
####################################
canopy_structure_object$fnames <- list(
  sys            = 'f_sys_leafdem',
  lai            = 'f_lai_leafdem_env',
  leafdem_upper  = 'f_leafdem_upper_wu2017',
  leafdem_lower  = 'f_leafdem_lower_wu2017',
  leafdem_vcmax0 = 'f_leafdem_vcmax0_constant',
  water_status   = 'f_water_status_none',
  fwdw           = 'f_fwdw_wtd_lin'
)


# environment
####################################
canopy_structure_object$env <- list(
  temp       = 25,      
  par        = 1000,      
  ca_conc    = 400,
  vpd        = 1,
  zenith     = 0,
  lai_young  = 1, 
  lai_mature = 2,
  lai_old    = 3,
  water_td   = numeric(1),
  sphag_h    = numeric(1)
)


# state
####################################
canopy_structure_object$state <- list(
  # External
  lai            = numeric(1),    
  lai_young      = numeric(1),    
  lai_mature     = numeric(1),    
  lai_old        = numeric(1),    
  upper_can_prop = numeric(3), 
  lower_can_prop = numeric(3), 
  
  # leaf demography associated traits
  leafdem_traits <- list(
    vcmax0 = numeric(3)
  ),
 
  # Calculated state
  # canopy_structure layer vectors
  vert    = list(
    # canopy_structure leaf demography 
    young = list( 
      apar    = numeric(1),        # canopy absorbed PAR
      A       = numeric(1),
      rd      = numeric(1),
      Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
      Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
      cb      = numeric(1),
      ci      = numeric(1),
      cc      = numeric(1),
      g       = numeric(1),  
      gb      = numeric(1),
      gs      = numeric(1),
      gi      = numeric(1)
    ),
    mature = list( 
      apar    = numeric(1),        # canopy absorbed PAR
      A       = numeric(1),
      rd      = numeric(1),
      Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
      Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
      cb      = numeric(1),
      ci      = numeric(1),
      cc      = numeric(1),
      g       = numeric(1),  
      gb      = numeric(1),
      gs      = numeric(1),
      gi      = numeric(1)
    ),
    old = list( 
      apar    = numeric(1),        # canopy absorbed PAR
      A       = numeric(1),
      rd      = numeric(1),
      Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
      Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
      cb      = numeric(1),
      ci      = numeric(1),
      cc      = numeric(1),
      g       = numeric(1),  
      gb      = numeric(1),
      gs      = numeric(1),
      gi      = numeric(1)
    ),
    layer = list( 
      apar    = numeric(1),        # canopy absorbed PAR
      A       = numeric(1),
      rd      = numeric(1),
      Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
      Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
      Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
      cb      = numeric(1),
      ci      = numeric(1),
      cc      = numeric(1),
      g       = numeric(1),  
      gb      = numeric(1),
      gs      = numeric(1),
      gi      = numeric(1)
    )
  ),
  
  # integrated canopy values
  integrated = list(
    apar    = numeric(1),        # canopy absorbed PAR
    A       = numeric(1),        # canopy assimilation rate                         (umol m-2s-1)
    rd      = numeric(1),        # canopy respiration rate                          (umol m-2s-1)        
    Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
    Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
    Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
    cb      = numeric(1),        # canopy mean boundary layer CO2                   (Pa)
    ci      = numeric(1),        # canopy mean canopy internal CO2                  (Pa) 
    cc      = numeric(1),        # canopy mean chloroplast CO2                      (Pa)
    g       = numeric(1),        # canopy conductance                               (mol H2O m-2 s-1 )
    gb      = numeric(1),        # canopy leaf boundary conductance                 (mol H2O m-2 s-1 )
    gs      = numeric(1),        # canopy stomatal conductance                      (mol H2O m-2 s-1 )
    gi      = numeric(1)         # canopy leaf internal conductance                 (mol CO2 m-2 s-1 ) 
  )
)


# state parameters (i.e. calculated parameters)
####################################
canopy_structure_object$state_pars <- list(
  alb_dir_can  = numeric(1),
  alb_diff_can = numeric(1)
)


# parameters
####################################
canopy_structure_object$pars   <- list(
  canopy_cores     = 1,
  layers           = 30,
  lai              = 10,
  lai_upper        = 2.5,
  lai_ftop         = 0.7,
  vcmax0 = list(
    young  = 7,
    mature = 35,
    old    = 20
  ),
  lai_max          = 4,
  lai_curve        = 0.5,
  fwdw_wl_slope    = -0.022, # delta sphagnum fwdw ratio per mm of decrease in water level      (mm-1), currently from Adkinson & Humpfries 2010, Rydin 1985 has similar intercept but slope seems closer to -0.6 
  fwdw_wl_sat      = 16,     # sphagnum fwdw ratio at 0 water level, currently from Adkinson & Humpfries 2010     
  fwdw_wl_exp_a    = -0.037, # decrease in sphagnum fwdw ratio as an exponential f of water level (cm), currently from Strack & Price 2009
  fwdw_wl_exp_b    = 3.254   # decrease in sphagnum fwdw ratio as an exponential f of water level (cm) 
)


# run control parameters
####################################
canopy_structure_object$cpars <- list(
  verbose       = F,          # write diagnostic output during runtime 
  cverbose      = F,          # write configuration output during runtime 
  output        = 'run'       # type of output from run function
)



# output functions
#######################################################################        

canopy_structure_object$output <- function(.){
  if(.$cpars$output=='run') {
    c(A=.$state$integrated$A, gs=.$state$integrated$gs, rd=.$state$integrated$rd)
    
  } else if(.$cpars$output=='canopy') {
    c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
      gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd )
    
  } else if(.$cpars$output=='all_lim') {
    c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
      gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd,
      Acg_lim=.$state$integrated$Acg_lim, 
      Ajg_lim=.$state$integrated$Ajg_lim, 
      Apg_lim=.$state$integrated$Apg_lim 
    )
    
  } else if(.$cpars$output=='full') {
    c(unlist(.$state$integrated), unlist(.$state_pars) )
  }
}    



# test functions
#######################################################################        

canopy_structure_object$.test <- function(., verbose=T, par=1320, ca_conc=400  ){
  
  # Child Objects
  .$build()

  # parameter settings
  .$cpars$verbose       <- verbose
  .$canopy$cpars$verbose  <- F
  
  .$env$par        <- par 
  .$env$ca_conc    <- ca_conc
  
  .$run()
}


canopy_structure_object$.test_aca <- function(., verbose=F, verbose_loop=F, 
                                              canopy_structure.par=c(100,1000), 
                                              canopy_structure.ca_conc=seq(50,1200,50),
                                              rs = 'f_r_zero' ) {
  
  # Child Objects
  .$build()
  .$canopy$leaf$fnames$rs <- rs
  .$canopy$leaf$configure_test

  .$cpars$verbose         <- verbose
  .$canopy$cpars$verbose  <- F
  
  if(verbose) str.proto(canopy_structure_object)
 
  # generate met data 
  .$dataf       <- list()
  .$dataf$met   <- expand.grid(mget(c('canopy_structure.ca_conc','canopy_structure.par')))
  
  # run maat
  .$dataf$out  <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
  
  # run output
  print(cbind(.$dataf$met,.$dataf$out))
  p1 <- 
    if(length(canopy_structure.ca_conc) > 2) 
      xyplot(A~.$dataf$met$canopy_structure.ca_conc|as.factor(.$dataf$met$canopy_structure.par),
             .$dataf$out,abline=0,
             ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),
             xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'))
    else if(length(canopy_structure.par) > 2) 
      xyplot(A~.$dataf$met$canopy_structure.par|as.factor(.$dataf$met$canopy_structure.ca_conc),
             .$dataf$out,abline=0,
             ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),
             xlab=expression('PAR ['*mu*mol*' '*m^-2*s^-1*']'))
  
  print(p1)
}



### END ###
