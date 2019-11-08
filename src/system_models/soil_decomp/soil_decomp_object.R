################################
#
# MAAT soil_decomp model object 
#
# AWalker, Matt Craig, October 2019 
#
################################

library(proto)
source('soil_decomp_functions.R')
source('soil_decomp_SoilR_functions.R')
source('soil_decomp_system_functions.R')



# soil_decomp OBJECT
###############################################################################

# use generic soil_decomp
setwd('..')
source('generic_model_object.R')
soil_decomp_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('soil_decomp')



# assign object functions
###########################################################################
soil_decomp_object$name <- 'soil_decomp'


# function to configure unique elements of the object
# - adds functions to fns that are not in fnames
# - or functions that are derivations of other functions, see leaf_object.R for an example case 
####################################
soil_decomp_object$configure_unique <- function(., init=F, flist=NULL ) {
  if(init) {
    source('../../functions/packagemod_functions.R')
    .$fns$plsoda         <- plsoda
    .$fns$inputrates     <- f_inputrates
    .$fns$DotO           <- f_DotOi
    .$fns$transfermatrix <- f_transfermatrixi
    .$fns$solver_func    <- f_solver_func
  }

  if(any(names(flist)=='xyz')) {
   .$fns$xyz_fe <- get(paste0(.$fnames$xyz,'_fe'), pos=1 )
  }
}


# assign unique run & init functions
####################################

soil_decomp_object$init1 <- init_state

soil_decomp_object$init  <- function(.) {
  .$init1()
  .$state$cpools   = matrix(unlist(.$pars$cstate0), ncol=1 )
  .$state$cpools_n = length(.$pars$cstate0)
}



# assign object variables 
###########################################################################

# function names
####################################
soil_decomp_object$fnames <- list(
  sys   = 'f_sys_npools',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_MM_microbe',
    d2 = 'f_decomp_lin',
    d3 = 'f_decomp_MM_microbe'
  ),
  
  # transfer list
  transfer = list(
    t1_to_2 = 'f_transfer_cue',
    #t1_to_3 = 'f_transfer_zero',
    #t2_to_1 = 'f_transfer_zero',
    t2_to_3 = 'f_transfer_cue',
    #t3_to_1 = 'f_transfer_zero',
    t3_to_2 = 'f_transfer_cue'
  )
)


# environment
####################################
soil_decomp_object$env <- list(
  litter = 3.2,
  temp   = 10,
  swc    = 0.3
)


# state
####################################
soil_decomp_object$state <- list(
  cpools   = matrix(1:3, ncol=1 ),
  cpools_n = 3 
)


# state parameters (i.e. calculated parameters)
####################################
soil_decomp_object$state_pars <- list(
  solver_out = matrix(1)
)


# parameters
####################################
soil_decomp_object$pars <- list(
  cstate0 = list(
    c1 = 0.6,
    c2 = 0.1,
    c3 = 0.1
  ),
  ks      = 1.8e-05,
  kb      = 0.007,
  Km      = 900,
  r       = 0.6,
  Af      = 1,
  silt    = 0.2,
  clay    = 0.2,
  cuec1   = 0.47,       # currently, all cue values are the same (at the value in MEND), might need a different cue_max for density dependent function... 
  h       = 0.566,      # humification constant
  cuec3   = 0.47,     
  vmax1   = 0.2346,     
  vmax3   = 0.0777,   
  km1     = 101,       
  km3     = 250,       
  k23     = 0.00672,    # microbial turnover constant
  mbcmax  = 2,          # microbial biomass max value
  beta    = 1.5,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  maommax = 26.725,     # max maom capacity (calculated using Hassink formula assuming 15% clay)
  cue = list(
    cue1 = 0.47,       # currently, all cue values are the same (at the value in MEND), might need a different cue_max for density dependent function... 
    cue2 = 0.566,      # humification constant
    cue3 = 0.47   
  ),  
  vmax = list(   
    vmax1 = 0.2346,     
    vmax2 = 0.2346,     
    vmax3 = 0.0777   
  ),       
  km = list(   
    km1 = 101,       
    km2 = 101,       
    km3 = 250
  ),       
  k = list(   
    k1 = 0.007,       
    k2 = 0.00672,    # microbial turnover constant
    k3 = 0.006
  ),
  poolmax = list(       
    pmax1 = 2,       # POM value
    pmax2 = 2,       # microbial biomass max value
    pmax3 = 26.725   # max maom capacity (calculated using Hassink formula assuming 15% clay)
  )
)


# run control parameters
####################################
soil_decomp_object$cpars <- list(
  verbose  = F,          # write diagnostic output during runtime 
  cverbose = F,          # write diagnostic output from configure function 
  output   = 'run'       # type of output from run function
)



# output functions
#######################################################################        

f_output_soil_decomp_run <- function(.) {
  unlist(.$state)
}

f_output_soil_decomp_state <- function(.) {
  unlist(.$state)
}

f_output_soil_decomp_full <- function(.) {
  c(unlist(.$state),unlist(.$state_pars))
}



# test functions
#######################################################################        

soil_decomp_object$.test <- function(., verbose=F, metdf=F, litter=3.2, ntimes=100 ) {

  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test() # if only used in test functions should begin with a .

  if(metdf) {
    .$dataf       <- list()
    if(length(litter)==1) litter <- rep(litter, ntimes )   
    .$dataf$metdf <- as.data.frame(matrix(litter, ncol=1 ))
    #colnames(.$dataf$metdf) <- 'soil_decomp.litter'  
    names(.$dataf$metdf) <- c('soil_decomp.litter')  
    .$dataf$lm    <- length(.$dataf$metdf[,1])
    .$dataf$mout  <- .$output()
    .$run_met()
  } else {
    .$env$litter  <- litter
    .$run()
  }
}


soil_decomp_object$.test_3pool <- function(., verbose=F, metdf=F, litter=0.00384, ntimes=36500, time=T ) {

  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test() # if only used in test functions should begin with a .

  # initialise boundary data 
  .$dataf       <- list()
  if(length(litter)==1) litter <- rep(litter, ntimes )   
  .$dataf$metdf <- as.data.frame(matrix(litter, ncol=1 ))
  #colnames(.$dataf$metdf) <- 'soil_decomp.litter'  
  names(.$dataf$metdf) <- c('soil_decomp.litter')  
  .$dataf$lm    <- length(.$dataf$metdf[,1])
  .$dataf$mout  <- .$output()


  ### Run models
  olist <- list()
  # run default no saturation or DD model
  olist$noSaturation  <- .$run_met()

  # saturating MAOM
  .$fnames$transfer$t2_to_3 <- 'f_transfer_cue_sat'
  .$configure_test() 
  olist$MaomMax       <- .$run_met()

  # denisty dependent microbial turnover 
  .$fnames$transfer$t2_to_3 <- 'f_transfer_cue'
  .$fnames$decomp$d2        <- 'f_decomp_dd'
  .$configure_test() 
  olist$DDturnover    <- .$run_met()

  # denisty dependent microbial cue 
  .$fnames$transfer$t1_to_2 <- 'f_transfer_cue_sat'
  .$fnames$transfer$t3_to_2 <- 'f_transfer_cue_sat'
  .$fnames$decomp$d2        <- 'f_decomp_lin'
  .$configure_test() 
  olist$DDcue         <- .$run_met()

  # denisty dependent microbial turnover and cue 
  .$fnames$decomp$d2        <- 'f_decomp_dd'
  .$configure_test() 
  olist$DDturnover.DDcue <- .$run_met()

  # denisty dependent microbial turnover and cue and MAOM saturation 
  .$fnames$transfer$t2_to_3 <- 'f_transfer_cue_sat'
  .$configure_test() 
  olist$DDturnover.DDcue.MaomMax <- .$run_met()


  # plotting functions
  thp_plot_time <- function(mod) {
    ylab <- expression('Pool C mass ['*gC*' '*m^-2*']')
    matplot(1:dim(mod)[1], mod[,1:3], type='l', ylab=ylab, xlab='Days', lty=1,
            ylim=c(0,max(mod)*1.2), col=1:3,main=deparse(substitute(mod)) )
    legend('topleft', c('POM','MB','MAOM'), lty=1, col=1:3, bty='n')
  }
 
  thp_plot_MBC <- function(mod) {
    matplot(mod[,2], mod[,c(1,3)], type='l', ylab=ylab, xlab='MB C mass', lty=1,
            xlim=c(0,max(mod[,2])),
            ylim=c(0,max(mod)*1.2), col=1:2, main=deparse(substitute(mod)) )
    legend('topleft', c('POM','MAOM'), lty=1, col=1:2, bty='n')
  }
  
  par(mfrow = c(2,3))
  if(time) { 
    ## plotting versus time
    thp_plot_time(olist$noSaturation)
    thp_plot_time(olist$MaomMax)
    thp_plot_time(olist$DDturnover)
    thp_plot_time(olist$DDcue)
    thp_plot_time(olist$DDturnover.DDcue)
    thp_plot_time(olist$DDturnover.DDcue.MaomMax)
  } else {  
    # plotting pools vs MBC
    thp_plot_MBC(olist$noSaturation)
    thp_plot_MBC(olist$MaomMax)
    thp_plot_MBC(olist$DDturnover)
    thp_plot_MBC(olist$DDcue)
    thp_plot_MBC(olist$DDturnover.DDcue)
    thp_plot_MBC(olist$DDturnover.DDcue.MaomMax)
  }

  olist
}
  


### END ###
