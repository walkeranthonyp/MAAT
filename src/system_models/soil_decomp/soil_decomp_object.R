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
source('../../functions/packagemod_functions.R')



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
soil_decomp_object$name      <- 'soil_decomp'

# parameter names that have a value per pool
soil_decomp_object$pool_pars <- c('cstate0', 'cue', 'vmax', 'km', 'k', 'poolmax', 'input_coefs' )



# function to configure unique elements of the object
# - adds functions to fns that are not in fnames
# - or functions that are derivations of other functions, see leaf_object.R for an example case 
####################################
soil_decomp_object$configure_unique <- function(., init=F, flist=NULL ) {
  if(init) NULL 

  if(any(names(flist)=='xyz')) {
   .$fns$xyz_fe <- get(paste0(.$fnames$xyz,'_fe'), pos=1 )
  }
}


# assign unique run & init functions
####################################

soil_decomp_object$init1 <- init_state

soil_decomp_object$init  <- function(.) {
  .$init1()
  .$state$cpools = matrix(unlist(.$pars$cstate0)[1:.$pars$n_pools], ncol=1 )
}



# assign object variables 
###########################################################################

# function names
####################################
soil_decomp_object$fnames <- list(

  sys            = 'f_sys_npools',
  solver         = 'plsoda',
  solver_func    = 'f_solver_func_corpse',
  input          = 'f_input',
  DotO           = 'f_DotO',
  transfermatrix = 'f_transfermatrix',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_rmm_sulman', 
    d2 = 'f_decomp_rmm_sulman',
    d3 = 'f_decomp_rmm_sulman',
    d4 = 'f_micturn_sulman',
    d5 = 'f_zero',
    d6 = 'f_zero',
    d7 = 'f_zero'
  ),
  
  desorp = list(
    ds1 = 'f_zero',
    ds2 = 'f_zero',
    ds3 = 'f_zero',
    ds4 = 'f_zero',
    ds5 = 'f_decomp_lin',
    ds6 = 'f_decomp_lin',
    ds7 = 'f_decomp_lin'
  ),
  
  sorp = list(
    s1 = 'f_decomp_lin',
    s2 = 'f_decomp_lin',
    s3 = 'f_decomp_lin',
    s4 = 'f_zero',
    s5 = 'f_zero',
    s6 = 'f_zero',
    s7 = 'f_zero'
  ),
  
  scor = 'f_scor_sulman',
  wcor = 'f_wcor_sulman',
  
  tcor = list(
    t1 = 'f_tcor_arrhenius',
    t2 = 'f_tcor_arrhenius',
    t3 = 'f_tcor_arrhenius',
    t4 = 'f_zero',
    t5 = 'f_zero',
    t6 = 'f_zero',
    t7 = 'f_zero'
  ),
  
  # transfer list
  transfer = list(
    t1_to_3 = 'f_transfer_cue_remainder',
    t1_to_5 = 'f_transfer_cue',
    t2_to_1 = 'f_transfer_mend21',
    t2_to_5 = 'f_transfer_mend25',
    t2_to_6 = 'f_transfer_mend26',
    t2_to_7 = 'f_transfer_mend27',
    t3_to_5 = 'f_transfer_all',
    t4_to_5 = 'f_transfer_all',
    t5_to_2 = 'f_transfer_mend52',
    t5_to_4 = 'f_transfer_mend54',
    t6_to_5 = 'f_transfer_all',
    t7_to_5 = 'f_transfer_all'
    #######left off here#######
  )
)


# environment
####################################
soil_decomp_object$env <- list(
  litter = .5/365,
  temp   = 20,
  vwc    = .25,
  porosity = 0.5,
  clay = 5
)


# state
####################################
soil_decomp_object$state <- list(
  cpools   = matrix(1:3, ncol=1 ) 
)


# state parameters (i.e. calculated parameters)
####################################
soil_decomp_object$state_pars <- list(
  solver_out = matrix(1)
  
  #scor
  #wcor
  #tcor
)


# parameters
####################################
soil_decomp_object$pars <- list(

  n_pools = 7,          # number of pools in model  
  beta    = 1.5,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  silt    = 0.2,        # soil silt content (proportion)
  clay    = 0.2,        # soil clay content (proportion)
  clayref = 5,
  mr      = 0.00028,     # specific maintenance factor (MEND)
  pep     = 0.01,        #Fraction of mr allocated for production of EP
  pem     = 0.01,        #Fraction of mr allocated for production of EM
  Kads    = 0.006,       #Specific adsorption rate (could make this k for the DOC pool)
  qslope_mayes = 0.4833,
  reftemp = 20,
  R = 8.314472,
  minmic = 1e-3,
  q10 = 2,
  
  ea = list(
    ea1 = 37e3,
    ea2 = 54e3,
    ea3 = 50e3,
    ea4 = NA,
    ea5 = NA,
    ea6 = NA,
    ea7 = NA
  ),
 
  #input coefficients (allocaties proportions of inputs into different pools)
  input_coefs = list(
    input_coef1 = .1,
    input_coef2 = .9,
    input_coef3 = 0,
    input_coef4 = 0,
    input_coef5 = 0,
    input_coef6 = 0,
    input_coef7 = 0
  ),
  
  # initial pool mass for each pool
  cstate0 = list(
    cstate01 = 0.1, 
    cstate02 = 0.1,          
    cstate03 = 0.1,      
    cstate04 = 0.1,       
    cstate05 = 0.1,     
    cstate06 = 0.1,    
    cstate07 = 0.1     
  ),

  # Carbon use or transfer efficiency from pool i to any another 
  # - if this varies by the 'to' pool we need another function / parameters
  # MC: I've created a function cue_remainder that can allocate the remainder to another pool rather than CO2
  cue = list(
    cue1 = 0.6,       
    cue2 = 0.05,       
    cue3 = 0.6,
    cue4 = 0.5,
    cue5 = 0,      
    cue6 = 0,
    cue7 = 0
  ),  

  # max turnover rate per unit microbial biomass for pool i 
  vmax = list(   
    vmax1 = 1500/365,     
    vmax2 = 50/365,     
    vmax3 = 600/365,
    vmax4 = 0,
    vmax5 = 0,
    vmax6 = 0,
    vmax7 = 0
  ),
 
  # half-saturation constant for microbial d3ecomnp of pool i      
  km = list(   
    km1 = 0.01,       
    km2 = 0.01,       
    km3 = 0.01,
    km4 = 0,
    km5 = 0,
    km6 = 0,
    km7 = 0
  ),       

  # turnover rate for linear decomposition
  k = list(   
    k1 = 1.0/365,       
    k2 = 0.00005/365,   
    k3 = 1.0/365,
    k4 = .15*365,
    k5 = 3.652968e-5,
    k6 = 3.652968e-5,
    k7 = 3.652968e-5
  ),

  # maximum size for pool i 
  poolmax = list(       
    poolmax1 = 2,       # POM value
    poolmax2 = 2,       # microbial biomass max value
    poolmax3 = 26.725,   # max maom capacity (calculated using Hassink formula assuming 15% clay)
    poolmax4 = 1.7,
    poolmax5 = 26.725,
    poolmax6 = 26.725,
    poolmax7 = 26.725
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

soil_decomp_object$.test <- function(., verbose=F, metdf=F, litter=.001369863, ntimes=100 ) {

  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test() # if only used in test functions should begin with a .

  if(metdf) {
    .$dataf       <- list()
    if(length(litter)==1) litter <- rep(litter, ntimes )   
    .$dataf$metdf <- matrix(litter, nrow=1 )
    rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
    .$dataf$lm    <- length(.$dataf$metdf[1,])
    .$dataf$mout  <- .$output()
    .$run_met()
  } else {
    .$env$litter  <- litter
    .$run()
  }
}

soil_decomp_object$.test_coprpse <- function(., verbose=F, metdf=F, litter=.001369863, ntimes=100 ) {
  
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test() # if only used in test functions should begin with a .
  
  if(metdf) {
    .$dataf       <- list()
    if(length(litter)==1) litter <- rep(litter, ntimes )   
    .$dataf$metdf <- matrix(litter, nrow=1 )
    rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
    .$dataf$lm    <- length(.$dataf$metdf[1,])
    .$dataf$mout  <- .$output()
    .$run_met()
  } else {
    .$env$litter  <- litter
    .$run()
  }
}

soil_decomp_object$.test_changepool <- function(., verbose=F, metdf=F, litter=.00016, ntimes=100, n_pool=3) {
  
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$pars$n_pools = n_pool #.$build_pool_structure 
  .$fnames$transfer$t1_to_2 <- 'f_transfer_cue'
  .$fnames$transfer$t1_to_3 <- 'f_transfer_zero'
  .$fnames$transfer$t2_to_1 <- 'f_transfer_zero'
  .$fnames$transfer$t2_to_3 <- 'f_transfer_all'
  .$fnames$transfer$t3_to_1 <- 'f_transfer_zero'
  .$fnames$transfer$t3_to_2 <- 'f_transfer_cue'
  .$fnames$decomp$d1 <- 'f_decomp_MM_microbe'
  .$fnames$decomp$d2 <- 'f_decomp_lin'
  .$fnames$decomp$d3 <- 'f_decomp_MM_microbe'
  .$configure_test() # if only used in test functions should begin with a .
  
  if(metdf) {
    .$dataf       <- list()
    if(length(litter)==1) litter <- rep(litter, ntimes )   
    .$dataf$metdf <- matrix(litter, nrow=1 )
    rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
    .$dataf$lm    <- length(.$dataf$metdf[1,])
    .$dataf$mout  <- .$output()
    .$run_met()
  } else {
    .$env$litter  <- litter
    .$run()
  }
}


soil_decomp_object$.test_3pool <- function(., verbose=F, metdf=F, litter=0.00384, ntimes=365, time=T ) {

  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  soil_decomp_object$pars$n_pools = 3 

  .$configure_test() # if only used in test functions should begin with a .

  # initialise boundary data 
  .$dataf       <- list()
  if(length(litter)==1) litter <- rep(litter, ntimes )   
  .$dataf$metdf <- matrix(litter, nrow=1 )
  rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
  .$dataf$lm    <- length(.$dataf$metdf[1,])
  .$dataf$mout  <- .$output()

  ### Run models
  olist <- list()
  # run default no saturation or DD model
  print('')
  print('')
  print('Config: 1')
  print('')
  olist$noSaturation  <- .$run_met()

  # saturating MAOM
  print('')
  print('')
  print('Config: 2')
  print('')
  .$fnames$transfer$t2_to_3 <- 'f_transfer_cue_sat'
  .$configure_test() 
  olist$MaomMax       <- .$run_met()

  # denisty dependent microbial turnover 
  print('')
  print('')
  print('Config: 3')
  print('')
  .$fnames$transfer$t2_to_3 <- 'f_transfer_cue'
  .$fnames$decomp$d2        <- 'f_decomp_dd_georgiou'
  .$configure_test() 
  olist$DDturnover    <- .$run_met()

  # denisty dependent microbial cue 
  print('')
  print('')
  print('Config: 4')
  print('')
  .$fnames$transfer$t1_to_2 <- 'f_transfer_cue_sat'
  .$fnames$transfer$t3_to_2 <- 'f_transfer_cue_sat'
  .$fnames$decomp$d2        <- 'f_decomp_lin'
  .$configure_test() 
  olist$DDcue         <- .$run_met()

  # denisty dependent microbial turnover and cue 
  print('')
  print('')
  print('Config: 5')
  print('')
  .$fnames$decomp$d2        <- 'f_decomp_dd_georgiou'
  .$configure_test() 
  olist$DDturnover.DDcue <- .$run_met()

  # denisty dependent microbial turnover and cue and MAOM saturation 
  print('')
  print('')
  print('Config: 6')
  print('')
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
