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
  solver_func    = 'f_solver_func',
  input          = 'f_input',
  DotO           = 'f_DotO',
  transfermatrix = 'f_transfermatrix',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_MM_enzpom',
    d2 = 'f_decomp_mbc_mend',
    d3 = 'f_decomp_MM_enzmaom',
    d4 = 'f_decomp_dd_mend',
    d5 = 'f_decomp_doc_mend',
    d6 = 'f_decomp_lin',
    d7 = 'f_decomp_lin'
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
  litter = 3.2,
  temp   = 10,
  swc    = 0.3
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
)


# parameters
####################################
soil_decomp_object$pars <- list(

  n_pools = 7,          # number of pools in model  
  beta    = 1.5,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  silt    = 0.2,        # soil silt content (proportion)
  clay    = 0.2,        # soil clay content (proportion)
  mr      = 0.00028,     # specific maintenance factor (MEND)
  pep     = 0.01,        #Fraction of mr allocated for production of EP
  pem     = 0.01,        #Fraction of mr allocated for production of EM
  Kads    = 0.006,       #Specific adsorption rate (could make this k for the DOC pool)
 
  #input coefficients (allocaties proportions of inputs into different pools)
  input_coefs = list(
    input_coef1 = 1,
    input_coef2 = 0,
    input_coef3 = 0,
    input_coef4 = 0,
    input_coef5 = 0,
    input_coef6 = 0,
    input_coef7 = 0
  ),
  
  # initial pool mass for each pool
  cstate0 = list(
    cstate01 = 10,         #Initial POM pool size
    cstate02 = 2,          #Initial MBC pool size
    cstate03 = 5,          #Initial MAOM pool size
    cstate04 = 0.1,        #Initial Q pool size
    cstate05 = 1,          #Initial DOC pool size
    cstate06 = 0.00001,    #Initial EP pool size
    cstate07 = 0.00001     #Initial EM pool size
  ),

  # Carbon use or transfer efficiency from pool i to any another 
  # - if this varies by the 'to' pool we need another function / parameters
  # MC: I've created a function cue_remainder that can allocate the remainder to another pool rather than CO2
  cue = list(
    cue1 = 0.5,       #Fd from MEND
    cue2 = 0.5,       #Gd from MEND
    cue3 = 0.47,
    cue4 = 0.47,
    cue5 = 0.47,      #Ec from MEND
    cue6 = 0.47,
    cue7 = 0.47
  ),  

  # max turnover rate per unit microbial biomass for pool i 
  vmax = list(   
    vmax1 = 2.5,     
    vmax2 = 1,     
    vmax3 = 1,
    vmax4 = 1,
    vmax5 = 0.26,
    vmax6 = 1,
    vmax7 = 1
  ),
 
  # half-saturation constant for microbial d3ecomnp of pool i      
  km = list(   
    km1 = 50,       
    km2 = 1,       
    km3 = 250,
    km4 = 1,
    km5 = .26,
    km6 = 1,
    km7 = 1
  ),       

  # turnover rate for linear decomposition
  k = list(   
    k1 = 0.007,       
    k2 = 0.00672,    # microbial turnover constant
    k3 = 0.001,
    k4 = 0.001,
    k5 = 0.006,
    k6 = 0.001,
    k7 = 0.001
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

soil_decomp_object$.test <- function(., verbose=F, metdf=F, litter=.00016, ntimes=100 ) {

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
