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
source('../../functions/packagemod_functions_deSolve.R')
source('../../functions/packagemod_functions_rootSolve.R')


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
soil_decomp_object$pool_pars <- c('cstate0', 'cue', 'vmax', 'km', 'k', 'poolmax' )



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
  
  # call steady state system model
  # APW Matt: this is the additional function call to initialise the model at steady state
  #           if you want to turn this off just set the fnames$steadystate value to f_steadystate_null (a dummy function that just returns NULL)
  .$fns$steadystate()
}



# assign object variables 
###########################################################################

# function names
####################################
soil_decomp_object$fnames <- list(

  sys                = 'f_sys_npools',
  solver             = 'plsoda',
  solver_func        = 'f_solver_func',
  input              = 'f_input',
  DotO               = 'f_DotO',
  transfermatrix     = 'f_transfermatrix',
  steadystate        = 'f_steadystate_npools',
  solver_steadystate = 'pstode',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_MM_microbe',
    d2 = 'f_decomp_lin',
    d3 = 'f_decomp_MM_microbe'
  ),
  
  # transfer list
  transfer = list(
    t1_to_2 = 'f_transfer_cue',
    t2_to_3 = 'f_transfer_cue',
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
  cpools = matrix(1:3, ncol=1 )
)


# state parameters (i.e. calculated parameters)
####################################
soil_decomp_object$state_pars <- list(
  solver_out             = matrix(1),
  solver_steadystate_out = matrix(1)
)


# parameters
####################################
soil_decomp_object$pars <- list(

  n_pools = 3,          # number of pools in model  
  beta    = 1.5,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  silt    = 0.2,        # soil silt content (proportion)
  clay    = 0.2,        # soil clay content (proportion)
 
  # initial pool mass for each pool
  cstate0 = list(
    cstate01 = 0.6,
    cstate02 = 0.1,
    cstate03 = 0.1
  ),

  # Carbon use or transfer efficiency from pool i to any another 
  # - if this varies by the 'to' pool we need another function / parameters
  cue = list(
    cue1 = 0.47,       # currently, all cue values are the same (at the value in MEND), might need a different cue_max for density dependent function... 
    cue2 = 0.566,      # humification constant
    cue3 = 0.47   
  ),  

  # max turnover rate per unit microbial biomass for pool i 
  vmax = list(   
    vmax1 = 0.2346,     
    vmax2 = 0.2346,     
    vmax3 = 0.0777   
  ),
 
  # half-saturation constant for microbial d3ecomnp of pool i      
  km = list(   
    km1 = 101,       
    km2 = 101,       
    km3 = 250
  ),       

  # turnover rate for linear decomposition
  k = list(   
    k1 = 0.007,       
    k2 = 0.00672,    # microbial turnover constant
    k3 = 0.006
  ),

  # maximum size for pool i 
  poolmax = list(       
    poolmax1 = 2,       # POM value
    poolmax2 = 2,       # microbial biomass max value
    poolmax3 = 26.725   # max maom capacity (calculated using Hassink formula assuming 15% clay)
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

f_output_soil_decomp_eval <- f_output_eval 


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
    .$dataf$metdf <- matrix(litter, nrow=1 )
    rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
    .$dataf$lm    <- dim(.$dataf$met)[2]
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
  .$dataf$metdf <- matrix(litter, nrow=1 )
  rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
  .$dataf$lm    <- dim(.$dataf$met)[2]
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
  .$fnames$decomp$d2        <- 'f_decomp_dd'
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
  .$fnames$decomp$d2        <- 'f_decomp_dd'
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


soil_decomp_object$.test_var_kinetics_yearly <- function(., verbose=F, litter=1.4, ntimes=200, kinetics = 'mm', 
                                                vmax1 = 88, vmax3 = 171, km1 = 144, km3 = 936, k1 = 1, k3 = .01) {
  
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test() # if only used in test functions should begin with a .
  
  # initialise boundary data 
  .$dataf       <- list()
  if(length(litter)==1) litter <- rep(litter, ntimes )   
  .$dataf$metdf <- matrix(litter, nrow=1 )
  rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
  .$dataf$lm    <- length(.$dataf$metdf[1,])
  .$dataf$mout  <- .$output()
  
  
  ### Run model
  # change parameters to default parms
  if(kinetics == 'rmm'){
    .$fnames$decomp$d1        <- 'f_decomp_RMM_microbe'
    .$fnames$decomp$d3        <- 'f_decomp_RMM_microbe'
  } else if (kinetics == 'lin') {
    .$fnames$decomp$d1        <- 'f_decomp_lin'
    .$fnames$decomp$d3        <- 'f_decomp_lin'
    .$pars$k$k1 = k1
    .$pars$k$k3 = k3
  } else if (kinetics == 'mm'){
    .$fnames$decomp$d1        <- 'f_decomp_MM_microbe'
    .$fnames$decomp$d3        <- 'f_decomp_MM_microbe'
  }
  

    .$pars$cue$cue1 = .47  #CUE from MEND (Wang et al. 2013)
    .$pars$cue$cue2 = 1    #NA; assuming mbc turnover is entirely transferred to MAOM pool
    .$pars$cue$cue3 = .47   #CUE from MEND (Wang et al. 2013)
  
  

    .$pars$vmax$vmax1 = vmax1  #MIMICS average of two microbial groups assuming 15 degC (Wieder et al. 2014)   
    .$pars$vmax$vmax3 = vmax3  #MIMICS average of two microbial groups assuming 15 degC (Wieder et al. 2014)  
  
  

    .$pars$km$km1 = km1  #MIMICS average of two microbial groups for sturctural litter and biochemically protected SOC assuming 15 degC and 15%clay (Wieder et al. 2014)        
    .$pars$km$km3 = km3   #MIMICS average of two microbial groups assuming 15 degC and lignin:N = 10 (Wieder et al. 2014)
  
  
  .$pars$k$k2 = 2.5    #Li et al. 
  
  
  .$pars$poolmax = list(       
    pmax1 = 2,      #NA
    pmax2 = 1.5,    #about 97.5% quantile of mbc synthesis data  
    pmax3 = 27      #Calculated from Hassink and Whitmore 1997 assuming 20% clay
  )
  
  .$pars$beta = 1.5
  
  
  .$configure_test()
  out  <- .$run_met()
  
  par(mfrow=c(1,1))
  ylab <- expression('Pool C mass ['*gC*' '*m^-2*']')
  matplot(1:dim(out)[1], out[,1:3], type='l', ylab=ylab, xlab='Years', lty=1,
          ylim=c(0,max(out)*1.2), col=1:3 )
  legend('topleft', c('POM','MB','MAOM'), lty=1, col=1:3, bty='n')
  
  print(tail(out, n=1))
}

### END ###
