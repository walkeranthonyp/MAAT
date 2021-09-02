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

  sys                = 'f_sys_npools', # does this need to change when working with steadystate?
  solver             = 'plsoda',
  solver_func        = 'f_solver_func_mimics',
  input              = 'f_input_mimics',
  DotO               = 'f_DotO',
  transfermatrix     = 'f_transfermatrix',
  steadystate        = 'f_steadystate_npools',#'f_steadystate_npools',
  solver_steadystate = 'pstode',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_rmm_wieder', #reverse michaelis menten (function accounts for two potential catalyst pools)
    d2 = 'f_decomp_rmm_wieder', 
    d3 = 'f_decomp_lin', #linear turnover of r-selected microbial biomass 
    d4 = 'f_decomp_lin', #linear turnover of k-selected microbial biomass
    d5 = 'f_zero',  #no decomp, just desorption
    d6 = 'f_decomp_rmm_wieder',
    d7 = 'f_decomp_rmm_wieder'
  ),
  
  desorp = list(
    ds5 = 'f_decomp_lin',
    ds6 = 'f_zero',
    ds7 = 'f_zero' 
  ),
  
  sorp = list(
    s1 = 'f_decomp_lin',
    s2 = 'f_decomp_lin',
    s3 = 'f_decomp_lin'    
  ),
  
  aggform = list(
    NA
  ),
  
  docuptake = 'f_docuptake_abramoff',
  
  scor = 'f_zero',
  wcor = 'f_wcor_abramoff',
  # tcor = 'f_tcor_wieder', #use this structure for millennial or models with only one scalar
  tcor = list(                  #this structure for corpse
    t1 = 'f_tcor_arrhenius',
    t2 = 'f_tcor_arrhenius',
    t3 = 'f_tcor_arrhenius'
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
  )
)


# environment
####################################
soil_decomp_object$env <- list(
  litter = 172/365,
  temp   = 20,
  vwc    = .24,
  porosity = 0.5,
  clay = .4, #MIMICS takes proportion, CORPSE and MILLENNIAL takes percentage
  lignin = 16.6,
  N = 1.37,
  anpp = 500,
  depth = 20,
  pH = 7,
  BD = 1350  #Bulk density 
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
    input_coef1 = .66666,
    input_coef2 = 0,
    input_coef3 = 0,
    input_coef4 = .33334,
    input_coef5 = 0
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
    cue1 = .5,       
    cue2 = .25,       
    cue3 = 0,
    cue4 = 0,
    cue5 = 0, 
    cue6 = 0,
    cue7 = .5
  ),  
  
  cue2 = list(
    cue1 = .7,       
    cue2 = .35,       
    cue3 = 0,
    cue4 = 0,
    cue5 = 0, 
    cue6 = 0,
    cue7 = .7
  ),  

  # max turnover rate per unit microbial biomass for pool i
  # commented out old list (updating to allow for multiple mic groups or catalysts)
  vmax = list(
    vmax1 = 10, #Vpd or Vpl in Millennial, decomp of POM toward DOC
    vmax2 = 2,
    vmax3 = 0,
    vmax4 = 0,
    vmax5 = 0,
    vmax6 = 0,
    vmax7 = 10
  ),
  
  vmax2 = list( #vmax for second catalyst
    vmax1 = 3, #Vpd or Vpl in Millennial, decomp of POM toward DOC
    vmax2 = 3,
    vmax3 = 0,
    vmax4 = 0,
    vmax5 = 0,
    vmax6 = 0,
    vmax7 = 2
  ),
  
  # idea for how to handle one pool being decomposed by multiple catalysts (e.g. multiple microbial pools or multiple enzyme pools)
  # vmax = list(
  #   #vmax of pool 1
  #   vmax1 = list(
  #     #specific to first catalyst (mic or enz pool)
  #     cat1 = 0,
  #     #specific to second catalyst (mic or enz pool)
  #     cat2 = 0),
  #   vmax2 = list(cat1 = 0, cat2 = 0),
  #   vmax3 = list(cat1 = 0, cat2 = 0),
  #   vmax4 = list(cat1 = 0, cat2 = 0),
  #   vmax5 = list(cat1 = 0, cat2 = 0)
  # ),
 
  # michaelis-menten half-saturation constant for microbial decomnp of pool i      
  km = list(   
    km1 = 8, #Kpd in MILLENNIAL double michaelis-menten equation     
    km2 = 2,     
    km3 = 0,
    km4 = 0,
    km5 = 0,
    km6 = 0,
    km7 = 4
  ),
  
  km2 = list(   #km for second catalyst
    km1 = 2, #Kpd in MILLENNIAL double michaelis-menten equation     
    km2 = 4,     
    km3 = 0,
    km4 = 0,
    km5 = 0,
    km6 = 0,
    km7 = 6
  ),
  
  # reverse michaelis-menten half-saturation constant for microbial decomnp of pool i
  # This is only neccesary for models that use km in reverse and forward mm decomp for the same pool (i.e. MILLENNIAL)
  rkm = list(   
    rkm1 = 12, #Kpe in MILLENNIAL double michaelis-menten equation (the reverse part i.e. the microbial limited part) 
    rkm2 = 0,     
    rkm3 = 0,
    rkm4 = 0,
    rkm5 = 0
  ),   

  # turnover rate for linear decomposition
  k = list(    
    k1 = 0,       
    k2 = 0,   
    k3 = 0, #calculated by MIMICS solver function 
    k4 = 0, #calculated by MIMICS solver function 
    k5 = 0, 
    k6 = 0,
    k7 = 0  #calculated by MIMICS solver function 
  ),

  # maximum size for pool i 
  poolmax = list(       
    poolmax1 = 0,       # POM value
    poolmax2 = 0,       # microbial biomass max value
    poolmax3 = 0,   # max maom capacity (calculated using Hassink formula assuming 15% clay)
    poolmax4 = 0,
    poolmax5 = 500
  ),
  
  #mimics-specific parameters
  mimics = list(
    fmet_p1 = .5,
    fmet_p2 = .85,
    fmet_p3 = .013,
    tau_mod1_p1 = 100,
    tau_mod1_p2 = .6,
    tau_mod1_p3 = 1.3,
    tau_r_p1 = .00052,
    tau_r_p2 = .3,
    tau_mod2 = 2,
    tau_k_p1 = .00024,
    tau_k_p2 = .1,
    desorb_p1 = .00002,
    desorb_p2 = -4.5,
    fSOMp_r_p1 = .15,
    fSOMp_r_p2 = 1.3,
    fSOMc_r_p1 = .1,
    fSOMc_r_p2 = -3,
    fSOMc_r_p3 = 1,
    fSOMp_k_p1 = .1,
    fSOMp_k_p2 = .8,
    fSOMc_k_p1 = .3,
    fSOMc_k_p2 = -3,
    fSOMc_k_p3 = 1,
    V_slope = .063,
    V_int = 5.47,
    aV = .000000125,
    pscalar_p1 = 3,
    pscalar_p2 = -2,
    ko_r = 6,
    ko_k = 6,
    K_slope = .02,
    K_int = 3.19,
    aK = .15625,
    fi_LITm = .005,
    fi_LITs = .005
  ),
  
  #millennial-specific parameters
  millennial = list(
    cuet = -0.012, #slope relating assimilation efficiency to temperature
    Taeref = 15,   #ref temp for CUE temp-dependence equation
    Vpa = 0.002,   #max aggregation rate of POM
    Kpa = 50,      #half-saturation constant for aggregation of POM
    c1 = 0.297,    #parameter relating clay to Qmax (i.e. MAOM poolmax)
    c2 = 3.355,    #parameter relating clay to Qmax (i.e. MAOM poolmax)
    Vma = 0.07,    #max aggregation rate of MAOM
    Kma = 200,     #half-saturation constant for aggregation of MAOM
    Vdm = 0.35,    #max DOC turnover rate (for microbial uptake, not leaching or sorption)
    Kdb = 7.2,     #half-saturation constant for microbial uptake of doc
    kmm = 0.025,   #rate constant for microbial turnover (and sorption in published eq version)
    pa = 0.333
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

soil_decomp_object$.test <- function(., verbose=F, metdf=F, litter=172/365, ntimes=100 ) {

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

soil_decomp_object$.test.ss <- function(., verbose=F, metdf=F, litter=172/365, ntimes=100 ) {
  
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

soil_decomp_object$.test_corpse <- function(., verbose=F, metdf=F, litter=.001369863, ntimes=100 ) {
  
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

# soil_decomp_object$.test_changepool <- function(., verbose=F, metdf=F, litter=.00016, ntimes=100, n_pool=3) {
#   
#   if(verbose) str(.)
#   .$build(switches=c(F,verbose,F))
#   .$pars$n_pools = n_pool #.$build_pool_structure 
#   .$fnames$transfer$t1_to_2 <- 'f_transfer_cue'
#   .$fnames$transfer$t1_to_3 <- 'f_transfer_zero'
#   .$fnames$transfer$t2_to_1 <- 'f_transfer_zero'
#   .$fnames$transfer$t2_to_3 <- 'f_transfer_all'
#   .$fnames$transfer$t3_to_1 <- 'f_transfer_zero'
#   .$fnames$transfer$t3_to_2 <- 'f_transfer_cue'
#   .$fnames$decomp$d1 <- 'f_decomp_MM_microbe'
#   .$fnames$decomp$d2 <- 'f_decomp_lin'
#   .$fnames$decomp$d3 <- 'f_decomp_MM_microbe'
#   .$configure_test() # if only used in test functions should begin with a .
#   
#   if(metdf) {
#     .$dataf       <- list()
#     if(length(litter)==1) litter <- rep(litter, ntimes )   
#     .$dataf$metdf <- matrix(litter, nrow=1 )
#     rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
#     .$dataf$lm    <- length(.$dataf$metdf[1,])
#     .$dataf$mout  <- .$output()
#     .$run_met()
#   } else {
#     .$env$litter  <- litter
#     .$run()
#   }
# }


# soil_decomp_object$.test_3pool <- function(., verbose=F, metdf=F, litter=0.00384, ntimes=365, time=T ) {
# 
#   if(verbose) str(.)
#   .$build(switches=c(F,verbose,F))
#   soil_decomp_object$pars$n_pools = 3 
# 
#   .$configure_test() # if only used in test functions should begin with a .
# 
#   # initialise boundary data 
#   .$dataf       <- list()
#   if(length(litter)==1) litter <- rep(litter, ntimes )   
#   .$dataf$metdf <- matrix(litter, nrow=1 )
#   rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
#   .$dataf$lm    <- length(.$dataf$metdf[1,])
#   .$dataf$mout  <- .$output()
# 
#   ### Run models
#   olist <- list()
#   # run default no saturation or DD model
#   print('')
#   print('')
#   print('Config: 1')
#   print('')
#   olist$noSaturation  <- .$run_met()
# 
#   # saturating MAOM
#   print('')
#   print('')
#   print('Config: 2')
#   print('')
#   .$fnames$transfer$t2_to_3 <- 'f_transfer_cue_sat'
#   .$configure_test() 
#   olist$MaomMax       <- .$run_met()
# 
#   # denisty dependent microbial turnover 
#   print('')
#   print('')
#   print('Config: 3')
#   print('')
#   .$fnames$transfer$t2_to_3 <- 'f_transfer_cue'
#   .$fnames$decomp$d2        <- 'f_decomp_dd_georgiou'
#   .$configure_test() 
#   olist$DDturnover    <- .$run_met()
# 
#   # denisty dependent microbial cue 
#   print('')
#   print('')
#   print('Config: 4')
#   print('')
#   .$fnames$transfer$t1_to_2 <- 'f_transfer_cue_sat'
#   .$fnames$transfer$t3_to_2 <- 'f_transfer_cue_sat'
#   .$fnames$decomp$d2        <- 'f_decomp_lin'
#   .$configure_test() 
#   olist$DDcue         <- .$run_met()
# 
#   # denisty dependent microbial turnover and cue 
#   print('')
#   print('')
#   print('Config: 5')
#   print('')
#   .$fnames$decomp$d2        <- 'f_decomp_dd_georgiou'
#   .$configure_test() 
#   olist$DDturnover.DDcue <- .$run_met()
# 
#   # denisty dependent microbial turnover and cue and MAOM saturation 
#   print('')
#   print('')
#   print('Config: 6')
#   print('')
#   .$fnames$transfer$t2_to_3 <- 'f_transfer_cue_sat'
#   .$configure_test() 
#   olist$DDturnover.DDcue.MaomMax <- .$run_met()
# 
# 
#   # plotting functions
#   thp_plot_time <- function(mod) {
#     ylab <- expression('Pool C mass ['*gC*' '*m^-2*']')
#     matplot(1:dim(mod)[1], mod[,1:3], type='l', ylab=ylab, xlab='Days', lty=1,
#             ylim=c(0,max(mod)*1.2), col=1:3,main=deparse(substitute(mod)) )
#     legend('topleft', c('POM','MB','MAOM'), lty=1, col=1:3, bty='n')
#   }
#  
#   thp_plot_MBC <- function(mod) {
#     matplot(mod[,2], mod[,c(1,3)], type='l', ylab=ylab, xlab='MB C mass', lty=1,
#             xlim=c(0,max(mod[,2])),
#             ylim=c(0,max(mod)*1.2), col=1:2, main=deparse(substitute(mod)) )
#     legend('topleft', c('POM','MAOM'), lty=1, col=1:2, bty='n')
#   }
#   
#   par(mfrow = c(2,3))
#   if(time) { 
#     ## plotting versus time
#     thp_plot_time(olist$noSaturation)
#     thp_plot_time(olist$MaomMax)
#     thp_plot_time(olist$DDturnover)
#     thp_plot_time(olist$DDcue)
#     thp_plot_time(olist$DDturnover.DDcue)
#     thp_plot_time(olist$DDturnover.DDcue.MaomMax)
#   } else {  
#     # plotting pools vs MBC
#     thp_plot_MBC(olist$noSaturation)
#     thp_plot_MBC(olist$MaomMax)
#     thp_plot_MBC(olist$DDturnover)
#     thp_plot_MBC(olist$DDcue)
#     thp_plot_MBC(olist$DDturnover.DDcue)
#     thp_plot_MBC(olist$DDturnover.DDcue.MaomMax)
#   }
# 
#   olist
# }
#   
# 

### END ###
