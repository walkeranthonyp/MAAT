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
source('soil_decomp_solver_functions.R')
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

  sys                = 'f_sys_npools', 
  solver             = 'plsoda',
  solver_func        = 'f_solver_func_millennialV2',
  input              = 'f_input',
  DotO               = 'f_DotO',
  transfermatrix     = 'f_transfermatrix',
  steadystate        ='f_steadystate_npools', # options: 'f_steadystate_null','f_steadystate_npools',
  solver_steadystate = 'pstode',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_rmm',
    d2 = 'f_decomp_dd_georgiou',  
    d3 = 'f_zero', 
    d4 = 'f_decomp_lin', 
    d5 = 'f_decomp_lin',
    d6 = NA,
    d7 = NA
  ),
  
  desorp = list(
    ds1 = NA,
    ds2 = NA,
    ds3 = 'f_desorp_millennialv2',
    ds4 = NA,
    ds5 = NA,
    ds6 = NA,
    ds7 = NA
  ),
  
  sorp = list(
    s1 = NA,
    s2 = NA,
    s3 = NA,
    s4 = 'f_sorp_sat',
    s5 = NA,
    s6 = NA,
    s7 = NA
  ),
  
  aggform = list( 
    a1 = 'f_decomp_lin', #Aggregate formation from POM
    a3 = 'f_decomp_lin'  #Aggregate formation from MAOM
  ),
  
  docuptake = 'f_decomp_mm', #MM uptake of DOC
  
  growthresp = NA,
  
  maintresp = NA,
  
  scor = NA,
  wcor = 'f_wcor_ghezzehei_diffusion',
  wcor2 = 'f_wcor_ghezzehei_biological',

  #tcor = 'f_tcor_wieder', #use this structure for millennial or models with only one scalar
  tcor = list(                  #this structure for corpse
    t1 = 'f_tcor_arrhenius_millennialv2',
    t2 = NULL,
    t3 = NULL,
    t4 = 'f_tcor_arrhenius_millennialv2',
    t5 = NULL
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
  litter = 172.8978/365, #forc_npp in MILLENNIALv2
  temp   = 11.21961, #default for MILLENNIALv2 (sum across mean year)
  vwc    = .2422044, #default for MILLENNIALv2 (average across mean year)
  porosity = 0.6, #default for MILLENNIALv2
  clay = .0, 
  lignin = 0,
  N = 0,
  anpp = 0,
  depth = 0,
  pH = 7, #default for MILLENNIALv2
  BD = 1000,  #default for MILLENNIALv2 #Bulk density 
  matpot = 15, #default for MILLENNIALv2
  lambda = 2.1000e-04, #default for MILLENNIALv2
  kamin = .2, #default for MILLENNIALv2
  claysilt = 80  #default for MILLENNIALv2 #MILLENNIAL uses clay+silt% to calculate Qmax
)


# state
####################################
soil_decomp_object$state <- list(
  cpools = matrix(1:3, ncol=1 ),
  respiration = vector("double",1)
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

  n_pools = 7,          # number of pools in model  #need to change and re-create XMLs when this changes for wrapper runs
  beta    = 2,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  silt    = NA,      
  clay    = NA,        
  clayref = NA,
  mr      = NA,     
  pep     = NA,        
  pem     = NA,        
  Kads    = NA,       
  qslope_mayes = NA,
  reftemp = NA,
  R = 8.31446,#72,        #used in MILLENNIALv2
  minmic = NA,
  q10 = NA,
  
  ea = list(
    ea1 = 6.3909e+04, #MILLENNIALv2 #activiation energy for temp sensitivity of POM decay
    ea2 = NA,
    ea3 = NA,
    ea4 = 5.7865e+04, #MILLENNIALv2 #activiation energy for temp sensitivity of DOC uptake
    ea5 = NA
  ),
 
  #input coefficients (allocaties proportions of inputs into different pools)
  input_coefs = list(
    input_coef1 = .66,
    input_coef2 = 0,
    input_coef3 = 0,
    input_coef4 = .34,
    input_coef5 = 0,
    input_coef6 = NA,
    input_coef7 = NA
  ),
  
  # initial pool mass for each pool
  cstate0 = list( 
    cstate01 = 1, 
    cstate02 = 1,          
    cstate03 = 1,      
    cstate04 = 1,       
    cstate05 = 1,
    cstate06 = 1,
    cstate07 = 1
  ),

  # Carbon use or transfer efficiency from pool i to any another 
  # - if this varies by the 'to' pool we need another function / parameters
  # MC: I've created a function cue_remainder that can allocate the remainder to another pool rather than CO2
  cue = list(
    cue1 = NA,       
    cue2 = NA,       
    cue3 = NA,
    cue4 = .19,
    cue5 = NA,
    cue6 = NA,
    cue7 = NA
  ),  
  
  cue2 = list( 
    cue1 = NA,       
    cue2 = NA,       
    cue3 = NA,
    cue4 = NA,
    cue5 = NA,
    cue6 = NA,
    cue7 = NA
  ),  

  # max turnover rate per unit microbial biomass for pool i
  # commented out old list (updating to allow for multiple mic groups or catalysts)
  vmax = list(
    vmax1 = 1.8000e+12,#alpha_pl; pre-exponential constant for temp sensitivity of POM decay
    vmax2 = NA,   #Vm
    vmax3 = NA,
    vmax4 = 2.3000e+12, #alpha_lb #pre-exponential constant for temp sensitivity of DOC uptake
    vmax5 = NA,
    vmax6 = NA,
    vmax7 = NA
  ),
  
  vmax2 = list(
    vmax1 = NA, 
    vmax2 = NA,
    vmax3 = NA,
    vmax4 = NA,
    vmax5 = NA,
    vmax6 = NA,
    vmax7 = NA
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
    km1 = 6443, #MillennialV2 #half sat const for POM breakdown 
    km2 = NA,
    km3 = NA,
    km4 = 774.6, #DOC uptake half sat constant,
    km5 = NA,
    km6 = NA,
    km7 = NA
  ),
  
  km2 = list( 
    km1 = NA, 
    km2 = NA,     
    km3 = NA,
    km4 = NA,
    km5 = NA,
    km6 = NA,
    km7 = NA
  ),
  
  # reverse michaelis-menten half-saturation constant for microbial decomnp of pool i
  # This is only neccesary for models that use km in reverse and forward mm decomp for the same pool (i.e. MILLENNIALv1)
  rkm = list(  
    rkm1 = NA, 
    rkm2 = NA,     
    rkm3 = NA,
    rkm4 = NA,
    rkm5 = NA,
    rkm6 = NA,
    rkm7 = NA
  ),   

  # turnover rate for linear decomposition
  k = list(    
    k1 = .018, #aggregate formation from POM       
    k2 = 4.5000e-03, #microbial turnover rate   
    k3 = 4.8000e-03, #aggregate formation from maom 
    k4 = .0015, #leaching rate (kd)
    k5 = .02, #aggregrate breakdown rate (kb)
    k6 = NA,
    k7 = NA
  ),

  # maximum size for pool i 
  poolmax = list(       
    poolmax1 = NA,      
    poolmax2 = NA,       
    poolmax3 = NA, #this is calculated in solver function 
    poolmax4 = NA,
    poolmax5 = NA,
    poolmax6 = NA,
    poolmax7 = NA
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
  ),
  
  #mend-specific parameters
  mend = list(
    Mr =  .00028,
    Pep = .01,
    Pem = .01,
    Gd =  .5,
    Fd =  .5
  ),
  
  millennialV2 = list(
    param_pc = .86, ##slope of mineral C - clay relationship from Georgiou et al. in review
    kld = 1, #desorption coefficient
    sorp_p1 = .12, #sorption affinity parameter
    sorp_p2 = .216, #sorption affinity parameter
    cue_t = 0.012, #slope of CUE temp sensitivity
    Taeref = 15,   #ref temp for CUE temp-dependence equation
    pa = 0.33,      #proportion agg breakdown into pom
    param_pb = .5 #fraction of mb turnover to maom vs doc
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

soil_decomp_object$.test <- function(., verbose=F, metdf=F, litter=.00016, ntimes=100 ) {

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



soil_decomp_object$.test.mimics.ss <- function(., verbose=F, metdf=F, litter=172/365, ntimes=100 ) {
  
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test() # if only used in test functions should begin with a .
  
    .$fnames <- list(
        sys            = 'f_sys_npools',
        solver         = 'plsoda',
        solver_func    = 'f_solver_func_mimics',
        input          = 'f_input_mimics',
        steadystate        = 'f_steadystate_npools',
        solver_steadystate = 'pstode',
        decomp = list(
          d1 = 'f_decomp_rmm_wieder', #reverse michaelis menten (function accounts for two potential catalyst pools)
          d2 = 'f_decomp_rmm_wieder', 
          d3 = 'f_decomp_lin', #linear turnover of r-selected microbial biomass 
          d4 = 'f_decomp_lin', #linear turnover of k-selected microbial biomass
          d6 = 'f_decomp_rmm_wieder',
          d7 = 'f_decomp_rmm_wieder'
        ),
        desorp = list(
          ds5 = 'f_decomp_lin'
        ),
        tcor = 'f_tcor_wieder'
      )
    .$pars <- list(
      n_pools = 7,          # number of pools in model  
      
      # initial pool mass for each pool
      # values are equilibrium calculated with stode function in DeSolve script
      cstate0 = list(
        cstate01 = 0.7, 
        cstate02 = 4,          
        cstate03 = .09,      
        cstate04 = .2,       
        cstate05 = 2.54,
        cstate06 = 2.32,
        cstate07 = 1.5
      ),
      
      # Carbon use efficiency from pool i r-microbes (cue) or k-microbes (cue2)
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
      
      # max turnover rate per unit r-(vmax) or k-(vmax) microbial biomass for pool i
      vmax = list(
        vmax1 = 10, 
        vmax2 = 2,
        vmax3 = 0,
        vmax4 = 0,
        vmax5 = 0,
        vmax6 = 0,
        vmax7 = 10
      ),
      
      vmax2 = list( #vmax for second catalyst (i.e. k-selected microbes)
        vmax1 = 3, 
        vmax2 = 3,
        vmax3 = 0,
        vmax4 = 0,
        vmax5 = 0,
        vmax6 = 0,
        vmax7 = 2
      ),
      
      # michaelis-menten half-saturation constant for microbial decomp of pool i for r-(km) and
      # k-(km2) selected microbes
      km = list(   
        km1 = 8,   
        km2 = 2,
        km3 = 0,
        km4 = 0,
        km5 = 0,
        km6 = 0,
        km7 = 4
      ),
      
      km2 = list(   #km for second catalyst (i.e. k-selected microbes)
        km1 = 2,   
        km2 = 4,
        km3 = 0,
        km4 = 0,
        km5 = 0,
        km5 = 0,
        km7 = 6
      ),
      
      # turnover rate for linear decomposition
      k = list(
        k1 = 0,
        k2 = 0,
        k3 = NULL, #calculated by MIMICS solver function 
        k4 = NULL, #calculated by MIMICS solver function 
        k5 = NULL,  #calculated by MIMICS solver function 
        k6 = 0,
        k7 = 0
      ),
      
      #mimics-specific parameters
      mimics = list(
        fmet_p1 = .5,       # These three parameters used to calculate fmet 
        fmet_p2 = .85,      # fmet partitions litter between structural and metabolic pools
        fmet_p3 = .013,     # These parameters relate lignin/N of litter to fmet
        tau_mod1_p1 = 100,  # "tau_" parameters used to calculate turnover of microbial biomass pools
        tau_mod1_p2 = .6,   #
        tau_mod1_p3 = 1.3,  #
        tau_r_p1 = .00052,  #
        tau_r_p2 = .3,      #
        tau_mod2 = 2*24,       #
        tau_k_p1 = .00024,  #
        tau_k_p2 = .1,      #
        desorb_p1 = .00002*24, # "desorb_" parmeters caluclate desorption of SOMp pool based on clay content 
        desorb_p2 = -4.5,   #
        fSOMp_r_p1 = .15,   # "fSOM..." parameters control transfers of microbial necromass to the three SOM pools
        fSOMp_r_p2 = 1.3,   # based on clay or fmet parameter
        fSOMc_r_p1 = .1,    #
        fSOMc_r_p2 = -3,    #
        fSOMc_r_p3 = 1,     #
        fSOMp_k_p1 = .1,    #
        fSOMp_k_p2 = .8,    #
        fSOMc_k_p1 = .3,    #
        fSOMc_k_p2 = -3,    #
        fSOMc_k_p3 = 1,     #
        V_slope = .063,     # Vmax temp sensitivity--slope of relationship between temp and ln(Vmax)
        V_int = 5.47,       # Vmax temp sensitivity--intercept of relationship between temp and ln(Vmax)
        aV = .000000125*24,    # Tuning coefficient for Vmax 
        pscalar_p1 = 3,     # "pscalar" determines effect of clay on Km for decomp of SOMa pool
        pscalar_p2 = -2,    #
        ko_r = 6,           # tunes Km for decomp of SOMc pool
        ko_k = 6,           #
        K_slope = .02,      # temp sensitivity of Km parameter (slope)
        K_int = 3.19,       # temp sensitivity of Km parameter (intercept)
        aK = .15625,        # temp sensitivity of Km parameter (tuning coefficient)
        fi_LITm = .005,     # fraction of metabolic litter input transferred to SOMp
        fi_LITs = .005      # fraction of structural litter input transferred to SOMc
      )
    )
    .$env <- list(
      temp   = 20,          # soil temperature
      clay = .05,            # proportion clay (i.e. .4 = 40%) #default = .4, changed to .05 to match CORPSE
      lignin = 16.6,        # Lignin concentration of litter inputs (units must be same as N (%))
      N = 1.37,             # N concentration of litter inputs (units must be same as lignin (%))
      anpp = 500*24,           # ANPP (gC / m^2 / y) 
      depth = 20           # depth (cm)
    )
    
    .$run()
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
