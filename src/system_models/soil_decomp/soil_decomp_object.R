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
#source('soil_decomp_solver_functions.R')
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
soil_decomp_object$name <- 'soil_decomp'



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
  .$state$cpools  = matrix(unlist(.$pars$cstate0)[1:.$pars$n_pools], ncol=1 )
  .$state$outflux = as.list(numeric(.$pars$n_outflux))
  names(.$state$outflux) = paste0('of',1:.$pars$n_outflux)
}



# assign object variables 
###########################################################################

# function names
####################################
soil_decomp_object$fnames <- list(

  # generic functions
  sys                    = 'f_sys_npools',                    # if f_steadystate_npools selected, solver_steadystate automatically used  
  solver                 = 'plsoda',                          # solver to use for below solver function 
  solver_func            = 'f_solver_func_soilR_outfluxes',   # solver function to solve matrix equations
  steadystate            = 'f_steadystate_npools',            # solver function to solve for steady state, options: 'f_steadystate_null','f_steadystate_npools',
  solver_steadystate     = 'pstode',                          # solver to use for steady state solution 
  input                  = 'f_input',                         # input function
  DotO                   = 'f_DotO',                          # pool decomp matrix  
  transfermatrix         = 'f_transfermatrix',                # transfer matrix
  transfer_fluxsum_prop  = 'f_transfer_fluxsum_prop',         # generic function to calculate flux proportions
  outfluxes              = 'f_outfluxes',                     # generic function to calculate each outflux
  calc_state_pars        = 'f_calc_state_pars_generic',       # generic function to calculate state paramters
  calc_state_pars_unique = NA,                                # function to calculate model-specific state parameters 
  calc_loss_fluxes       = 'f_calc_loss_fluxes',              # function to calculate loss fluxes, i.e. respiration & leaching currently
  mass_balance           = 'f_mass_balance',                  # function to calculate mass balance error
 
  # decay/decomposition/flux functions, one per pool
  # - these functions represent the total flux out of a given pool
  # - can be the sum of multiple terms in the "outflux" list
  # - see pars$decomp_outflux{1...n} lists for ouflux assignments to each pool 
  decomp = list(
    d1 = 'f_decomp_fluxsum',
    d2 = 'f_decomp_fluxsum',
    d3 = 'f_decomp_fluxsum',
    d4 = 'f_decomp_fluxsum',
    d5 = 'f_decomp_fluxsum'
  ),

  # functions for state parameters indexed by pool  
  poolmax = list(
    poolmax3 = 'f_poolmax_texture_abramoff'
  ),
 
  # output fluxes from pools functions 
  # - one function per flux, designed for when there are more than one output flux from a pool
  outflux = list(
    of1 = 'f_decomp_rmm',
    of2 = 'f_decomp_lin',
    of3 = 'f_decomp_dd_georgiou',  
    of4 = 'f_desorp_sat', 
    of5 = 'f_decomp_lin', 
    of6 = 'f_decomp_mm',
    of7 = 'f_sorp_sat',
    of8 = 'f_decomp_lin',
    of9 = 'f_decomp_lin'
  ),

  # functions that align with each outflux 
  # - temperature scalars 
  tcor = list(             
    tcor1 = 'f_tcor_arrhenius',
    tcor2 = 'f_tcor_none',
    tcor3 = 'f_tcor_none',
    tcor4 = 'f_tcor_none',
    tcor5 = 'f_tcor_none',
    tcor6 = 'f_tcor_arrhenius',
    tcor7 = 'f_tcor_none',
    tcor8 = 'f_tcor_none',
    tcor9 = 'f_tcor_none'
  ),

  # - water scalars 
  wcor = list(             
    wcor1 = 'f_wcor_ghezzehei_diffusion', 
    wcor2 = 'f_wcor_ghezzehei_diffusion', 
    #APW: wcor2 = 'f_wcor_ghezzehei_biological',
    wcor3 = 'f_wcor_none',
    wcor4 = 'f_wcor_none', 
    wcor5 = 'f_wcor_ghezzehei_diffusion',
    wcor6 = 'f_wcor_ghezzehei_biological',
    wcor7 = 'f_wcor_ghezzehei_diffusion',
    wcor8 = 'f_wcor_ghezzehei_diffusion',
    wcor9 = 'f_wcor_ghezzehei_diffusion'
  ),

  # - soil or litter scalars 
  scor = list(             
    scor1 = 'f_scor_none'
  ),
 
  # turnover rate functions
  k = list(
    k7 = 'f_k_abramoff_sorp_efficiency' 
  ),

  # vmax and km functions for Michaelis-Menten and reverse M-M fluxes
  vmax = list(NA),

  km = list(NA),
 
  # transfer list
  # - a set of functions that calculate the transfer coefficients from one pool to another
  # - numbers in names refer to pools in order of decomp lsit and state  
  # - the first number in the name is the pool 'from' which the flux is coming 
  # - the second number in the name is the pool 'to' which the flux is going 
  transfer = list(
    t1_to_4 = 'f_transfer_fluxsum_prop_one',
    t1_to_5 = 'f_transfer_fluxsum_prop_two', 
    t2_to_3 = 'f_transfer_fluxid_cue', 
    t2_to_4 = 'f_transfer_fluxid_cue_remainder', 
    t3_to_4 = 'f_transfer_fluxsum_prop_one', 
    t3_to_5 = 'f_transfer_fluxsum_prop_two', 
    t4_to_2 = 'f_transfer_fluxsum_prop_one_cue', 
    t4_to_3 = 'f_transfer_fluxsum_prop_two',
    t5_to_1 = 'f_transfer_fluxid_cue', 
    t5_to_3 = 'f_transfer_fluxid_cue_remainder' 
  ),

  # functions to calculate use/transfer efficiencies 
  cue = list(
    cue6 = 'f_cue_temp_mod_abramoff' 
  ),

  cue2 = list(NA),
  cue3 = list(NA),

  # additional functions 
  texture_mod   = 'f_scor_century_texture',
  quality_mod   = 'f_quality_mod_wieder',
  anpp_mod      = 'f_anpp_mod_wieder',
  matpot_sat    = 'f_matpot_sat_elm_cosby1984_tab5',
  conv_water_vwc_to_mp = 'f_conv_water_vwc_to_mp_vanGenuchten'
)


# environment
####################################
soil_decomp_object$env <- list(
  # litter input, ANPP, and litter quality  
  litter     = 172.8978/365/1000, # litter input     (g C cm-3 soil)   
  anpp       = 0,            # aboveground NPP  (??)
  lignin     = 0,            # litter lignin fraction (proportion)
  N          = 0,            # litter N content (g N per ???)   
  # physical environment
  temp       = 11.21961,     # soil temperature (oC) 
  vwc        = 0.2422044,    # soil volumetric water content (proportion)
  matpot     = 15,           # soil matric potential, should be linked to vwc default for MILLENNIALv2   # APW: units? doesn't look like MPa, but that's what I'm using it for for ELM for now 
  clay       =  0,           # clay (unitless, proportion)
  sand       = 0.2,          # sand (unitless, proportion)
  porosity   = 0.6,          # soil volumetric water content at saturation (proportion)
  matpot_min = 10,           # minimum matric potential, i.e. wilting point (MPa maybe?) 
  pH         = 7,            # soil pH (pH units) 
  depth      = 0             # soil layer depth (cm) APW: increment? lower bound?
  #BD         = 1000,         # bulk density (mg cm-3) 
)


# state parameters (i.e. calculated parameters)
####################################

# state parameter names (& therefore also fnames) that have a value per pool
soil_decomp_object$pool_state_pars    <- c('poolmax')

# state parameter names (& therefore also fnames) that have a value per outflux
soil_decomp_object$outflux_state_pars <- c('k', 'vmax', 'km', 'cue', 'cue2', 'cue3' ) 

soil_decomp_object$state_pars <- list(

  # various solver outputs
  solver_out             = matrix(1),           # output from solver
  solver_steadystate_out = matrix(1),           # output from steady state solver
  transfer_matrix        = matrix(1),           # transfer matrix
  flux_matrix            = matrix(1),           # matrix of internal fluxes from/to pools (diagonal is gross loss from pool, off-diag is flux from col pool to row pool, row-sums are net change in pool from internal fluxes i.e. not including inputs but including outputs) 
  previous_state         = matrix(1),           # previous timestep state
  mass_balance           = numeric(1),          # result from mass balance calculation

  # state pars lists of calculated parameters
  # - see above fnames lists of the same name for description
  poolmax     = list(NA),
  k           = list(NA),
  vmax        = list(NA),
  km          = list(NA),
  cue         = list(NA),
  cue2        = list(NA),
  cue3        = list(NA),

  # additional calculated parameters
  texture_mod = numeric(1),          # intermediate for texture modifier functions 
  quality_mod = numeric(1),          # intermediate for litter quality modifier functions, e.g. for k in MIMICs, Wieder et al. DATE 
  anpp_mod    = numeric(1),          # intermediate for ANPP modifier functions, e.g. for k in MIMICs, Wieder et al. DATE  
  matpot_sat  = numeric(1)           # max soil matric potential 
)


# parameters
####################################

# parameter names that have a value per pool
soil_decomp_object$pool_pars <- 
  c(soil_decomp_object$pool_state_pars, 'cstate0', 'input_coef' )

# parameter names that have a value per outflux
soil_decomp_object$outflux_pars <- 
  c(soil_decomp_object$outflux_state_pars, 
  'k_exp', 'vmax2', 'km2', 'ea', 'cat_pool', 'cue_texture_exp', 'cue_quality_exp' )

# parameters
soil_decomp_object$pars <- list(
  
  # model structure parameters
  pool_init_value  = 0.1,      # initialisation value for all state pools (same units as state pools)
  n_pools          = 5,        # number of pools in model  # APW: need to rebuild model when this changes in wrapper runs, likely other issues
  n_outfluxes      = 9,        # number of outfluxes (arrow bases) in model 
  # APW: sat_pool might need to be n_outfluxes long if there are several
  sat_pool         = 3,        # pool which saturates
  leaching_outflux = 8,        # ouflux which goes to leaching 
  error_tolerance  = 1e-6,     # error tolerance for mass balance check (same units as pools, mg C g-1 soil) 

  # general scalar parameters
  R                = 8.31446,  # universal gas constant (J K-1 mol-1)
  conv_mm_to_MPa   = -9.8e-6,  # convert water pressure from mm to MPa (MPa mm-1)
  beta             = 2,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  # used in CORPSE
  clayref          = NA,
  minmic           = NA,
  qslope_exp       = NA,
  # used in MEND 2020
  mr               = NA,     
  Kads             = NA,       
  # used in Millenial v2
  max_poolmax      = 86.0,     # maximum pool size when soil is 100 % clay (same units as pools)
  sorp_exp_a       = -0.216,   # exponent parameter intercept for calculating sorption k from pH (unitless)
  sorp_exp_b       = -0.12,    # exponent parameter slope for calculating sorption k from pH (1/pH units)
  # used in MIMICS
  quality_a        = 0.85,     # intercept parameter for quality modifier 
  quality_b        = 0.013,    # slope parameter for quality modifier 
  quality_c        = 0.5,      # scalar parameter for quality modifier 
  anpp_mod_p1      = NA,       #
  anpp_mod_p2      = NA,       #
  anpp_mod_p3      = NA,       #
  anpp_mod_p4      = NA,       #
  km_temp_slope    = 0.02,
  km_temp_int      = 3.19,
  km_texture_norm  = 3,
  km_texture_exp   = -2,
  k_desorp_tuning  = 0.1,
  km_tuning        = 0.15625,
  cue_texture_tuning = 0.1,
  cue_quality_tuning = 1,

  # temperature correction parameters
  reftemp                = 15,      # reference temperature for paramter to be scaled (oC) 
  reftemp_cue            = 15,      # reference temperature for CUE scaling (oC)
  q10                    = NA,      # Q10 (factor of increase per 10 oC)
  cue_temp_scalar        = 0.012,   # CUE scaling parameter for temp difference (oC-1) 
  tcor_mod               = 1,       # modifier on the scalar returned from the Q10 function (0-1, unitless), APW: could be moved into outflux calc function
  tcor_exp_norm          = 3e-06,   # normalisation constant for exponential scaling, MIMICS values (unitless)
  tcor_exp_exp_int       = 5.47,    # intercept linear relationship to temp in exponenet for exponential scaling, MIMICS values (unitless) 
  tcor_exp_exp_slope     = 0.063,   # slope linear relationship to temp in exponenet for exponential scaling, MIMICS values (oC-1) 
  tcor_century1          = 15.4,    # parameter 1 for CENTURY type logistic/saturating exponential scaling (oC) 
  tcor_century2          = 11.75,   # parameter 2 for CENTURY type logistic/saturating exponential scaling (unitless) 
  tcor_century3          = 29.7,    # parameter 3 for CENTURY type logistic/saturating exponential scaling (unitless) 
  tcor_century4          = 0.031,   # parameter 4 for CENTURY type logistic/saturating exponential scaling (unitless) 
  
  # water correction parameters
  wcor_century1          = 30,      # parameter 1 for CENTURY type inverse exponential scaling (oC) 
  wcor_century2          = 9,       # parameter 2 for CENTURY type inverse exponential scaling (unitless) 
  wcor_unimodal_exp      = 3,       # exponent for Sulman unimodal function
  wcor_unimodal_exp_diff = -0.5,    # difference in exponents for Sulman unimodal function
  wcor_diffusion_exp     = 0.5,     # exponent for Ghezzehei et al. 2018 diffusion power law 
  wcor_bio_exp           = 0.5,     # exponent for Ghezzehei et al. 2018 biological water limitation 
  wcor_bio_lambda        = 2.1e-4,  # dependence of rate on matric potential for Ghezzehei et al. 2018 biological water limitation 
  wcor_bio_kamin         = 0.2,     # minimum relative rate in saturated soil for Ghezzehei et al. 2018 biological water limitation 

  # soil property correction parameters
  scor_texture_int       = 0.85,    # intercept linear relationship to sand, from CENTURY (unitless)    
  scor_texture_slope     = -0.68,   # slope linear relationship to sand, from CENTURY (unitless)    
  scor_quality_exp       = -3,      # exponent exponential relationship to lignin, from CENTURY (unitless)    
  
  
  # Pool-specific parameters
  # - lists 1:n_pools long 
  ################################################
  
  # outfluxes to decomp pool assignment
  decomp_outflux1 = list(dof11=1, dof12=2 ),
  decomp_outflux2 = list(dof21=3),
  decomp_outflux3 = list(dof31=4, dof32=5 ),
  decomp_outflux4 = list(dof41=6, dof42=7, dof43=8),
  decomp_outflux5 = list(dof51=9),
  
  # initial pool mass for each pool
  # APW: set automatically during build, could modify to use this value
  cstate0 = list(NA),

  # input coefficients (allocates proportions of inputs into different pools)
  input_coef = list(
    input_coef1 = .66, 
    input_coef2 = 0, 
    input_coef3 = 0,
    input_coef4 = .34,
    input_coef5 = 0
  ),
  
  # maximum size for pool i 
  poolmax = list(NA),
  

  # Decomposition / Outflux parameters 
  # - these lists are n_outfluxes long and when there is just a single outflux per pool this collapses to n_pools long
  ################################################ 

  # turnover rates for linear and some non-linear decomposition 
  # -- linear decomp units are: day-1
  # -- however, some more funky decomp functions require different units, e.g. Giorgiou DD function: g soil mg-1 C day-1 
  k = list(    
    k1 = NA,
    k2 = 0.018,        
    k3 = 4.5,          
    k4 = 1e-3,
    k5 = 4.8e-3,  
    k6 = NA,     
    k7 = 1.0, 
    k8 = 1.5e-3,
    k9 = 0.02    
  ),

  # exponent in k calculation Wieder
  k_exp = list(),  

  # catalyst pool for a given outflux
  cat_pool = list(
    cat_pool1 = 2, 
    cat_pool2 = NA, 
    cat_pool3 = NA,
    cat_pool4 = NA,
    cat_pool5 = NA,
    cat_pool6 = 2,
    cat_pool7 = NA,
    cat_pool8 = NA, 
    cat_pool9 = NA
  ),
  
  # max rates for outflux 
  vmax = list(
    vmax1 = 4.680971,  
    vmax2 = NA,       
    vmax3 = NA,
    vmax4 = NA,
    vmax5 = NA,
    vmax6 = 74.54208, 
    vmax7 = NA,
    vmax8 = NA,
    vmax9 = NA
  ),
  
  # max rates for outflux when two fluxes are combined 
  vmax2 = list(NA),
 
  # activiation energy for temp sensitivity of flux 
  ea = list(
    ea1 = 6.3909e+04, 
    ea2 = NA,
    ea3 = NA,
    ea4 = NA,
    ea5 = NA,
    ea6 = 5.7865e+04, 
    ea7 = NA,
    ea8 = NA,
    ea9 = NA
  ),
 
  # Michaelis-Menten half-saturation constant for flux      
  km = list(   
    km1 = 6.443,  
    km2 = NA,
    km3 = NA,
    km4 = NA,
    km5 = NA,
    km6 = 0.7746,
    km7 = NA,
    km8 = NA,
    km9 = NA
  ),
  
  # Michaelis-Menten half-sat const for outflux when two fluxes are combined 
  km2 = list(NA),
  
  
  # Transfer parameters 
  # carbon use or transfer efficiency for outflux 
  # - if this varies by the 'to' pool we need another function / parameters
  cue = list(
    cue1 = NA,       
    cue2 = NA,
    cue3 = 0.5,      
    cue4 = NA,
    cue5 = NA, 
    cue6 = 0.19,
    cue7 = NA,
    cue8 = NA,
    cue9 = 0.33
  ),  
  
  # additional carbon use or transfer efficiency for outflux 
  cue2 = list( 
    cue1 = NA,       
    cue2 = NA,       
    cue3 = NA,
    cue4 = NA,
    cue5 = NA,
    cue6 = NA,
    cue7 = NA,
    cue8 = NA,
    cue9 = NA
  ),  

  # another additional carbon use or transfer efficiency for outflux 
  cue3 = list(NA), 
  cue_texture_exp = list(NA),      # multiplier in exponent of CUE calculation as a function of clay, Wieder
  cue_quality_exp = list(NA)       # multiplier in exponent of CUE calculation as a function of litter quality, Wieder
)



# state
####################################
soil_decomp_object$state <- list(
  # cpools and outflux extents are set during initialization
  cpools      = matrix(1),
  outflux     = list(NA), 
  respiration = numeric(1), 
  leaching    = numeric(1) 
)



# run control parameters
####################################
soil_decomp_object$cpars <- list(
  verbose  = F,          # write diagnostic output during runtime 
  cverbose = F,          # write diagnostic output from configure function 
  output   = 'pools'     # type of output from run function
)



# output functions
#######################################################################        

f_output_soil_decomp_eval <- f_output_eval 


f_output_soil_decomp_run <- function(.) {
  unlist(.$state)
}

f_output_soil_decomp_pools <- function(.) {
  c(unlist(.$state$cpools),unlist(.$state$respiration))
}

f_output_soil_decomp_state <- function(.) {
  unlist(.$state)
}

f_output_soil_decomp_full <- function(.) {
  c(unlist(.$state),unlist(.$state_pars))
}



# test functions
#######################################################################        

soil_decomp_object$.test <- function(., 
  verbose=F, metdf=F, ntimes=100, sigfig=3,
  steadystate=F, 
  litter=.00016, anpp=500*24 
  ) {

  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  if(steadystate) .$fnames$sys <- 'f_steadystate_npools'
  .$configure_test() # if only used in test functions should begin with a .

  if(metdf) {
    .$dataf       <- list()
    if(length(litter)==1) litter <- rep(litter, ntimes ) 
    .$dataf$metdf <- matrix(litter, nrow=1 )
    rownames(.$dataf$metdf) <- 'soil_decomp.litter'  
    .$dataf$lm    <- dim(.$dataf$met)[2]
    .$dataf$mout  <- .$output()
    print('')  
    signif(.$run_met(), sigfig )
  } else {
    .$env$litter  <- litter
    .$env$anpp    <- anpp 
    print('')  
    signif(.$run(), sigfig )
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
        tcor = list(
          t1 = 'f_tcor_wieder'
        )
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


soil_decomp_object$.test_ctc <- function(., verbose=F, metdf=F, 
                                         litter=1, ntimes=1 ) {
  
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  
  .$fnames <- list(
    sys                = 'f_sys_npools',
    solver             = 'plsoda',
    solver_func        = 'f_solver_func_elmv2ctc',
    # input              = 'f_input_mimics',
    steadystate        = 'f_steadystate_npools',
    solver_steadystate = 'pstode',
    decomp = list(
      d1 = 'f_decomp_lin', 
      d2 = 'f_decomp_lin', 
      d3 = 'f_decomp_lin', 
      d4 = 'f_decomp_lin', 
      d5 = 'f_decomp_lin',
      d6 = 'f_decomp_lin',
      d7 = 'f_decomp_lin'
    )
  )
  
  .$pars <- list(
    n_pools = 7,          # number of pools in model  
    
    # initial pool mass for each pool
    # values are equilibrium calculated with stode function in DeSolve script
    cstate0 = list(
      cstate01 = 0.1, 
      cstate02 = 0.1,          
      cstate03 = 0.1,      
      cstate04 = 0.1,       
      cstate05 = 0.1,
      cstate06 = 0.1,
      cstate07 = 0.1
    ),
    
    # Carbon use efficiency from pool i r-microbes (cue) or k-microbes (cue2)
    cue = list(
      cue1 = 1.0,       
      cue2 = 0.61,  
      cue3 = 0.45,
      cue4 = 0.71,
      cue5 = 0.72,
      cue6 = 0.54,
      cue7 = 0
    ),  

    # turnover rate for linear decomposition
    k = list(
      k1 = .001,
      k2 = .7,
      k3 = .07,
      k4 = .014,
      k5 = .07,
      k6 = .014,
      k7 = .0005
    )
  )
  
  .$env <- list(
    litter = 1
  )
  
  .$configure_test() # if only used in test functions should begin with a .
  
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
