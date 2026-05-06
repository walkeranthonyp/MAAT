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
  #sys                   = 'f_sys_npools', 
  sys                   = 'f_steadystate_npools', # if f_steadystate_npools selected, solver_steadystate automatically used  
  solver                = 'plsoda',
  solver_func           = 'f_solver_func_soilR_outfluxes',
  steadystate           = 'f_steadystate_npools', # options: 'f_steadystate_null','f_steadystate_npools',
  solver_steadystate    = 'pstode',
  input                 = 'f_input',
  DotO                  = 'f_DotO',
  transfermatrix        = 'f_transfermatrix',
  transfer_fluxsum_prop = 'f_transfer_fluxsum_prop',
  outfluxes             = 'f_outfluxes',
  
  # decay/decomposition/flux functions, one per pool
  # - these functions represent the total flux out of a given pool
  # - can be the sum of multipl terms in the "outflux" list
  # - pars$decomp_outflux{1-n} lists for ouflux assignments to each pool 
  decomp = list(
    d1 = 'f_decomp_fluxsum',
    d2 = 'f_decomp_fluxsum',
    d3 = 'f_decomp_fluxsum',
    d4 = 'f_decomp_fluxsum',
    d5 = 'f_decomp_fluxsum',
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
 
  # DOC uptake, another flux term 
  docuptake  = 'f_decomp_mm', #MM uptake of DOC
 
 
  # output fluxes from pools functions 
  # - one function per flux, designed for when there are more than one output flux from a pool
  outflux = list(
    of1 = 'f_decomp_lin',
    of2 = 'f_decomp_rmm',
    of3 = 'f_decomp_dd_georgiou',  
    of4 = 'f_desorp_millennialv2', 
    of5 = 'f_decomp_lin', 
    of6 = 'f_decomp_lin',
    of7 = 'f_sorp_sat',
    of8 = 'f_decomp_mm',
    of9 = 'f_decomp_lin'
  ),

  # currently use non-list tcor structure for model-specific solvers for millennial or models with only one scalar
  # - will be deprecated when generic solver used for all models   
  # - tcor = 'f_tcor_wieder', 
  # - tcor = 'f_tcor_daycent2_abramoff', #use this structure for millennial or models with only one scalar #century

  # functions that align with each outflux 
  # - temperature scalars 
  tcor = list(             
    t1 = 'f_tcor_none',
    t2 = 'f_tcor_arrhenius_millennialv2',
    t3 = 'f_tcor_none',
    t4 = 'f_tcor_none',
    t5 = 'f_tcor_none',
    t6 = 'f_tcor_none',
    t7 = 'f_tcor_none',
    t8 = 'f_tcor_arrhenius_millennialv2',
    t9 = 'f_tcor_none',
    t10 = 'f_tcor_none'
  ),

  # - water scalars 
  wcor = list(             
    w1 = 'f_wcor_ghezzehei_diffusion', 
    w2 = 'f_wcor_ghezzehei_diffusion', 
    #w2 = 'f_wcor_ghezzehei_biological',
    w3 = 'f_wcor_none',
    w4 = 'f_wcor_none', 
    w5 = 'f_wcor_ghezzehei_diffusion',
    w6 = 'f_wcor_ghezzehei_diffusion',
    w7 = 'f_wcor_ghezzehei_diffusion',
    w8 = 'f_wcor_ghezzehei_biological',
    w9 = 'f_wcor_ghezzehei_diffusion',
    w10 = 'f_wcor_ghezzehei_diffusion'
  ),

  # - soil or litter scalars 
  scor = list(             
    s1 = 'f_scor_none',
    s2 = 'f_scor_none',
    s3 = 'f_scor_none',
    s4 = 'f_scor_none',
    s5 = 'f_scor_none',
    s6 = 'f_scor_none',
    s7 = 'f_scor_none',
    s8 = 'f_scor_none',
    s9 = 'f_scor_none',
    s10 = 'f_scor_none'
  ),

  # transfer list
  # - a set of functions that calculate the transfer coefficients from one pool to another
  # - numbers in names refer to pools in order of decomp lsit and state  
  # - the first number in the name is the pool 'from' which the flux is coming 
  # - the second number in the name is the pool 'to' which the flux is going 
  transfer = list(
    t1_to_4 = 'f_transfer_fluxsum_prop_two',
    t1_to_5 = 'f_transfer_fluxsum_prop_one', 
    t2_to_3 = 'f_transfer_cue', 
    t2_to_4 = 'f_transfer_cue_remainder', 
    t3_to_4 = 'f_transfer_fluxsum_prop_one', 
    t3_to_5 = 'f_transfer_fluxsum_prop_two', 
    t4_to_2 = 'f_transfer_fluxsum_prop_three_cue', 
    t4_to_3 = 'f_transfer_fluxsum_prop_two',
    #t4_to_resp = uptake 'frac of three out fluxes from pool 4 * (1- CUE)',
    #t4_to_leached = leached 'frac of three out fluxes from pool 4',
    t5_to_1 = 'f_transfer_cue', 
    t5_to_3 = 'f_transfer_cue_remainder' 
    # t1_to_3 = 'f_transfer_cue_remainder',
    # t1_to_4 = 'f_transfer_cue_remainder',
    # t2_to_5 = 'f_transfer_mend27',
    # t3_to_6 = 'f_transfer_all',
    # t4_to_7 = 'f_transfer_all',
    # t5_to_6 = 'f_transfer_mend52',
    # t6_to_7 = 'f_transfer_all',
    # t7_to_8 = 'f_transfer_all',
    # t2_to_9 = 'f_transfer_mend27',
    # t3_to_9 = 'f_transfer_all',
    # t4_to_9 = 'f_transfer_all',
    # t5_to_9 = 'f_transfer_mend52',
    # t6_to_9 = 'f_transfer_all',
    # t7_to_9 = 'f_transfer_all',
    # t8_to_9 = 'f_transfer_all'
    # t1_to_3 = 'f_transfer_cue_remainder',
    # t1_to_5 = 'f_transfer_cue',
    # t2_to_1 = 'f_transfer_mend21',
    # t2_to_5 = 'f_transfer_mend25',
    # t2_to_6 = 'f_transfer_mend26',
    # t2_to_7 = 'f_transfer_mend27',
    # t3_to_5 = 'f_transfer_all',
    # t4_to_5 = 'f_transfer_all',
    # t5_to_2 = 'f_transfer_mend52',
    # t5_to_4 = 'f_transfer_mend54',
    # t6_to_5 = 'f_transfer_all',
    # t7_to_5 = 'f_transfer_all'
  ),

  # individual functions not directly associated with a specific pool or flux  
  growthresp = NA,
  maintresp  = NA,
  water_unit_converter = 'f_SWC2SWP_vanGenuchten'
)


# environment
####################################
soil_decomp_object$env <- list(
  litter   = 172.8978/365, # forc_npp in MILLENNIALv2
  temp     = 11.21961,     # default for MILLENNIALv2 (sum across mean year)
  vwc      = .2422044,     # default for MILLENNIALv2 (average across mean year)
  porosity = 0.6,          # default for MILLENNIALv2
  clay     = .0, 
  lignin   = 0,
  N        = 0,
  anpp     = 0,
  depth    = 0,
  pH       = 7,            # default for MILLENNIALv2
  BD       = 1000,         # default for MILLENNIALv2   # Bulk density 
  matpot   = 15,           # default for MILLENNIALv2
  lambda   = 2.1000e-04,   # default for MILLENNIALv2
  kamin    = .2,           # default for MILLENNIALv2
  claysilt = 80            # default for MILLENNIALv2   # MILLENNIAL uses clay+silt% to calculate Qmax
)


# state parameters (i.e. calculated parameters)
####################################
soil_decomp_object$state_pars <- list(
  solver_out             = matrix(1),
  solver_steadystate_out = matrix(1)
)


# parameters
####################################
# parameter names that have a value per pool
# APW: decomp_ouflux needs a way to configure each list
soil_decomp_object$pool_pars <- c('cstate0', 'decomp_outflux', 'cue', 'cue2', 'poolmax', 'input_coefs' )

# parameter names that have a value per outflux
soil_decomp_object$outflux_pars <- c('k', 'vmax', 'vmax2', 'km', 'km2', 'ea' ) 

# parameters
soil_decomp_object$pars <- list(
  
  # model structure parameters
  n_pools      = 5,        # number of pools in model  #need to change and re-create XMLs when this changes for wrapper runs
  n_outfluxes  = 9,        # number of pools in model  #need to change and re-create XMLs when this changes for wrapper runs
  # APW: these might need to be n_outfluxes long
  cat_pool     = 2,        # pool which catalyses reactions 
  sat_pool     = 3,        # pool which saturates

  # general scalar parameters
  # APW: looks like some of these might be env vars
  beta         = 2,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  silt         = NA,      
  clay         = NA,        
  clayref      = NA,
  mr           = NA,     
  pep          = NA,        
  pem          = NA,        
  Kads         = NA,       
  qslope_mayes = NA,
  reftemp      = 30,       # CENTURY
  R            = 8.31446,  # used in MILLENNIALv2
  minmic       = NA,
  q10          = NA,

  
  # Pool-specific parameters
  # - lists 1:n_pools long 
  ################################################
  
  # outfluxes to decomp pool assignment
  # APW: list format didn't work due to nesting lists within a list, configure functions cannot handle 
  decomp_outflux1 = list(dof11=1, dof12=2),
  decomp_outflux2 = list(dof21=1, dof22=2),
  decomp_outflux3 = list(dof31=1, dof32=2),
  decomp_outflux4 = list(dof41=7),
  decomp_outflux5 = list(dof51=8),
  decomp_outflux6 = list(dof61=9),
  decomp_outflux7 = list(dof71=10),
#  decomp_outflux1 = list(dof11=1),
#  decomp_outflux2 = list(dof21=3),
#  decomp_outflux3 = list(dof31=4),
#  decomp_outflux4 = list(dof41=6),
#  decomp_outflux5 = list(dof51=9),
#  decomp_outflux1 = list(dof11=1, dof12=2),
#  decomp_outflux2 = list(dof21=3),
#  decomp_outflux3 = list(dof31=4, dof32=5),
#  decomp_outflux4 = list(dof41=6, dof42=7, dof43=8),
#  decomp_outflux5 = list(dof51=9),
#  decomp_outflux = list(
#    decomp_outflux1 = list(dof1.1=1, dof1.2=2),
#    decomp_outflux2 = list(dof2.1=3),
#    decomp_outflux3 = list(dof3.1=4, dof3.2=5),
#    decomp_outflux4 = list(dof4.1=6, dof4.2=7, dof4.3=8),
#    decomp_outflux5 = list(dof5.1=9),
#  ),
  
  # initial pool mass for each pool
  cstate0 = list( 
    cstate01 = 1, 
    cstate02 = 1,          
    cstate03 = 1,      
    cstate04 = 1,       
    cstate05 = 1       
  ),

  # input coefficients (allocates proportions of inputs into different pools)
  input_coefs = list(
    input_coefs1 = .66, 
    input_coefs2 = 0, 
    input_coefs3 = 0,
    input_coefs4 = .34,
    input_coefs5 = 0
  ),
  
  # maximum size for pool i 
  poolmax = list(       
    poolmax1 = NA,      
    poolmax2 = NA,       
    poolmax3 = 68.8, 
    poolmax4 = NA,
    poolmax5 = NA
  ),
  

  # Decomposition / Outflux parameters 
  # - these lists are n_outfluxes long and when there is just a single outflux per pool this collapses to n_pools long
  ################################################ 

  # turnover rates for linear and some non-linear decomposition 
  # -- linear decomp units are: day-1
  # -- however, some more funky decomp functions require different units, e.g. Giorgiou DD function: g soil mg-1 C day-1 
  k = list(    
    k1 = 0.018,       # aggregate formation from POM       
    k2 = NA,
    k3 = 4.5,         # microbial turnover rate   
    k4 = NA,
    #k5 = 4.8000e-03,  # aggregate formation from maom 
    #k6 = 0.0015,      # leaching rate (kd)
    k7 = NA, 
    k8 = NA,
    k9 = 0.02,         # aggregrate breakdown rate (kb)
    k10 = 0.02         # aggregrate breakdown rate (kb)
  ),

  # max turnover rates per unit microbial biomass for pool i
  # commented out old list (updating to allow for multiple mic groups or catalysts)
  vmax = list(
    vmax1 = NA,       # Vm
    vmax2 = 1.8e+12,  # alpha_pl; pre-exponential constant for temp sensitivity of POM decay
    vmax3 = NA,
    vmax4 = NA,
    vmax5 = NA,
    vmax6 = NA,
    vmax7 = NA,
    vmax8 = 2.3e+12, # alpha_lb #pre-exponential constant for temp sensitivity of DOC uptake
    vmax9 = NA,
    vmax10 = NA
  ),
  
  vmax2 = list(
    vmax1 = NA, 
    vmax2 = NA,
    vmax3 = NA,
    vmax4 = NA,
    vmax5 = NA,
    vmax6 = NA,
    vmax7 = NA,
    vmax8 = NA,
    vmax9 = NA,
    vmax10 = NA
  ),
  
  ea = list(
    ea1 = NA,
    ea2 = 6.3909e+04, # MILLENNIALv2 #activiation energy for temp sensitivity of POM decay
    ea3 = NA,
    ea4 = NA,
    ea5 = NA,
    ea6 = NA,
    ea7 = NA,
    ea8 = 5.7865e+04, # MILLENNIALv2 #activiation energy for temp sensitivity of DOC uptake
    ea9 = NA,
    ea10 = NA
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
    km1 = NA,
    km2 = 6.443,  # MillennialV2, half sat const for POM breakdown 
    km3 = NA,
    km4 = NA,
    km5 = NA,
    km6 = NA,
    km7 = NA,
    km8 = 0.7746, # DOC uptake half sat constant
    km9 = NA,
    km10 = NA
  ),
  
  km2 = list( 
    km1 = NA, 
    km2 = NA,     
    km3 = NA,
    km4 = NA,
    km5 = NA,
    km6 = NA,
    km7 = NA,
    km8 = NA,
    km9 = NA,
    km10 = NA
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
    rkm7 = NA,
    rkm8 = NA,
    rkm9 = NA,
    rkm10 = NA
  ),   

  
  # Transfer parameters 
  # Carbon use or transfer efficiency from pool i to any another 
  # - if this varies by the 'to' pool we need another function / parameters
  # MC: I've created a function cue_remainder that can allocate the remainder to another pool rather than CO2
  # APW: these lists should probably be n_outfluxes long 
  cue = list(
    cue1 = NA,       
    cue2 = 0.5,
    cue3 = NA,      
    cue4 = 0.19, # what is this for in Mv2 context?
    cue5 = 0.33 
  ),  
#  cue = list(
#    cue1 = NA,       
#    cue2 = NA,
#    cue3 = 0.5,      
#    cue4 = NA,
#    cue5 = NA, 
#    cue6 = NA,
#    cue7 = NA,
#    cue8 = 0.19,
#    cue9 = 0.33
#  ),  
  
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


  # Model specific parameters
  ################################################

  # MIMICS-specific parameters
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
  
  # MEND-specific parameters
  mend = list(
    Mr =  .00028,
    Pep = .01,
    Pem = .01,
    Gd =  .5,
    Fd =  .5
  ),
  
  # MILLENNIAL-specific parameters
  # v1
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
  
  # v2
  millennialV2 = list(
    param_pc = .86, # slope of mineral C - clay relationship from Georgiou et al. in review
    kld = 1,        # desorption coefficient
    sorp_p1 = .12,  # sorption affinity parameter
    sorp_p2 = .216, # sorption affinity parameter
    cue_t = 0.012,  # slope of CUE temp sensitivity
    Taeref = 15,    # ref temp for CUE temp-dependence equation
    pa = 0.33,      # proportion agg breakdown into pom
    param_pb = .5   # fraction of mb turnover to maom vs doc
  ),
  
  # CENTURY-specific parameters
  century = list(
    c1= 0.85, #century constant in texture function
    c2= 0.68, #century constant in texture function
    slow_to_active= 0.42, #century transfer coefficient
    slow_to_passive= 0.03,#century transfer coefficient
    passive_to_active= 0.45,#century transfer coefficient
    active_to_passive= 0.004,#century transfer coefficient
    metlitter_to_active= 0.45,#century transfer coefficient
    strlitter_to_active= 0.5,#century transfer coefficient
    strlitter_to_slow= 0.7#century transfer coefficient
  )
)

# state
####################################
soil_decomp_object$state <- list(
  #cpools      = matrix(1:5, ncol=1 ),     # APW: would ideally be 1:n_pools 
  cpools      = matrix(1:soil_decomp_object$pars$n_pools, ncol=1 ),   
  respiration = vector("double",1),
  # APW: must also be set during initialization
  outflux     = as.list(numeric(soil_decomp_object$pars$n_outfluxes)) 
#  outflux     = list(   # APW: would ideally use 1:n_oufluxes to set up  
#    of1 = numeric(1),
#    of2 = numeric(1),
#    of3 = numeric(1),
#    of4 = numeric(1),
#    of5 = numeric(1),
#    of6 = numeric(1),
#    of7 = numeric(1),
#    of8 = numeric(1),
#    of9 = numeric(1)
#  )
)
names(soil_decomp_object$state$outflux) <- paste0('of',1:soil_decomp_object$pars$n_outfluxes)



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
    signif(.$run(), 3 )
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
    
    input_coefs = list(
      input_coef1 = 0.1,
      input_coef2 = 0.2,
      input_coef3 = 0.5,
      input_coef4 = 0.3,
      input_coef5 = 0.5,
      input_coef6 = 0.5,
      input_coef7 = 0
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
