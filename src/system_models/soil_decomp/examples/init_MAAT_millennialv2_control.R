################################
#
# MAAT soil_decomp model - initialisation script
# 
# AWalker (walkerap@ornl.gov) 
# February 2018
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script initialises the soil_decomp version of MAAT
# - by setting the values of the init_static and init_dynamic lists
# - each of these lists elements are sublists named with the names of the model objects to which the variables in the sublists belong  
# - each of these model object sublist elements are sublists named with the type of variable they comprise - fnames, pars, or env  
# - each of these variables sublist elements are sublists named with the variable to which they refer in the model object  

# for example, to setup a leaf model object simulation:
# set static variables during runtime (init_static) by setting:
# - fnames.static
# - pars.static     
# - env.static


# set dynamic variables during runtime (init_dynamic) by setting:
# - fnames.var
# - pars.var     
# - env.var
# for a process SA two more sublists need to be added as elements to the init_dynamic$leaf list:
# - pars_proc and pars_eval, both with the same named elements
#   - pars_proc is a list of strings that associate each parameter with the process it is associated with, these should be one of the element names in the list init_dynamic$X$fnames
#   - pars_eval is a string that once evaluated gives a vector of parameter samples  


# This script must set the values of the above 6(8) lists, 
# if you do not want to set any values using these lists they must be specified as NA

###################################################################



### Define the functions, parameters, and environment variables that are to remain static through a simulation 
###############################

# define lists
fnames.static <- list(
  
  sys            = 'f_sys_npools',
  solver         = 'plsoda',
  solver_func    = 'f_solver_func_millennialV2',
  input          = 'f_input',
  steadystate        ='f_steadystate_npools', # options: 'f_steadystate_null','f_steadystate_npools',
  solver_steadystate = 'pstode',
  
  decomp = list(
    d1 = 'f_decomp_rmm',
    d2 = 'f_decomp_dd_georgiou',  
    d3 = 'f_zero', 
    d4 = 'f_decomp_lin', 
    d5 = 'f_decomp_lin'
  ),
  
  desorp = list(
    ds3 = 'f_desorp_millennialv2'
  ),
  
  sorp = list(
    s4 = 'f_sorp_sat'
  ),
  
  aggform = list( 
    a1 = 'f_decomp_lin', #Aggregate formation from POM
    a3 = 'f_decomp_lin'  #Aggregate formation from MAOM
  ),
  
  docuptake = 'f_decomp_mm', #MM uptake of DOC
  
  wcor = 'f_wcor_ghezzehei_diffusion',
  
  wcor2 = 'f_wcor_ghezzehei_biological',
  
  tcor = list(                  #this structure for corpse
    t1 = 'f_tcor_arrhenius_millennialv2',
    t2 = NULL,
    t3 = NULL,
    t4 = 'f_tcor_arrhenius_millennialv2',
    t5 = NULL
  )
)

pars.static <- list(
  #converted all parameter values to yearly rates
  
  
  n_pools = 5,        # number of pools in model
  beta    = 2,        # density dependent turnover, biomass exponent (can range between 1 and 2)
  R = 8.31446,       #used in MILLENNIALv2
  
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
    input_coef5 = 0
  ),
  
  # initial pool mass for each pool
  # values from 1000-yr spinup in desolve version of model (i.e. external to MAAT)
  cstate0 = list( 
    cstate01 = 1, 
    cstate02 = 1,          
    cstate03 = 1,      
    cstate04 = 1,       
    cstate05 = 1
  ),
  
  # Carbon use or transfer efficiency from a given pool
  cue = list(
    cue1 = NA,       
    cue2 = NA,       
    cue3 = NA,
    cue4 = .19,
    cue5 = NA
  ),  
  
  # max turnover rate per unit microbial biomass for pool i 
  vmax = list(
    vmax1 = 1.8000e+12,#alpha_pl; pre-exponential constant for temp sensitivity of POM decay
    vmax2 = NA,   
    vmax3 = NA,
    vmax4 = 2.3000e+12, #alpha_lb #pre-exponential constant for temp sensitivity of DOC uptake
    vmax5 = NA
  ),
  
  # michaelis-menten half-saturation constant for microbial decomnp of pool i      
  km = list(   
    km1 = 6443, #MillennialV2 #half sat const for POM breakdown 
    km2 = NA,
    km3 = NA,
    km4 = 774.6, #DOC uptake half sat constant,
    km5 = NA
  ),

  # turnover rate for linear decomposition
  k = list(    
    k1 = .018, #aggregate formation from POM       
    k2 = 4.5000e-03, #microbial turnover rate   
    k3 = 4.8000e-03, #aggregate formation from maom 
    k4 = .0015, #leaching rate (kd)
    k5 = .02 #aggregrate breakdown rate (kb)
  ),
  
  # maximum size for pool i 
  poolmax = list(       
    poolmax1 = NA,      
    poolmax2 = NA,       
    poolmax3 = NA, #this is calculated in solver function 
    poolmax4 = NA,
    poolmax5 = NA
  ),
  
  #millennialV2-specific parameters
  millennialV2 = list(
    param_pc = .86, ##slope of mineral C - clay relationship from Georgiou et al. in review
    kld = 1, #desorption coefficient
    sorp_p1 = .12, #sorption affinity parameter
    sorp_p2 = .216, #sorption affinity parameter
    kld = 1, #desorption coefficient
    cue_t = 0.012, #slope of CUE temp sensitivity
    Taeref = 15,   #ref temp for CUE temp-dependence equation
    pa = 0.33,      #proportion agg breakdown into pom
    param_pb = .5 #fraction of mb turnover to maom vs doc
  )
)

env.static  <- list(
  litter = 172.8978/365, #forc_npp in MILLENNIALv2
  temp   = 11.21961, #default for MILLENNIALv2 (sum across mean year)
  vwc    = .2422044, #default for MILLENNIALv2 (average across mean year)
  porosity = 0.6, #default for MILLENNIALv2
  pH = 7, #default for MILLENNIALv2
  BD = 1000,  #default for MILLENNIALv2 #Bulk density 
  matpot = 15, #default for MILLENNIALv2
  lambda = 2.1000e-04, #default for MILLENNIALv2
  kamin = .2, #default for MILLENNIALv2
  claysilt = 80  #default for MILLENNIALv2 #MILLENNIAL uses clay+silt% to calculate Qmax
)



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NA where variation is not required  

# define lists
fnames.var <- NULL

pars.var <- NULL

pars_proc.var <- NULL

pars_eval.var <- NULL

env.var <- NULL



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  fnames = list( soil_decomp = fnames.static),
  pars   = list( soil_decomp = pars.static),
  env    = list( soil_decomp = env.static)
)

init_dynamic <- list(
  fnames    = list( soil_decomp = fnames.var),
  pars      = list( soil_decomp = pars.var),
  pars_proc = list( soil_decomp = pars_proc.var),
  pars_eval = list( soil_decomp = pars_eval.var),
  env       = list( soil_decomp = env.var)
)



### END ###
