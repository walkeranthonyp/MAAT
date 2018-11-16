################################
#
# MAAT Leaf Model - initialisation script
# 
# AWalker (walkerap@ornl.gov) 
# December 2015
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script initialises the leaf photosynthesis version of MAAT
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
  solver  = 'f_A_r_leaf_analytical',
  gstar   = 'f_scalar_none',
  tcor_asc = list(
  	Kc     = 'f_scalar_none',
  	Ko 		 = 'f_scalar_none',
  	vcmax  = 'f_scalar_none',
  	jmax   = 'f_scalar_none'
  ),
    tcor_des = list(
  	vcmax  = 'f_scalar_none',
  	jmax   = 'f_scalar_none'
  ),
  vcmax = 'f_vcmax_constant',
  jmax  = 'f_jmax_lin',
  rd    = 'f_rd_lin_vcmax', 
  ri    = 'f_r_zero',
  rb    = 'f_r_zero',
  rs    = 'f_rs_medlyn2011'
)

pars.static <- list(
  atref = list(vcmax  = 50),
  diag  = T
  )

env.static  <- list(
  par     = 1000,
  ca_conc = 400
  )



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NA where variation is not required  

# define lists
fnames.var <- list(
  tcor_asc = list(vcmax = c('f_scalar_none','f_tcor_asc_Q10','f_tcor_asc_Arrhenius')),
  tcor_des = list(vcmax = c('f_scalar_none','f_tcor_des_modArrhenius','f_tcor_des_collatz1991','f_tcor_des_cox2001'))
)

pars.var      <- list(
  reftemp    = list(vcmax = 25),          # reference temperature at which Vcmax scalar = 1         (oC) 
  Ha         = list(vcmax = 54000),       # activation energy of Vcmax                              (J mol-1)
  Hd         = list(vcmax = 200000),      # deactivation energy of Vcmax                            (J mol-1)
  Topt       = list(vcmax = 32),          # temperature optimum of Vcmax                            (oC)
  deltaS     = list(vcmax  = numeric(1)), # 
  q10        = list(vcmax = 2),           # Q10 of Vcmax                                            (-)
  tupp_cox   = list(vcmax= 36),           # upper leaf T for Vcmax temp scaling from Cox 2001       (oC)
  tlow_cox   = list(vcmax= -20),          # lower leaf T for Vcmax temp scaling from Cox 2001       (oC)
  exp_cox    = list(vcmax = 0.3),         # exponent for Vcmax temp scaling from Cox 2001           (-)
  gstar_bf_a = 0.012,                     # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
  gstar_bf_b = 1.68,                      # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
  gstar_bf_c = 42.7                       # quadratic temperature dependence of gamma star from Brooks & Farquhar 1985 
)

pars_proc.var <- NA

pars_eval.var <- NA

env.var <- list(
  temp = seq(0,45,2)
)



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  leaf = list(
    fnames    = fnames.static,
    pars      = pars.static,
    env       = env.static
  ))

init_dynamic <- list(
  leaf = list(
    fnames    = fnames.var,
    pars      = pars.var,
    pars_proc = pars_proc.var,
    pars_eval = pars_eval.var,
    env       = env.var
  ))



### END ###