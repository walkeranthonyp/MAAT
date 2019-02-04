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
# if you do not want to set any values using these lists they must be specified as NULL

###################################################################



### Define the functions, parameters, and environment variables that are to remain static through a simulation 
###############################

# define lists
fnames.static <- list(
  solver_func     = 'f_A_r_leaf',
  vcmax           = 'f_vcmax_constant',
  jmax            = 'f_jmax_power',
  tcor_asc = list(  
    vcmax   = 'f_scalar_none',
    jmax    = 'f_scalar_none',
    Kc      = 'f_scalar_none',
    Ko      = 'f_scalar_none'
  ),
  tcor_des = list(  
    vcmax   = 'f_scalar_none',
    jmax    = 'f_scalar_none'
  ),
  rd     = 'f_rd_lin_vcmax', 
  ri     = 'f_r_zero',
  rb     = 'f_r_zero'
  )

pars.static <- list(
  atref = list(vcmax = 50 )
  )

env.static  <- list(
  par     = 1000,
  temp    = 25
  )



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NULL where variation is not required  

# define lists
fnames.var    <- list(
  etrans = c('f_j_farquharwong1984','f_j_harley1992','f_j_collatz1991')
)

pars.var      <- NULL

pars_proc.var <- NULL

pars_eval.var <- NULL

env.var       <- list(
  ca_conc = seq(50,1500,50)
)



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  fnames = list( leaf = fnames.static),
  pars   = list( leaf = pars.static),
  env    = list( leaf = env.static)
)

init_dynamic <- list(
  fnames    = list( leaf = fnames.var),
  pars      = list( leaf = pars.var),
  pars_proc = list( leaf = pars_proc.var),
  pars_eval = list( leaf = pars_eval.var),
  env       = list( leaf = env.var)
)



### END ####
