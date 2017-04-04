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
# - leaf.fnames.static
# - leaf.pars.static     
# - leaf.env.static


# set dynamic variables during runtime (init_dynamic) by setting:
# - leaf.fnames.var
# - leaf.pars.var     
# - leaf.env.var
# for a process SA two more sublists need to be added as elements to the init_dynamic$leaf list:
# - pars_proc and pars_eval, both with the same named elements
#   - pars_proc is a list of strings that associate each parameter with the process it is associated with, these should be one of the elemnet names in the list init_dynamic$X$fnames
#   - pars_eval is a string that once evaluated gives a vector of parameter samples  


# This script must set the values of the above 6 lists, even if their value is NA

###################################################################



### Define the functions, parameters, and environment variables that are to remain static through a simulation 
###############################

# define lists
leaf.fnames.static <- list(
  solver_func = 'f_A_r_leaf',
  gstar       = 'f_scalar_none',
  Kc_tcor     = 'f_scalar_none',
  Ko_tcor     = 'f_scalar_none',
  vcmax       = 'f_constant_vcmax',
  jmax        = 'f_jmax_walker2014',
  vcmax_tcor  = 'f_scalar_none',
  jmax_tcor   = 'f_scalar_none',
  respiration = 'f_rd_lin_vcmax', 
  ri          = 'f_r_zero',
  rs          = 'f_r_zero',
  rb          = 'f_r_zero'
  )

leaf.pars.static <- list(
  atref.vcmax  = 50
  )

leaf.env.static  <- list(
  par     = 1000,
  temp    = 25
  )



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NA where variation is not required  

# define lists
leaf.fnames.var <- list(
  etrans = c('f_j_farquharwong1984','f_j_harley1992','f_j_collatz1991')
)

leaf.pars.var <- NA

leaf.pars_proc.var <- NA

leaf.pars_eval.var <- NA

leaf.env.var <- list(
  ca_conc = seq(50,1500,50)
)



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  leaf = list(
    fnames = leaf.fnames.static,
    pars   = leaf.pars.static,
    env    = leaf.env.static
  ))

init_dynamic <- list(
  leaf = list(
    fnames    = leaf.fnames.var,
    pars      = leaf.pars.var,
    pars_proc = leaf.pars_proc.var,
    pars_eval = leaf.pars_eval.var,
    env       = leaf.env.var
  ))

