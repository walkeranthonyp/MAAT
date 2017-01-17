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
# - by setting the values of these lists:

# Static Variables:

# fnames.static
# pars.static     
# env.static


# Dynamic Variables:

# fnames.var
# pars.var     
# env.var

# This script must set the values of the above 6 lists, even if their value is NA

###################################################################



### Define the functions, parameters, and environment variables that are to remain static through a simulation 
###############################

# define lists
leaf.fnames.static <- list(
  solver_func = 'f_A_r_leaf',
  gstar       = 'f_temp_scalar_no_response',
  Kc          = 'f_temp_scalar_no_response',
  Ko          = 'f_temp_scalar_no_response',
  vcmax       = 'f_constant_vcmax',
  jmax        = 'f_jmax_walker2014',
  vcmax_tcor  = 'f_temp_scalar_no_response',
  jmax_tcor   = 'f_temp_scalar_no_response',
  respiration = 'f_rd_collatz1991', 
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

# define lists
leaf.fnames.var <- list(
  etrans = c('f_j_farquharwong1984','f_j_harley1992','f_j_collatz1991')
)

leaf.pars.var <- list(
  e_ajv_25 = rnorm(n,1,0.02)
)

leaf.env.var <- list(
  ca_conc = c(280,400,600)
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
    fnames = leaf.fnames.var,
    pars   = leaf.pars.var,
    env    = leaf.env.var
  ))

# for an process SA two more elements need to be added to the init_dynamic$leaf list
# - pars_proc and pars_eval, both with the same named elements
# - pars_proc is a list of strings that associate each parameter with the process it is associated with, these should be one of the elemnet names in the list init_dynamic$X$fnames
# - pars_eval is a string that once evaluated gives a vector of parameter samples  
