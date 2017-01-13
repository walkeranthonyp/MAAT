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
# parsB.var - for process UQ     
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

leaf.parsB.var <- list(
  e_bjv_25 = rnorm(3*n^2,0.89,0.05)
)

leaf.env.var <- list(
  ca_conc = c(280,400,600)
)



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_ls <- list(
  lfs  = leaf.fnames.static,
  lps  = leaf.pars.static,
  les  = leaf.env.static,
  lfv  = leaf.fnames.var,
  lpv  = leaf.pars.var,
  lpBv = leaf.parsB.var,
  lev  = leaf.env.var
)


