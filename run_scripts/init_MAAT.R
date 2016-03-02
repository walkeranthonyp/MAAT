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
  gstar       = 'f_temp_scalar_Arrhenius',
  Kc          = 'f_temp_scalar_Arrhenius',
  Ko          = 'f_temp_scalar_Arrhenius',
  vcmax       = 'f_vcmax_lin',
  jmax        = 'f_jmax_walker2014',
  vcmax_tcor  = 'f_temp_scalar_modArrhenius',
  jmax_tcor   = 'f_temp_scalar_modArrhenius',
  respiration = 'f_rd_collatz1991', 
  rs          = 'f_r_zero',
  rb          = 'f_r_zero'
)

leaf.pars.static <- NA

leaf.env.static  <- list(
  ca_conc = 400,
  par     = 500
)



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths

# define lists
leaf.fnames.var <- list(
  etrans = c('f_j_farquhar1980','f_j_farquharwong1984','f_j_harley1992','f_j_collatz1991'),
  ri     = c('f_r_zero','f_ri_constant')
)

leaf.pars.var <- list(
  avn_25 = rnorm(10,10,1)
)

leaf.env.var <- list(
  temp = c(5,20)
)



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_ls <- list(
  lfs = leaf.fnames.static,
  lps = leaf.pars.static,
  les = leaf.env.static,
  lfv = leaf.fnames.var,
  lpv = leaf.pars.var,
  lev = leaf.env.var
)


