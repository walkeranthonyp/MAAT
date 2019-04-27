################################
#
# MAAT gwater model - initialisation script
# 
# AWalker (walkerap@ornl.gov) 
# November 2017
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script initialises the gwater MAAT system model 
# by setting the values of the init_static and init_dynamic lists
# - each of these lists elements are sublists named with the names of the model objects to which the variables in the sublists belong  
# - each of these model object sublist elements are sublists named with the type of variable they comprise - fnames, pars, or env  
# - each of these variables sublist elements are sublists named with the variable to which they refer in the model object  

# for example, to setup a gwater model simulation:
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
fnames.static <- list()

pars.static <- list()

env.static  <- list()



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NA where variation is not required  

# define lists
fnames.var <- list()

pars.var <- NA

pars_proc.var <- NA

pars_eval.var <- NA

env.var <- list()



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  fnames = list( gwater_rt = fnames.static),
  pars   = list( gwater_rt = pars.static),
  env    = list( gwater_rt = env.static)
)

init_dynamic <- list(
  fnames    = list( gwater_rt = fnames.var),
  pars      = list( gwater_rt = pars.var),
  pars_proc = list( gwater_rt = pars_proc.var),
  pars_eval = list( gwater_rt = pars_eval.var),
  env       = list( gwater_rt = env.var)
)



### END ###
