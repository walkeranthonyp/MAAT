################################
#
# MAAT gwater_rt model - initialisation script
# 
# AWalker (walkerap@ornl.gov) 
# January 2018
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
  recharge = 'f_rechrg_lin', 
  geology  = 'f_trans_single'
)

pars.static <- list(
  out_hsub = 10:14
) 

env.static  <- 1524 



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NA where variation is not required  

# define lists
fnames.var <- list(
  recharge = c('f_rechrg_lin','f_rechrg_power'),  # recharge function
  geology  = c('f_trans_single','f_trans_double') # geology/transport function
)

pars.var <- NULL

pars_proc.var <- list(
  a   =  'gwater_rt.recharge',
  b   =  'gwater_rt.recharge',
  K   =  'gwater_rt.geology',            
  K1  =  'gwater_rt.geology',            
  K2  =  'gwater_rt.geology'             
)

pars_eval.var <- list(
  a   =  "3.35 + rnorm(n)",
  b   =  "0.1  + 0.1*runif(n)",
  K   =  "15   + 1*rnorm(n)",            
  K1  =  "20   + 1*rnorm(n)",            
  K2  =  "10   + 1*rnorm(n)"            
)

env.var <- NULL



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  gwater_rt = list(
    fnames    = fnames.static,
    pars      = pars.static,
    env       = env.static
  ))

init_dynamic <- list(
  gwater_rt = list(
    fnames    = fnames.var,
    pars      = pars.var,
    pars_proc = pars_proc.var,
    pars_eval = pars_eval.var,
    env       = env.var
  ))



### END ###
