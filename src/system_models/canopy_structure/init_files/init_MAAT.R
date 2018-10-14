################################
#
# MAAT Canopy Structure Model - example initialisation script
# 
# AWalker (walkerap@ornl.gov) 
# October 2018
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script initialises the leaf and/or canopy photosynthesis version of MAAT
# - by setting the values of the init_static and init_dynamic lists
# - each of these lists elements are sublists named with the names of the model objects to which the variables in the sublists belong  
# - each of these model object sublist elements are sublists named with the type of variable they comprise - fnames, pars, or env  
# - each of these variables sublist elements are sublists named with the variable to which they refer in the model object  

# To setup a model object simulation,
#  set up the following lists, for OBJ leaf and canopy:

#  set static variables during runtime (init_static) by setting:
#  - OBJ.fnames.static
#  - OBJ.pars.static     
#  - OBJ.env.static

#  set dynamic variables during runtime (init_dynamic) by setting:
#  - OBJ.fnames.var
#  - OBJ.pars.var     
#  - OBJ.env.var

#  for a process SA two more sublists need to be added as elements to the init_dynamic lists:
#   OBJ.pars_proc.var and OBJ.pars_eval.var, both with the same named elements as OBJ.pars.var
#   - OBJ.pars_proc.var is a list of strings that associate each parameter with the process it is associated with, these should be one of the element names in the list init_dynamic$OBJ$fnames
#   - OBJ.pars_eval.var is a string that once evaluated gives a vector of parameter samples  


# This script must set the values of the above 6(8) lists per model object, 
# if you do not want to set any values using these lists they must be specified as NULL

###################################################################



### Define the functions, parameters, and environment variables that are to remain static through a simulation 
###############################

# define lists
leaf.fnames.static <- list(
  solver_func = 'f_A_r_leaf',
  gstar       = 'f_scalar_none',
  ri          = 'f_r_zero',
  rs          = 'f_r_zero',
  rb          = 'f_r_zero'
  )

leaf.pars.static <- NULL 

leaf.env.static  <- list(
  par     = 1000,
  temp    = 25
  )


canopy.fnames.static <- NULL   

canopy.pars.static   <- NULL

canopy.env.static    <- NULL


canopy_structure.fnames.static <- NULL   

canopy_structure.pars.static   <- NULL

canopy_structure.env.static    <- NULL



### Define the functions, parameters, and environment variables that are to be varied 
###############################
# if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length - need to put a check for this in the wrapper 
# if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
# set to NULL where variation is not required  

# define lists
leaf.fnames.var    <- NULL 

leaf.pars.var      <- NULL

leaf.pars_proc.var <- NULL

leaf.pars_eval.var <- NULL

leaf.env.var       <- list(
  ca_conc = seq(50,1500,50)
)


canopy.fnames.var    <- NULL

canopy.pars.var      <- NULL 

canopy.pars_proc.var <- NULL

canopy.pars_eval.var <- NULL

canopy.env.var       <- NULL


canopy_structure.fnames.var    <- NULL

canopy_structure.pars.var      <- NULL 

canopy_structure.pars_proc.var <- NULL

canopy_structure.pars_eval.var <- NULL

canopy_structure.env.var       <- NULL



### Combine the functions, parameters, and environment static and variable lists into a single list 
###############################

init_static <- list(
  leaf = list(
    fnames = leaf.fnames.static,
    pars   = leaf.pars.static,
    env    = leaf.env.static
  ),
  canopy = list(
    fnames = canopy.fnames.static,
    pars   = canopy.pars.static,
    env    = canopy.env.static
  ),
  canopy_structure = list(
    fnames = canopy_structure.fnames.static,
    pars   = canopy_structure.pars.static,
    env    = canopy_structure.env.static
  )
)

init_dynamic <- list(
  leaf = list(
    fnames    = leaf.fnames.var,
    pars      = leaf.pars.var,
    pars_proc = leaf.pars_proc.var,
    pars_eval = leaf.pars_eval.var,
    env       = leaf.env.var
  ),
  canopy = list(
    fnames    = canopy.fnames.var,
    pars      = canopy.pars.var,
    pars_proc = canopy.pars_proc.var,
    pars_eval = canopy.pars_eval.var,
    env       = canopy.env.var
  ),
  canopy_structure = list(
    fnames    = canopy_structure.fnames.var,
    pars      = canopy_structure.pars.var,
    pars_proc = canopy_structure.pars_proc.var,
    pars_eval = canopy_structure.pars_eval.var,
    env       = canopy_structure.env.var
  )
)



### END ####
