################################
#
# MAAT soil_decomp model object 
#
# AWalker, Matt Craig, October 2019 
#
################################

library(proto)
source('soil_decomp_functions.R')
source('soil_decomp_system_functions.R')



# soil_decomp OBJECT
###############################################################################

# use generic soil_decomp
setwd('..')
source('generic_model_object.R')
soil_decomp_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('soil_decomp')



# assign object functions
###########################################################################
soil_decomp_object$name <- 'soil_decomp'
# if this new object will contain nested child objects uncomment and edit the below lines of code and delete this comment 
# - otherwise delete all 5 lines
#soil_decomp_object$child_list      <- list('child_name1') 
#soil_decomp_object$build_child     <- build_child  
#soil_decomp_object$configure_child <- configure_child  



# function to configure unique elements of the object
# - adds functions to fns that are not in fnames
# - or functions that are derivations of other functions, in this case teh rs derived fuinctions like rs_r0 and rs_fe 
####################################
soil_decomp_object$configure_unique <- function(., init=F, flist=NULL ) {
  if(init) {
    source('../../functions/packagemod_functions.R')
    .$fns$plsoda         <- plsoda
    .$fns$inputrates     <- f_inputrates
    .$fns$DotO           <- f_DotO
    .$fns$DotC           <- f_DotC
    .$fns$transfermatrix <- f_transfermatrix
    .$fns$solver_func    <- f_solver_func
  }

  if(any(names(flist)=='xyz')) {
   .$fns$xyz_fe <- get(paste0(.$fnames$xyz,'_fe'), pos=1 )
  }
}


# assign unique run function
###########################################################################
# if run function needs to be modified - add new function here
soil_decomp_object$run <- function(.) {

  # call system model
  .$fns$sys()

  # print to screen
  if(.$cpars$verbose) print(.$state)

  # output
  .$output()
}


# functions unique to object that do not live in fnames/fns, i.e. do not vary ever
###########################################################################
# add structural functions (i.e not alternative process functions)  unique to model object here

# assign object variables 
###########################################################################

# function names
####################################
soil_decomp_object$fnames <- list(
  sys   = 'f_sys_m2pool',
  
  # decay/decomposition functions
  decomp = list(
    d1 = 'f_decomp_MM_microbe',
    d2 = 'f_decomp_lin'
  ),
  
  # transfer list
  transfer = list(
    t1_to_2 = 'f_transfer_resploss',
    t2_to_1 = 'f_transfer_all'
  )
)


# environment
####################################
soil_decomp_object$env <- list(
  litter = 3.2,
  times  = 2 
)


# state
####################################
soil_decomp_object$state <- list(
  c_pools    = matrix(c(c1=0.1, c2=0.1 ), ncol=1 )
)


# state parameters (i.e. calculated parameters)
####################################
soil_decomp_object$state_pars <- list(
  solver_out = matrix(1)
)


# parameters
####################################
soil_decomp_object$pars <- list(
  ks = 1.8e-05,
  kb = 0.007,
  Km = 900,
  r  = 0.6,
  Af = 1
)


# run control parameters
####################################
soil_decomp_object$cpars <- list(
  verbose  = F,          # write diagnostic output during runtime 
  cverbose = F,          # write diagnostic output from configure function 
  output   = 'run'       # type of output from run function
)



# output functions
#######################################################################        

f_output_soil_decomp_run <- function(.) {
  unlist(.$state)
}

f_output_soil_decomp_state <- function(.) {
  unlist(.$state)
}

f_output_soil_decomp_full <- function(.) {
  c(unlist(.$state),unlist(.$state_pars))
}



# test functions
#######################################################################        

soil_decomp_object$.test <- function(., verbose=F ) {
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  .$configure_test()

  .$run()
}



### END ###
