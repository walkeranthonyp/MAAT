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
    .$fns$gen_alpha      <- f_gen_alpha
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

  #.$fnames$transfer$t.1_to_2 <- 'f_eff12'
  #.$fnames$transfer$t.2_to_1 <- 'f_eff21'
  .$configure_test()

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
  #c1dec = 'f_c1dec',
  #c2dec = 'f_c2dec',
  decomposition = list(
    d.1 = 'f_c1dec',
    d.2 = 'f_c2dec'
  ),
  
  # transfer functions
  #eff12 = 'f_eff12',
  #eff21 = 'f_eff21'

  # transfer list
  transfer = list(
    t.1_to_2 = 'f_eff12',
    t.2_to_1 = 'f_eff21'
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
  c_pools    = matrix(c(c1 = .1, c2 = .1 ), ncol = 1),
  solver_out = matrix(1)
)

soil_decomp_object$alpha <- list()


# state parameters (i.e. calculated parameters)
####################################
soil_decomp_object$state_pars <- list(
  vcmax    = numeric(0)   
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
  c(unlist(.$state),unlist(.$statei_pars))
}



# test functions
#######################################################################        

soil_decomp_object$.test <- function(., verbose=F) {
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))

  .$run()
}



### END ###
