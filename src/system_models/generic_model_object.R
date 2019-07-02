################################
#
# Generic system model object 
# 
# AWalker Dec 2018
#
################################

library(proto)

source('generic_model_functions.R')



# GENERIC SYSTEM MODEL OBJECT
###############################################################################

system_model_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- character(1)
    
    # child objects 
    child_list <- NULL 
    
    # build function 
    build <- function(.) {}
    
    # build function 
    build_child <- NULL 
    
    
    ###########################################################################
    # main run function
    
    run <- function(.) {} 
    
    
    ###########################################################################
    # Output functions

    # -- returns a vector of outputs
    state_retrive <- function(.,snames,state='state') {
      lsubs <- match(snames,names(.[[state]]))
      unlist(.[[state]][lsubs])
    }
    
    output <- function(.) {}    

    
    ###########################################################################
    # Variables etc
    
    # function names
    fnames     <- list()
    
    # functions - created as a child proto object during first call to configure  
    #fns        <- proto(.=system_model_object)
    
    # environment
    env        <- list()

    # state
    state      <- list()
    
    # state parameters
    state_pars <- list()
    
    # parameters
    pars       <- list()
    
    # run control parameters
    cpars      <- list()
    
    
    ###########################################################################
    # configure & run_met functions

    # function that gets passed named values from MAAT wrapper object 
    # - and assigns them to the model object  
    configure <- function(.) {}
 
    # configure a list variable 
    configure_sublist <- function(.) {}

    # configure functions not found in fnames list
    # - either unique functions that never change or 
    # - sub-functions associated with functions that are in fnames  
    configure_unique <- function(.) {NULL}

    # run model over a meteorological (or boundary condition) dataset
    run_met <- function(.) {}
    
    
    ###########################################################################
    # Test functions
    # - not copied when the object is cloned 
    # - (signaled by the . at the beginning of the functions name)

    .test <- function(.) {}
    
    
#######################################################################        
# end object      
})



### END ###
