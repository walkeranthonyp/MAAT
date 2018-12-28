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
    fnames <- list()
    
    # leaf environment
    env <- list()

    # leaf state
    state <- list()
    
    #leaf parameters
    pars   <- list()
    
    # run control parameters
    cpars <- list()
    
    
    
    ###########################################################################
    # configure & run_met functions

    configure <- function(.) {}
 
    
    # configure a list variable 
    configure_sublist <- function(.) {}

    
    run_met <- function(.) {}
    
    
    
    #######################################################################        
    # Test functions
    # - not copied when the object is cloned

    .test <- function(.) {}
    
#######################################################################        
# end object      
})



### END ###
