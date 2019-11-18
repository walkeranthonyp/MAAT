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
    child_list  <- NULL 
    
    # build function 
    build       <- build
    
    # build function 
    build_child <- NULL 
    
    
    ###########################################################################
    # main run functions
    
    run  <- run 
    
    
    ###########################################################################
    # Output functions

    # -- returns a vector of outputs
    state_retrive <- function(.,snames,state='state') {
      lsubs <- match(snames,names(.[[state]]))
      unlist(.[[state]][lsubs])
    }
   
    # output is a function assigned during build function execution 
    output <- function(.) {}    

    
    ###########################################################################
    # Variable lists
   
    # calculated by model 
    # - 2 lists
    # state
    state      <- list()
    
    # state parameters
    state_pars <- list()
    
    # lists with values configurable by user
    # - 3 lists
    # - within each list, all elements must be of the same class  
    # function names
    fnames     <- list()
    
    # environment
    env        <- list()

    # parameters
    pars       <- list()
   
    # functions 
    # - place holder for fns proto object that lives within model object
    # - created as a child proto object during first call to configure  
    # - this object has methods named exactly as the elements in fnames are named
    # - each of these methods are the functions named in the elements of the fnames list 
    #fns        <- proto(.=system_model_object)
    
    # run control parameters
    cpars      <- list(
      output   = 'run', # output function
      verbose  = F,     # standard verbose output for debugging
      cverbose = F,     # more in-depth verbose output from configure function
      diag     = F      # model specific diagnostic output - written to state list
    )
    
    
    ###########################################################################
    # configure & run_met functions

    # function that gets passed named values from MAAT wrapper object 
    # - and assigns them to the model object  
    configure         <- configure 
 
    # configure a list variable 
    configure_sublist <- configure_sublist 

    # configure child objects 
    configure_child   <- NULL  

    # configure functions not found in fnames list
    # - either unique functions that never change or 
    # - sub-functions associated with functions that are in fnames  
    configure_unique  <- NULL 

    # configure for test functions 
    configure_test    <- configure_test 

    # run model over a meteorological (or boundary condition) dataset
    init              <- NULL
    run_met           <- run_met 
    run_met1          <- run_met1 
    
    # configure env from met dataset 
    configure_met     <- configure_met 

    # configure check 
    configure_check   <- configure_check 

    
    ###########################################################################
    # Test functions
    # - not copied when the object is cloned 
    # - (signaled by the . at the beginning of the functions name)

    .test <- function(.) {}
    
    
#######################################################################        
# end object      
})



### END ###
