################################
#
# mcmc_test object (SMO)
# 
# AWalker July 2018
#
################################

library(proto)
source('mcmc_test_functions.R')
source('mcmc_test_system_functions.R')



# TEMPLATE OBJECT
###############################################################################

mcmc_test_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'mcmc_test'
    
    # no expected child objects
    
    # build function
    build <- function(.) {
      as.proto(.$as.list())
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.) {
    
      # call system model 
      get(.$fnames$mcmc_testsys)(.)
      
      # print to screen
      if(.$cpars$verbose) {
        print(.$state)
      }
      
      # output
      .$output()
    } 
    
    
    
    ###########################################################################
    # Output functions

    #output processing function
    # -- returns a list of outputs
    output <- function(.) {
      if(.$fnames$mcmc_testsys == 'f_mcmc_testsys_mixture') {
        .$state$mixture_p
      } else if(.$fnames$mcmc_testsys == 'f_mcmc_testsys_regression') {
        .$state$linreg_y
      } else stop('Output type not defined')
    }    
    
    
    
    ###########################################################################
    # Variables etc
    
    # run control parameters
    cpars <- list(
      verbose       = F,      # write diagnostic output during runtime 
      output        = 'full'  # type of output from run function
    )
    
    # function names
    fnames <- list(
      mcmc_testsys = 'f_mcmc_testsys_mixture',
      reg_func     = 'f_reg_func_linear'
    )

    #parameters
    pars   <- list(
      # mixture model parameters
      mixture_scale = 1e12, 
      height1       = 0.5, 
      height2       = 0.2, 
      height3       = 0.3, 
      mu1           = -5, 
      mu2           = 0, 
      mu3           = 5, 
      sd1           = 2,
      sd2           = 2,
      sd3           = 2,
      proposal1     = 1,
      proposal2     = 1,
      proposal3     = 1,
      proposal4     = 1,
      
      # linear regresssion model parameters
      syn_a_mu = 2,
      syn_b_mu = 7,          
      syn_a_sd = 3,          
      syn_b_sd = 2,
      a        = 7,
      b        = 2          
    )

    # environment
    env <- list(
      dummy    = 1,
      linreg_x = numeric(1)    
    )

    # state parameters (i.e. calculated parameters)
    state_pars <- list(
      none      = numeric(0)   
    )
    
    # state
    state <- list(
      mixture_p = numeric(1), # probability output from mixture model
      linreg_y  = numeric(1)  # y value output from line function
    )
    
    
    ###########################################################################
    # Run & configure functions
    
    configure <- function(.,vlist,df,o=T) {
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
      
      # name and assign the UQ variables
      uqvars <- names(df)
      prefix <- vapply( strsplit(uqvars,'.', fixed=T), function(cv) cv[1], 'character' )
      modobj <- .$name
      dfss   <- which(prefix==modobj)
      vlss   <- match(uqvars[dfss], paste0(modobj,'.',names(.[[vlist]])) )
      
      # catch NAs in vlss
      if(any(is.na(vlss))) stop(paste('names mismatch between model object variables and input list variable:', uqvars[which(is.na(vlss))] ))

      # assign UQ variables
      .[[vlist]][vlss] <- df[dfss]
      
      if(.$cpars$verbose&o) {
        print('',quote=F)
        print('mcmc_test configure:',quote=F)
        print(prefix,quote=F)
        print(df,quote=F)
        print(.[vlist],quote=F)
      }
    }
    
    run_met <- function(.,l) {
      # This wrapper function is called from an lapply function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are sequential
      # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
      # typically used to run the model using data collected at a specific site and to compare against observations
      
      # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
      # any "env" variables specified in the "dataf$env" dataframe but also specified in .$dataf$met will be overwritten by the .$dataf$met values 
      
      # met data assignment
      # .$configure(vlist='env',df=.$dataf$met[l,],F)
      
      # run model
      .$run()              
    }
    
    
    
    #######################################################################        
    # Test functions
    # - not copied when the object is cloned

    .test_mixture <- function(.,verbose=T,verbose_loop=T) {
      
      if(verbose) {
        str(.)
      }

      #.$cpars$verbose       <- verbose
      #.$cpars$verbose_loop  <- verbose_loop
      #.$cpars$output        <-'full'
      
     .$fnames$mcmc_testsys <- 'f_mcmc_testsys_mixture' 

      .$run()
    }

    .test_linreg <- function(.,verbose=T,verbose_loop=T) {
                  
      
      if(verbose) {
        str(.)
      }

      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
     .$fnames$mcmc_testsys <- 'f_mcmc_testsys_regression' 
      
      .$run()
    }


      
#######################################################################        
# end object      
})



### END ###
