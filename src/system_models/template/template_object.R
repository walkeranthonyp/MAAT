################################
#
# template for MAAT object functions
# 
# AWalker March 2017
#
################################

library(proto)
source('template_functions.R')
source('template_system_functions.R')



# TEMPLATE OBJECT
###############################################################################

template_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'template'
    
    # no expected child objects
    
    # build function
    build <- function(.) {
      as.proto(.$as.list())
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.) {
      
      # call system model 
      get(.$fnames$templatesys)(.)
      
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
      if(.$cpars$output=='full') {
        print('done')
      else if(.$cpars$output=='none') { 
        print('')
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
      templatesys = 'f_templatesys_1',
      text        = 'f_text_combine',
      calcval     = 'f_calcval_product',
      print       = 'f_print_textonly'
    )

    #template parameters
    pars   <- list(
      text1 = 'hello',    
      text2 = 'world',    
      text3 = 'The answer is:',    
      val1  = 6,           
      val2  = 7           
    )

    # template environment
    env <- list(
      ca_conc   = numeric(0)    
      )

    #template state parameters (i.e. calculated parameters)
    state_pars <- list(
      vcmax    = numeric(0)   
    )
    
    # template state
    state <- list(
      text    = character(0),     # text to print
      calcval = numeric(0)        # value to print
    )
    
    
    
    ###########################################################################
    # Run & configure functions
    
    configure <- function(.,vlist,df,o=T){
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
      
      if(.$cpars$cverbose&o) {
        print('',quote=F)
        print('template configure:',quote=F)
        print(prefix,quote=F)
        print(df,quote=F)
        print(.[vlist],quote=F)
      }
    }
    
    run_met <- function(.,l){
      # This wrapper function is called from an lapply function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are sequential
      # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
      # typically used to run the model using data collected at a specific site and to compare against observations
      
      # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
      # any "env" variables specified in the "dataf$env" dataframe but also specified in .$dataf$met will be overwritten by the .$dataf$met values 
      
      # met data assignment
      .$configure(vlist='env',df=.$dataf$met[l,],F)
      
      # run model
      .$run()              
    }
    
    
    
    #######################################################################        
    # Test functions
    # - not copied when the object is cloned

    .test <- function(.,verbose=T,verbose_loop=T) {
      
      if(verbose) {
        str.proto(.)
      }

      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
      .$run()
    }

    .test_change_func <- function(.,verbose=T,verbose_loop=T,
                                  template.text='f_text_combine',template.calcval='f_calcval_product',template.print='f_print_textonly') {
      
      if(verbose) {
        str.proto(.)
      }

      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
      .$fnames$text         <- template.text
      .$fnames$calcval      <- template.calcval
      .$fnames$print        <- template.print
      
      .$run()
    }

    .test_change_pars <- function(.,verbose=T,verbose_loop=T,
                                  template.text1='hello',template.text2='world') {
      
      if(verbose) {
        str.proto(.)
      }

      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
      .$pars$text1          <- template.text1
      .$pars$text2          <- template.text2
      
      .$run()
    }
      
      
    #######################################################################        
    # end object      
})



### END ###
