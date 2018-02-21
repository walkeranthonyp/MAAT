################################
#
# gwater_rt for MAAT object functions
# from Dai etal 2017 WRR 
#
# AWalker November 2017
#
################################

library(proto)
library(stringr)

source('gwater_rt_functions.R')



# SIMPLE GROUND WATER REACTIVE TRANSPORT OBJECT
###############################################################################

gwater_rt_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'gwater_rt'
    
    # no expected child objects
    
    # build function
    build <- function(.) {
      as.proto(.$as.list())
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.) {
      
      # call system model 
      get(.$fnames$gwatersys)(.)

      # print to screen
      if(.$cpars$verbose) {
        print(.$state)
      }
      
      # output
      .$output()
    } 
    
    
    
    ###########################################################################
    # Output functions

    # output processing function
    # -- returns a vector of outputs
    output <- function(.) {
      if(.$cpars$output=='run') {
        lout <- .$state$h
      } else if(.$cpars$output=='slim') {
        lout <- c(recharge=.$state$recharge,h=.$state$h[.$pars$out_hsub])
      } else stop()
      
      lout
    }    
    
    
    
    ###########################################################################
    # Variables etc
    
    # run control parameters
    cpars <- list(
      verbose  = F,               # write diagnostic output during runtime 
      output   = 'run'            # type of output from run function
    )
    
    # function names
    fnames <- list(
      gwatersys = 'f_gwatersys_daiye', # system model
      rechrg    = 'f_rechrg_lin',     # recharge function
      geol      = 'f_trans_single'    # geology/transport function
    )

    # gwater_rt parameters
    pars   <- list(
      nx  = 21,                   # horizontal discretisation of geological domain   
      K   = 15,                   # hydraulic conductivity for single domain geology 
      K1  = 20,                   # hydraulic conductivity for double domain geology           
      K2  = 10,                   # hydraulic conductivity for double domain geology         
      a   = 3.35,                 # recharge scalar power law (unitless)
      b   = 0.15,                 # recharge scalar linear    (unitless)
      L   = 1e4,                  # horizontal domain length (m)
      out_hsub = 13               # subscript of .$state$h vector for 'slim' output 
    )

    # gwater_rt environment
    env <- list(
      h1     = 180,               # hydraulic head left  (x=0) boundary condition (m)
      h2     = 100,               # hydraulic head right (x=L) boundary condition (m)
      precip = 1524               # precipitation                                 (mm)
    )

    # gwater_rt state parameters (i.e. calculated parameters)
    state_pars <- list(
      x        = numeric(0)
    )
    
    # gwater_rt state
    state <- list(
      recharge = numeric(1),      # recharge rate                    (m d-1?)
      h        = numeric(21)      # hydraulic head across the domain (Pa)
    )
    
    
    
    ###########################################################################
    # Run & configure functions
    
    configure <- function(.,vlist,df,o=T){
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
      
      # name and assign the UQ variables
      uqvars <- names(df)
      prefix <- substr(uqvars,1,str_locate(uqvars,'\\.')[,2]-1)
      modobj <- .$name
      dfss   <- which(prefix==modobj)
      vlss   <- match(uqvars[dfss], paste0(modobj,'.',names(.[[vlist]])) )
      # could write a line to catch NAs in vlss
      .[[vlist]][vlss] <- df[dfss]
      
      # if(.$cpars$cverbose&o) {
      #   print('',quote=F)
      #   print('gwater_rt configure:',quote=F)
      #   print(prefix,quote=F)
      #   print(df,quote=F)
      #   print(.$fnames,quote=F)
      # }
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

    .test <- function(., verbose=F, gwater_rt.precip=1524, gwater_rt.recharge='f_rechrg_lin', gwater_rt.geology='f_trans_single' ){
      
      if(verbose) {
        str(.)
        print(.$env)
      }
      
      .$cpars$verbose   <- verbose
      .$cpars$output    <-'run'
      
      .$env$precip      <- gwater_rt.precip
      .$fnames$recharge <- gwater_rt.recharge
      .$fnames$geology  <- gwater_rt.geology
      
      .$run()
    }

    
    
    #######################################################################        
    # end object      
})





