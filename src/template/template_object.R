################################
#
# template for MAAT object functions
# 
# AWalker March 2017
#
################################

library(proto)
library(stringr)

source('template_functions.R')



# LEAF OBJECT
###############################################################################

template_object <- 
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'template'
    
    # no expected child objects
    
    # build function
    .build <- function(.) {
      as.proto(.$as.list())
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.) {
      

      # calculate state parameters
      # photosynthetic parameters
      .$state_pars$vcmax   <- get(.$fnames$vcmax)(.)
      .$state_pars$jmax    <- get(.$fnames$jmax)(.)
      .$state_pars$tpu     <- get(.$fnames$tpu)(.)
      .$state_pars$rd      <- get(.$fnames$respiration)(.)
      .$state_pars$alpha   <- 0.5 * (1-.$pars$f)

      .$state$respiration  <- get(.$fnames$rl_rd_scalar)(.) * .$state$respiration
      # determine rate limiting step - this is done based on carboxylation, not net assimilation (Gu etal 2010).
      .$state$A       <- get(.$fnames$solver)(.)      
      # assign the limitation state - assumes the minimum is the dominant limiting rate
      .$state$lim     <- c('wc','wj','wp')[which(c(.$state$wc,.$state$wj,.$state$wp)==min(c(.$state$wc,.$state$wj,.$state$wp),na.rm=T))]       
      
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
      if(.$cpars$output=='slim') {
        lout <- 
          list(A=.$state$A,ci=.$state$ci,lim=.$state$lim)
        
      } else if(.$cpars$output=='run') {
        lout <- 
          list(A=.$state$A,cc=.$state$cc,ci=.$state$ci,
               gi=1/.$state_pars$ri,gs=1/.$state_pars$rs,gb=1/.$state_pars$rb,
               respiration=.$state$respiration,lim=.$state$lim)
        
      } else stop()
      
    }    
    
    
    
    ###########################################################################
    # Variables etc
    
    # run control parameters
    cpars <- list(
      verbose       = F,          # write diagnostic output during runtime 
      output        = 'run'       # type of output from run function
    )
    
    # function names
    fnames <- list(
      gstar               = 'f_gstar_constref',
      Alim                = 'f_lim_farquhar1980'
    )

    #template parameters
    pars   <- list(
      a             = 0.80,       # fraction of PAR absorbed                               (unitless)  --- this should equal 1 - template scattering coefficient, there is potential here for improper combination of models
      #physical constants
      R   = 8.31446               # molar gas constant                                      (m2 kg s-2 K-1 mol-1  ==  Pa m3 mol-1K-1)
    )

    # template environment
    env <- list(
      ca_conc   = numeric(0),          # (umol mol-1)
      atm_press = 101325               # ( Pa)
      )

    #template state parameters (i.e. calculated parameters)
    state_pars <- list(
      vcmax    = numeric(0),   # umol m-2 s-1
      cica_chi = numeric(0)    # Ci:Ca ratio 
    )
    
    # template state
    state <- list(
      #environmental state
      oi = numeric(0),                 # atmospheric & internal O2  (kPa)
      transition   = numeric(0)        # cc at the transition point where wc = wj                        (Pa)
    )
    
    
    
    ###########################################################################
    # Run & configure functions
    
    configure <- function(.,func,df,o=T){
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
      
      # name and assign the UQ variables
      uqvars     <- names(df)
      prefix     <- substr(uqvars,1,str_locate(uqvars,'\\.')[,2]-1)
      lapply(uqvars[which(prefix=='template')], func, .=., df=df)

      if(.$cpars$cverbose&o) {
        print('',quote=F)
        print('template configure:',quote=F)
        print(prefix,quote=F)
        print(df,quote=F)
        print(.$fnames,quote=F)
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
      .$configure(func='write_env',df=.$dataf$met[l,],F)
      
      # run model
      .$run()              
    }
    
    
    
    #######################################################################        
    # Test functions
    # - not copied when the object is cloned

    .test_template <- function(.,verbose=T,verbose_loop=T,template.par=1000,template.ca_conc=300,rs='f_rs_medlyn2011',gd='f_ficks_ci'){
      
      if(verbose) {
        str.proto(.)
        print(.$env)
      }
      .$cpars$verbose       <- verbose
      .$cpars$verbose_loop  <- verbose_loop
      .$cpars$output        <-'full'
      
      .$fnames$ri          <- 'f_r_zero'
      .$fnames$rs          <- rs
      .$fnames$solver_func <- 'f_A_r_template'
      .$fnames$gas_diff    <- gd
      
      .$env$par     <- template.par
      .$env$ca_conc <- template.ca_conc
      
      .$run()
    }

    #######################################################################        
    # end object      
})





