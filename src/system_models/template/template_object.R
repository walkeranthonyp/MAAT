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
    
    # child objects - none expected 
    child_list <- NULL 
    
    # build function 
    build <- function(., mod_mimic=NULL, ... ) {
      # no expected child objects
      ## if template model requires child objects, see build function in canopy_object.R ## 

      # read default model setup for highest level model
      init_default <- readXML(paste(.$name,'default.xml',sep='_'))
     
      # read model mimic setup
      if(!is.null(mod_mimic)) {
        setwd('mimic_xmls')
        print(paste('template mimic:', mod_mimic ))
        init_mimic   <- readXML(paste(.$name,'_',mod_mimic,'.xml',sep=''))
        init_default <- fuselists(init_default,init_mimic)
        setwd('..')
      }

      init_default
    }
    
    
    
    ###########################################################################
    # main run function
    
    run   <- function(.) {
     
      print('', quote=F )
      print('Run Model:', quote=F )
      print('', quote=F )
 
      # call system model 
      get(.$fnames$templatesys)(.)
      
      # output
      .$output()
    } 
    
    
    
    ###########################################################################
    # Output functions

    #output processing function
    # -- returns a list of outputs
    output <- function(.) {
      if(.$cpars$output=='full') {
        ## delete the two below print statements to allow a single output vector ##  
        print('', quote=F )
        print('Output:', quote=F )
        unlist(.$state)
      } else if(.$cpars$output=='none') { 
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
    ## - this is a place holder for model inputs, not used in this example
    env <- list(
      ca_conc   = numeric(0)    
      )

    #template state parameters (i.e. calculated parameters)
    ## - this is a place holder for state variables taht are also parameters, not used in this example
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
    
    configure <- function(., vlist, df, o=T ) {
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames / .$pars / .$env / .$state to the values passed in df 

      # split variable names at . 
      listnames <- vapply( strsplit(names(df),'.', fixed=T), function(cv) {cv3<-character(3); cv3[1:length(cv)]<-cv; t(cv3)}, character(3) )

      modobj <- .$name
      # df subscripts for model object
      moss   <- which(listnames[1,]==modobj)
      # df subscripts for model object sublist variables (slmoss) and model object numeric variables (vlmoss) 
      slss   <- which(listnames[3,moss]!='') 
      if(length(slss)>0) {
        slmoss <- moss[slss] 
        vlmoss <- moss[-slss] 
      } else {
        slmoss <- NULL 
        vlmoss <- moss 
      }
      # variable list subscripts for non-sublist variable variables 
      vlss   <- match(listnames[2,vlmoss], names(.[[vlist]]) )

      # catch NAs in vlss
      # allows variables to be passed that belong to different lists, e.g. state and env when run_leaf is called by the canopy   
      if(any(is.na(vlss))) {
        vlmoss <- vlmoss[-which(is.na(vlss))]
        vlss   <- vlss[-which(is.na(vlss))]
      }

      # print configure setup if requested
      if(.$cpars$cverbose&o) {
        print('', quote=F )
        print('template configure:', quote=F )
        print(df, quote=F )
        print(listnames, quote=F )
        print(moss, quote=F )
        print(slmoss, quote=F )
        print(vlmoss, quote=F )
        print(vlss, quote=F )
        print(which(is.na(vlss)), quote=F )
        print(.[[vlist]], quote=F )
      }

      # assign UQ variables
      if(length(slss)>0)   vapply( slmoss, .$configure_sublist, numeric(1), vlist=vlist, df=df ) 
      if(length(vlmoss)>0) .[[vlist]][vlss] <- df[vlmoss]
    }
 
    
    # configure a list variable 
    configure_sublist <- function(., ss, vlist, df ) {
      lnames <- strsplit(names(df)[ss], '.', fixed=T )
      ss1    <- which(names(.[[vlist]])==lnames[[1]][2])
      ss2    <- which(names(.[[vlist]][[ss1]])==lnames[[1]][3])
      .[[vlist]][[ss1]][ss2] <- df[ss] 
      return(1) 
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

    .test <- function(., verbose=F ) {
      
      if(verbose) str(.)
      
      .$cpars$verbose       <- verbose
      .$cpars$output        <-'full'
      
      .$run()
    }

    .test_change_func <- function(., verbose=F,
                                  template.text='f_text_combine',
                                  template.calcval='f_calcval_product',
                                  template.print='f_print_textonly' ) {
      
      if(verbose) str(.)

      .$cpars$verbose       <- verbose
      .$cpars$output        <-'full'
      
      .$fnames$text         <- template.text
      .$fnames$calcval      <- template.calcval
      .$fnames$print        <- template.print
      
      .$run()
    }

    .test_change_pars <- function(., verbose=T,
                                  template.text1='hello',
                                  template.text2='world' ) {
      
      if(verbose) str(.)

      .$cpars$verbose       <- verbose
      .$cpars$output        <-'full'
      
      .$pars$text1          <- template.text1
      .$pars$text2          <- template.text2
      
      .$run()
    }
      
      
#######################################################################        
# end template object      
})



### END ###
