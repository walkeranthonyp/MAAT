################################
#
# Wrapper object for model objects 
# 
# AWalker December 2015
#
################################

library(proto)
library(parallel)
library(plyr)

source('general_functions.R')


### HIGH LEVEL WRAPPER FUNCTION
#####################################
wrapper_object <-
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'maat'
    
    # expected child objects
    # a proto object that is the model structure in which to run the SA/UQ 
    model <- NULL    
    
    # build function
    .build <- function(.,model) {
      maat       <- as.proto(.$as.list(),parent=.GlobalEnv)
      maat$model <- as.proto(model$as.list(),parent=maat)
      maat
    }
    
    
    
    ###########################################################################
    # Main run function

    run   <- function(.,verbose=T) {
      
      # Initialisation checks
      # need to add a check for equal par vector legnths if this is a UQ run, or add functioanlity that allows a function to be specified that will generate pars
      
      # need to add a checking function that allows the names of fnames, pars, and env to be checked and if they don't match up to print an error message 
      
      # configure initialisation lists - not called if unit testing (need to develop unit testing specifically for the model object init functions)
      if(!.$wpars$unit_testing) .$model$init()
      
      # initialise model with static variables
      .$model$configure(func='write_fnames',df=data.frame(.$static$fnames,stringsAsFactors=F) )
      .$model$configure(func='write_pars',  df=data.frame(.$static$pars,stringsAsFactors=F)   )
      .$model$configure(func='write_env',   df=data.frame(.$static$env,stringsAsFactors=F)    )      
      
      # create UQ dataframe  
      # expand the fnames and driving variables  
      .$dataf$fnames <- if(!is.na(.$vars$fnames[1])) expand.grid(.$vars$fnames,stringsAsFactors=F) else NULL
      .$dataf$env    <- if(!is.na(.$vars$env[1]))    expand.grid(.$vars$env,stringsAsFactors=F)    else NULL
      # bind the parameter vectors OR expand pre-determined parameter vectors 
      # - this could be modified to allow a distribution function to be specifed
      if(.$wpars$UQ) { 
        # yet to solve is how to fix parameters to a relevant function (not sure that would be necessary)
        .$dataf$pars <- if(!is.na(.$vars$pars[1]))   as.data.frame(do.call(cbind,.$vars$pars))     else NULL
      } else {
        .$dataf$pars <- if(!is.na(.$vars$pars[1]))   expand.grid(.$vars$pars,stringsAsFactors=F)   else NULL
      } 
      
      # add an extra column to the dataframes if they have only one column
      # - this prevents them being coerced to a vector in further functions
      if(!is.null(.$dataf$fnames)) if(dim(.$dataf$fnames)[2]==1) .$dataf$fnames <- data.frame(.$dataf$fnames,NA) 
      if(!is.null(.$dataf$pars))   if(dim(.$dataf$pars)[2]==1)   .$dataf$pars   <- data.frame(.$dataf$pars,NA) 
      if(!is.null(.$dataf$env))    if(dim(.$dataf$env)[2]==1)    .$dataf$env    <- data.frame(.$dataf$env,NA) 
      if(!is.null(.$dataf$met))    if(dim(.$dataf$met)[2]==1)    .$dataf$met    <- data.frame(.$dataf$met,NA)       
      
      # calculate dataframe lengths, if no dataframe return 1
      # - used to set the number of iterations in the run functions  
      .$dataf$lf <- ( if(is.null(.$dataf$fnames)|(sum(dim(.$dataf$fnames)==0)==2)) 1 else length(.$dataf$fnames[,1]) )
      .$dataf$lp <- ( if(is.null(.$dataf$pars)|(sum(dim(.$dataf$pars)==0)==2)    ) 1 else length(.$dataf$pars[,1])   )
      .$dataf$le <- ( if(is.null(.$dataf$env)|(sum(dim(.$dataf$env)==0)==2)      ) 1 else length(.$dataf$env[,1])    )
      .$dataf$lm <- ( if(is.null(.$dataf$met)|(sum(dim(.$dataf$met)==0)==2)      ) 1 else length(.$dataf$met[,1])    )
      
      # output summary of maat setup
      .$print_data()
      
      # run model
      .$print_run()
      .$dataf$out <- if(.$wpars$multic) as.data.frame(do.call(rbind,mclapply(1:.$dataf$lf,.$runf,mc.cores=.$wpars$procs)))
                     else               as.data.frame(do.call(rbind,  lapply(1:.$dataf$lf,.$runf)))
      print("output ::",quote=F)
      print(paste('length :',length(.$dataf$out[,1])))
      print(head(.$dataf$out),quote=F)
      print('',quote=F)      
      
    }
    
    
    
    ###########################################################################
    # nested run functions

    printc <- function(.,r1,r2) {
      print(r1,quote=F)
      print(r2,quote=F)
    }
    
    runf <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnames[i,]),F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[i,])
      
      # call next run function
      data.frame(do.call(rbind,lapply(1:.$dataf$lp,.$runp)))      
    }
    
    runp <- function(.,j) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=data.frame(.$dataf$pars[j,]),F)
      if(.$wpars$cverbose) .$printc('pars',.$dataf$pars[j,])
      
      # call next run function
      data.frame(do.call(rbind,lapply(1:.$dataf$le,.$rune)))      
    }
    
    rune <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure environment in the model
      if(!is.null(.$dataf$env)) .$model$configure(func='write_env',df=.$dataf$env[k,],F)
      if(.$wpars$cverbose) .$printc('env',.$dataf$env[k,])
      
      # call next run function
      if(is.null(.$dataf$met)){
        .$model$run()        
      } else {
        data.frame(do.call(rbind,lapply(1:.$dataf$lm,.$model$run_met)))
      }  
    }
    
    
    
    ###########################################################################
    # variables that are to be varied in the SA/UQ
    
    # initialisation lists
#     init_ls <- NA
    init_static  <- NA
    init_dynamic <- NA
    
    # static variables
    # all expected to be of class 'list'
    # each list in the below list comprise a single string or numeric value for each variable, labelled by the variable name
    static <- list( 
      fnames = NA,
      pars   = NA,
      env    = NA
    )
    
    # dynamic variables
    # all expected to be of class 'list'
    # each list in the below list comprise vectors of the values for each variable, labelled by the variable name
    # each of these list is expanded factorially by expand.grid and placed into the below dataframes
    vars <- list( 
      fnames = NA,
      pars   = NA,
      env    = NA
    )
    
    # input/output data frames
    # all expected to be of class 'dataframe'
    dataf  <- list( # currently dataframes but could be matrices
      fnames = NULL,
      pars   = NULL,
      env    = NULL,
      met    = NULL, # a dataframe of sequential meteorological driving data, for running the analysis at a particular site for example 
      obs    = NULL, # a dataframe of observations against which to valiadate/ calculate likelihood of model
      obsse  = NULL, # a dataframe of observation errors for the obs data, must exactly match the above dataframe
      out    = NULL  # output dataframe
    )
    
    # parameters specific to the wrapper object
    wpars <- list(
      multic   = F,  # multicore the simulation
      procs    = 6,  # number of processors to use if multic = T
      cverbose = F,  # write configuration output during runtime 
      UQ       = F,  # run a UQ analysis
      unit_testing = F
    )
    
    
    
    ###########################################################################
    # Output processing functions
    
    combine <- function(.,i,df) suppressWarnings(data.frame(.$dataf$met,df[i,]))
        
    output  <- function(.){
      # function that combines the "vars", "met", and "out" dataframes correctly for output 
      return(
        # if vars
        if(is.null(.$dataf$env)+is.null(.$dataf$pars)+is.null(.$dataf$fnames) < 3) {
          vpars <- if(is.null(.$dataf$pars))  NULL  else if(.$wpars$UQ) .$vars$pars[1] else .$vars$pars
          if(is.null(.$dataf$env))    .$vars$env    <- NULL
          if(is.null(.$dataf$fnames)) .$vars$fnames <- NULL
             
          vardf <- expand.grid(c(.$vars$env,vpars,.$vars$fnames),stringsAsFactors=F)        
          
          # if no met data
          if(is.null(.$dataf$met)) {
            cbind(vardf,.$dataf$out)

          # if met data
          } else  {
            odf <- cbind(do.call(rbind , lapply(1:length(vardf[,1]) , .$combine ,df=vardf )), .$dataf$out )            
            print(head(odf))
            if(dim(vardf)[2]==1) names(odf)[which(names(odf)=='df.i...')] <- names(vardf)
            rm(vardf)
            odf
          } 

        # if no vars  
        } else {
          # if met data
          if(!is.null(.$dataf$met)) cbind( .$dataf$met , .$dataf$out ) else .$dataf$out  
        }
      )
    }
    
    print_data <- function(.){
      print('',quote=F)
      print("MAAT :: summary of data",quote=F)
      print('',quote=F)
      if(is.null(.$dataf$met)) {
        print(paste('ensemble number:',.$dataf$lf*.$dataf$lp*.$dataf$le*.$dataf$lm),quote=F)        
      } else {
        print(paste('ensemble number:',.$dataf$lf*.$dataf$lp*.$dataf$le),quote=F)        
        print(paste('timesteps in met data:',.$dataf$lm),quote=F)                
        print(paste('total number of timesteps:',.$dataf$lf*.$dataf$lp*.$dataf$le*.$dataf$lm),quote=F)                
      }
      print('',quote=F)
      print('',quote=F)
      print("fnames ::",quote=F)
      print(summary(.$dataf$fnames),quote=F)
      print('',quote=F)
      print("pars ::",quote=F)
      print(summary(.$dataf$pars),quote=F)
      print('',quote=F)
      print("env ::",quote=F)
      print(summary(.$dataf$env),quote=F)
      print('',quote=F)
      print("met data ::",quote=F)
      print(summary(.$dataf$met),quote=F)
      print('',quote=F)
    }
    
    print_run <- function(.) {
      print('',quote=F)
      print("MAAT :: run model",quote=F)
      print('',quote=F)
      if(.$wpars$multic) {
        print(paste('parallel processing over ::',.$wpars$procs,'cores'),quote=F)        
        print('',quote=F)
      }
    }
    
    
    
    ######################################################################################
    # Unit testing functions
    
    .test <- function(.,metd=T,mc=T,pr=4,oconf=F) {
      
      # source directory
      source('leaf_object.R')
      library(lattice)
      
      # clone the model object
      .$model <- as.proto(leaf_object$as.list(),parent=.)
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- oconf      
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc  # multicore the ensemble
      .$wpars$procs        <- pr  # number of cores to use if above is true
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if faslse a fully factorial ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening 
      
      ### Define meteorological and environment dataset
      ###############################
      # can load a met dataset here
      # below a trivial met dataset is created to be used as an example
#       metdata <- expand.grid(list(leaf.par = 500,leaf.ca_conc = seq(10,1200,10)))      
      metdata <- expand.grid(list(leaf.par = seq(0,1000,100),leaf.ca_conc = 400))      
      if(metd) .$dataf$met <- metdata
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length,
      # if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
      
      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$vars$fnames <- list(
        leaf.etrans = c('f_j_farquhar1980','f_j_collatz1991'),
        leaf.rs     = c('f_r_zero','f_rs_medlyn2011')
      )
        
      .$vars$env <- list(
        leaf.vpd  = c(1,2),
        leaf.temp = c(5,20)
      )
        
      .$model$fnames$vcmax <- 'f_vcmax_lin'
      .$vars$pars <- list(
        leaf.avn_25 = 9:11,
        leaf.bvn_25 = 4:6
      )
      
      # Run model
      st <- system.time(
        .$run()
      )
      print('',quote=F)
      print('Run time:',quote=F)
      print(st)
      print('',quote=F)
      
      # process & record output
      df <- .$output()
#       p1 <- xyplot(A~ci|leaf.etrans*leaf.rs,df,type='l',auto.key=T)
      p1 <- xyplot(A~leaf.ca_conc|leaf.etrans*leaf.rs,df,groups=leaf.temp,type='l',auto.key=T,
                   panel=function(...) { panel.abline(h=seq(0,20,2.5)) ; panel.xyplot(...) })
      list(df,p1)
    }    
    
    ###########################################################################
    # end maat wrapper 
}) 











