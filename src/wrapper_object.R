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
    
    # clean function to reset object
    clean <- function(.) {
      for(d in 1:length(.$dataf)) if(!is.null(.$dataf[[d]])) .$dataf[[d]] <- NA
      gc()
    }
    
    ###########################################################################
    # Main run function

    run   <- function(.,verbose=T) {
      
      # Initialisation checks
      # need to add a check for equal par vector legnths if this is a UQ run, or add functionality that allows a function to be specified that will generate pars
      
      # need to add a checking function that allows the names of fnames, pars, and env to be checked and if they don't match up to print an error message 
      
      # configure initialisation lists - not called if unit testing (need to develop unit testing specifically for the model object init functions)
      if(!.$wpars$unit_testing) .$model$init()
      
      # initialise model with static variables
      # .$model$configure(func='write_fnames',df=data.frame(.$static$fnames,stringsAsFactors=F) )
      # .$model$configure(func='write_pars',  df=data.frame(.$static$pars,stringsAsFactors=F)   )
      # .$model$configure(func='write_env',   df=data.frame(.$static$env,stringsAsFactors=F)    )      
      .$model$configure(func='write_fnames',df=t(as.matrix(.$static$fnames,stringsAsFactors=F))[1,] )
      .$model$configure(func='write_pars',  df=t(as.matrix(.$static$pars,stringsAsFactors=F))[1,]   )
      .$model$configure(func='write_env',   df=t(as.matrix(.$static$env,stringsAsFactors=F))[1,]    )      
      
      # create dataframes of runtime variables  
      # expand the fnames and driving variables  
      # .$dataf$fnames  <- if(!is.na(.$vars$fnames[1])) expand.grid(.$vars$fnames,stringsAsFactors=F) else NULL
      # .$dataf$env     <- if(!is.na(.$vars$env[1]))    expand.grid(.$vars$env,stringsAsFactors=F   ) else NULL
      .$dataf$fnames  <- if(!is.na(.$vars$fnames[1])) as.matrix(expand.grid(.$vars$fnames,stringsAsFactors=F)) else NULL
      .$dataf$env     <- if(!is.na(.$vars$env[1]))    as.matrix(expand.grid(.$vars$env,stringsAsFactors=F   )) else NULL
      if(.$wpars$UQ) {
        # if an SA/UQ run
        if(.$wpars$UQtype=='saltelli') {
          # if Saltelli style Sobol
          if(.$wpars$eval_strings) {
            # sample parameters from character string code snippets
            n <- 2 * .$wpars$n * .$wpars$nmult
            .$vars$pars <- lapply(.$vars$pars_eval,function(cs) eval(parse(text=cs)))
          }
          if(is.na(.$vars$pars[1] )) stop('wrapper: pars list in vars list is empty')
          .$dataf$pars <- do.call(cbind,.$vars$pars )

        } else if(.$wpars$UQtype=='ye') {
          # ye method variable matrices are created in the first loop of the alogorithm
          
          # due to different parameter sample numbers in process A and B loops,
          # parameters samples must be generated from code snippets as strings
          .$wpars$eval_strings <- T
          if(is.na(.$vars$pars_eval)) {
            stop('wrapper: Ye method SA must draw parameter samples during runtime \n
                  from code snippets expressed as strings in vars$pars_eval') 
          }
          # check input vars$pars* are same length
          test_in <- length(.$vars$pars) - length(.$vars$pars_proc)
          if(test_in!=0) stop('wrapper: Parameter input vectors - pars & pars_proc - are not the same length')
          test_in <- length(.$vars$pars_eval) - length(.$vars$pars_proc)
          if(test_in!=0) stop('wrapper: Parameter input vectors - pars_eval & pars_proc - are not the same length')
          # check input vars$pars* elements have same names
          # - to be done
          
        } else {
          # alternative UQ methods
          stop(paste('wrapper: no method for SA/UQ type',.$wpars$UQtype))
        }
        
      } else {
        # not a formal SA/UQ run - just a factorial combination of variables specified in the vars lists  
        # .$dataf$pars  <- if(!is.na(.$vars$pars[1]))   expand.grid(.$vars$pars,stringsAsFactors=F)   else NULL
        .$dataf$pars <- if(!is.na(.$vars$pars[1])) as.matrix(expand.grid(.$vars$pars,stringsAsFactors=F))   else NULL
      } 
      
      # add an extra column to the dataframes if they have only one column
      # - this prevents them being coerced to a vector in further functions
      # if(!is.null(.$dataf$fnames)) if(dim(.$dataf$fnames)[2]==1) .$dataf$fnames <- data.frame(.$dataf$fnames,NA) 
      # if(!is.null(.$dataf$pars))   if(dim(.$dataf$pars)[2]==1)   .$dataf$pars   <- data.frame(.$dataf$pars,NA) 
      # if(!is.null(.$dataf$env))    if(dim(.$dataf$env)[2]==1)    .$dataf$env    <- data.frame(.$dataf$env,NA) 
      # if(!is.null(.$dataf$met))    if(dim(.$dataf$met)[2]==1)    .$dataf$met    <- data.frame(.$dataf$met,NA)       
      
      # calculate dataframe lengths, if no dataframe return 1
      # - used to set the number of iterations in the run functions  
      # - parameter lengths 
      if(.$wpars$UQ&.$wpars$UQtype=='ye') {
        # determine number of processes to be analaysed
        .$dataf$lf <- length(.$vars$fnames)
      } else {
        # any type of run other than Ye process sensitivity analysis 
        .$dataf$lf <- if(is.null(.$dataf$fnames)|(sum(dim(.$dataf$fnames)==0)==2)) 1 else length(.$dataf$fnames[,1]) 
        .$dataf$lp <- if(is.null(.$dataf$pars)|(sum(dim(.$dataf$pars)==0)==2) )    1 else length(.$dataf$pars[,1])
      }
      .$dataf$le <- if(is.null(.$dataf$env)|(sum(dim(.$dataf$env)==0)==2)      ) 1 else length(.$dataf$env[,1])    
      .$dataf$lm <- if(is.null(.$dataf$met)|(sum(dim(.$dataf$met)==0)==2)      ) 1 else length(.$dataf$met[,1])    
      
      # output summary of maat setup
      .$print_data()
      .$print_run()


      
      # Call initial run function in the hierarchy of nested run functions
      
      # process UQ run, or either general variable run or matrix A and B of Saltelli method
      if(.$wpars$UQ&.$wpars$UQtype=='ye') {
        # Ye process SA is not multicored at this stage as multicoring here messes with the processes A and B in the data structure  
        lapply(1:.$dataf$lf,.$run_general_process_SA)
      } else {
        # Saltelli SA matrix AB run, or factorial run
        .$dataf$out <-
          as.data.frame(
            do.call(
              rbind, {
                # multicore or not
                if(.$wpars$multic) mclapply(1:.$dataf$lf,.$runf,mc.cores=.$wpars$procs)
                else                 lapply(1:.$dataf$lf,.$runf)
              }
            ))
      }
        
      # if Saltellis SA generate output from parameter matrices ABi
      if(.$wpars$UQtype=='saltelli') {
        # prob should run output here and clear out df 
        print('Saltelli matrix AB completed', quote=F)
        print('run Saltelli array ABi',quote=F)
        print('',quote=F)
        
        .$dataf$out_saltelli <-
          if(.$wpars$multic) mclapply(1:.$dataf$lf,.$runf_saltelli,mc.cores=.$wpars$procs)
          else                 lapply(1:.$dataf$lf,.$runf_saltelli)          
        
        print('Saltelli array ABi completed',quote=F)
        print('',quote=F)
      }             
      
      # print summary of results
      # - output function needs to be called from run script
      print("output ::",quote=F)
      print(paste('length :',length(.$dataf$out[,1])))
      print(head(.$dataf$out),quote=F)
      print('',quote=F)      
      
    }
    
    
    
    ###########################################################################
    # nested run functions - 
    ###########################################################################

    printc <- function(.,r1,r2) {
      print(r1,quote=F)
      print(r2,quote=F)
    }

    # for factorial runs
    ###########################################################################
    
    runf <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call runp
      
      # configure function names in the model
      # if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnames[i,]),F)
      if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=.$dataf$fnames[i,],F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[i,])
      
      # call next run function
#       data.frame(do.call(rbind,lapply(1:.$dataf$lp,.$runp)))
      data.frame(
        do.call(
          rbind, {
            if(.$wpars$multic) mclapply(1:.$dataf$lp,.$runp,mc.cores=max(1,floor(.$wpars$procs/.$dataf$lf)) )
            else                 lapply(1:.$dataf$lp,.$runp)
          }
      ))
    }
    
    runp <- function(.,j) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
      # assumes that each row of the pars matrix are independent and non-sequential
      # call rune
      
      # configure parameters in the model
      # if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=data.frame(.$dataf$pars[j,]),F)
      if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=.$dataf$pars[j,],F)
      if(.$wpars$cverbose) .$printc('pars',.$dataf$pars[j,])
      
      # call next run function
      data.frame(do.call(rbind,lapply(1:.$dataf$le,.$rune)))      
    }
    
    rune <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$env matrix to the model
      # assumes that each row of the env matrix are independent and non-sequential
      # call .$model$run or .$model$run_met if met data are provided
      
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

    
    # for ABi array for Sobol SA using Saltelli method
    ###########################################################################
    
    runf_saltelli <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call rune_saltelli
      
      # configure function names in the model
      # if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnames[i,]),F)
      if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=.$dataf$fnames[i,],F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[i,])
      
      # call next run function
      if(.$wpars$multic) mclapply(1:.$dataf$le,.$rune_saltelli,mc.cores=max(1,floor(.$wpars$procs/.$dataf$lf)) )
      else                 lapply(1:.$dataf$le,.$rune_saltelli)
    }

    rune_saltelli <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$env matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call runpmat_saltelli
      
      # configure environment in the model
      if(!is.null(.$dataf$env)) .$model$configure(func='write_env',df=.$dataf$env[k,],F)
      if(.$wpars$cverbose) .$printc('env',.$dataf$env[k,])
      
      # call parameter matrix run function
      if(is.null(.$dataf$met)){
        omatl <- lapply(1:dim(.$dataf$pars)[2],.$runpmat_saltelli)
        # convert to 3 dim array
        array( unlist(omatl) , c(dim(omatl[[1]]),length(omatl)) )
      } else {
        print('Met data run not yet supported with Sobol',quote=F)
        stop
      }  
    }
    
    runpmat_saltelli <- function(.,p) {
      # This wrapper function is called from an lapply or mclappy function to run over each column of the dataf$pars matrix
      # call runp_saltelli

      # returns a numeric matrix
      do.call(rbind,lapply((.$wpars$n+1):.$dataf$lp,.$runp_saltelli,pk=p))      
    }
    
    runp_saltelli <- function(.,j,pk) {
      # This wrapper function is called from an lapply or mclappy function 
      # wrapper creates a parameter matrix ABi and passes every row of the matrix to the model
      # assumes that each row of the matrix are independent and non-sequential
      # call .$model$run
      
      # create index matrix to create row on matrix ABi for the .$dataf$par matrix (which is matrix A stacked on top of matrix B)
      sub     <- rep(j-.$wpars$n,dim(.$dataf$pars)[2])
      sub[pk] <- j
      smat    <- cbind(sub,1:dim(.$dataf$pars)[2])      
      
      # create a dataframe from vector and add names
      # psdf        <- data.frame(t(.$dataf$pars[smat]))
      psdf        <- t(.$dataf$pars[smat])
      names(psdf) <- names(.$dataf$pars)
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=psdf,F)
      if(.$wpars$cverbose) .$printc('pars',psdf)
      
      # run model
      # - rune_saltelli requires a numeric vector for the output from runpmat_saltelli
      # - warnings are suppressed because character strings in model output are converted to NAs which comes with a warning   
      suppressWarnings(as.numeric( .$model$run() ))        
    }
    
    
    
    # nested run functions for Process Sensitivity Analysis after Ye etal 201x 
    ###########################################################################

        # The Ye method is a nested system of loops implemented by apply type functions
    # Loop 1: switch process A process B loop         
    # Loop 2: process loop for process A               - use standard location for variable function values
    # Loop 3: parameter loop for process A             - use standard location for variable parameter values
    # Loop 4: process representation loop for proces B - use an additional non-standard location for variable function values
    # Loop 5: parameter loop for proces B              - use an additional non-standard location for variable parameter values
    # Loop 6: environment loop                         - use standard location for variable environment values
    
    run_general_process_SA <- function(.,f) {
      # this function is the overall wrapper function to run a generic process sensitivity analysis
      # The principle is to create a loop that runs the process SA nested loops once for each process to be analysed
      # This separates out the process in question - process A - from the other process(es) - process B (can be more than one process).
      # This function partitions the parameters to process A and and process B, 
      # then creates the fnames and pars dataframes for process A and process B,
      # calls the below run function,
      # outputs an .RDS for each process segregation

      # create the fnames dataframes for process A and process B
      .$dataf$fnames  <- if(!is.na(.$vars$fnames[f]))  as.matrix(expand.grid(.$vars$fnames[f] ,stringsAsFactors=F)) else stop()
      .$dataf$fnamesB <- if(!is.na(.$vars$fnames[-f])) as.matrix(expand.grid(.$vars$fnames[-f],stringsAsFactors=F)) else stop()
      
      # determine the number of the rows in process dataframes
      .$dataf$lfA <- if(is.null(.$dataf$fnames )|(sum(dim(.$dataf$fnames )==0)==2)) 1 else length(.$dataf$fnames[,1]) 
      .$dataf$lfB <- if(is.null(.$dataf$fnamesB)|(sum(dim(.$dataf$fnamesB)==0)==2)) 1 else length(.$dataf$fnamesB[,1]) 
      
      # partition the parameters to process A and and process B
      .$procA_name <- names(.$vars$fnames)[f]
      .$procA_subs <- which(unlist(.$vars$pars_proc)==.$procA_name)

      # evaluate parameter strings to sample vectors
      # - this allows a distribution function to be specifed
      # - also allows the dynamic calcuation of n for process A and B parameter samples
      if(.$wpars$eval_string) {
        n <- .$wpars$n
        .$vars$pars[.$procA_subs ] <- lapply(.$vars$pars_eval[.$procA_subs ],function(cs) eval(parse(text=cs)))
        n <- .$dataf$lfA * .$dataf$lfB * .$wpars$n^2
        .$vars$pars[-.$procA_subs] <- lapply(.$vars$pars_eval[-.$procA_subs],function(cs) eval(parse(text=cs)))
      }
      
      # bind the parameter vectors 
      # .$dataf$pars  <- if(!is.na(.$vars$pars[1])) as.data.frame(do.call(cbind,.$vars$pars[.$procA_subs] )) else stop()
      # .$dataf$parsB <- if(!is.na(.$vars$pars[2])) as.data.frame(do.call(cbind,.$vars$pars[-.$procA_subs])) else stop()
      .$dataf$pars  <- if(!is.na(.$vars$pars[1])) do.call(cbind,.$vars$pars[.$procA_subs] ) else stop()
      .$dataf$parsB <- if(!is.na(.$vars$pars[2])) do.call(cbind,.$vars$pars[-.$procA_subs]) else stop()

      # determine the number of the rows in parameter dataframes
      .$dataf$lp  <- .$wpars$n # convert these to be the row number of the actual matrices
      .$dataf$lpB <- .$dataf$lfA * .$dataf$lfB * .$wpars$n^2
      
      # # add an extra column to the dataframes if they have only one column
      # # - this prevents them being coerced to a vector in further functions
      # if(!is.null(.$dataf$fnames))  if(dim(.$dataf$fnames)[2]==1)  .$dataf$fnames  <- data.frame(.$dataf$fnames,NA)
      # if(!is.null(.$dataf$fnamesB)) if(dim(.$dataf$fnamesB)[2]==1) .$dataf$fnamesB <- data.frame(.$dataf$fnamesB,NA)
      # if(!is.null(.$dataf$pars))    if(dim(.$dataf$pars)[2]==1)    .$dataf$pars    <- data.frame(.$dataf$pars,NA)
      # if(!is.null(.$dataf$parsB))   if(dim(.$dataf$parsB)[2]==1)   .$dataf$parsB   <- data.frame(.$dataf$parsB,NA)
      
      # call the below run function
      .$dataf$out <-
        data.frame(
          do.call(rbind,
                  if(.$wpars$multic) mclapply(1:.$dataf$lfA, .$run_repA, mc.cores=.$wpars$procs)
                  else                 lapply(1:.$dataf$lfA, .$run_repA)
          ))
      
      # output an .RDS for each process AB combination
      # print(.$output())
      
      # process & record output
      if(.$wpars$unit_testing) { 
        hd     <- getwd()
        setwd('~/tmp')
        ofname <- 'Ye_test' 
      } else setwd(odir)

      write_to_file(.$output(),paste(ofname,'proc',f,sep='_'),type='rds')  
      setwd(hd)
    }
    
    run_repA <- function(.,g) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call run_parA
      
      # configure function names in the model
      # print(.$dataf$fnames)
      # print(.$dataf$fnames[g,])
      # if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnames[g,]),F)
      if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=.$dataf$fnames[g,],F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[g,])
      
      # call process A parameter run function
      # data.frame(do.call(rbind,lapply(1:.$dataf$lp,.$run_parA,offset=g)))      
      data.frame(
        do.call(rbind,
                if(.$wpars$multic) mclapply(1:.$dataf$lp,.$run_parA,offset=g,mc.cores=max(floor(.$wpars$procs/.$dataf$lfA),1) )
                else                 lapply(1:.$dataf$lp,.$run_parA,offset=g)
      ))
    }
        
    run_parA <- function(.,h,offset) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
      # assumes that each row of the pars matrix are independent and non-sequential
      # call run_repB
      
      # configure parameters in the model
      # if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=data.frame(.$dataf$pars[h,]),F)
      if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=.$dataf$pars[h,],F)
      if(.$wpars$cverbose) .$printc('pars',.$dataf$pars[h,])
      
      # call process B process representation run function
      # os is the row in the factorial ftypes and pars matrix
      os <- .$wpars$n*(offset - 1) + h     
      data.frame(do.call(rbind,lapply(1:.$dataf$lfB,.$run_repB,offset=os)))      
    }

    run_repB <- function(.,i,offset) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnamesB matrix to the model
      # assumes that each row of the fnamesB matrix are independent and non-sequential
      # call run_parB
      
      # configure function names in the model
      # if(!is.null(.$dataf$fnamesB)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnamesB[i,]),F)
      if(!is.null(.$dataf$fnamesB)) .$model$configure(func='write_fnames',df=.$dataf$fnamesB[i,],F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnamesB[i,])
      
      # call process B parameter run function
      # os is the row in the parsB matrix
      os <- .$dataf$lfB*(offset - 1) + i    
      data.frame(do.call(rbind,lapply(os:(os+.$wpars$n-1),.$run_parB)))      
    }
    
    run_parB <- function(.,j) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$parsB matrix to the model
      # assumes that each row of the parsB matrix are independent and non-sequential
      # call rune
      
      # configure parameters in the model
      # if(!is.null(.$dataf$parsB)) .$model$configure(func='write_pars',df=data.frame(.$dataf$parsB[j,]),F)
      if(!is.null(.$dataf$parsB)) .$model$configure(func='write_pars',df=.$dataf$parsB[j,],F)
      if(.$wpars$cverbose) .$printc('pars',.$dataf$parsB[j,])
      
      # call the environment run function to loop over CO2 values
      data.frame(do.call(rbind,lapply(1:.$dataf$le,.$rune)))
    }
        
    
    
    ###########################################################################
    # variables that are to be varied in the SA/UQ
    
    # initialisation list
    init_ls <- NA
    
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
    # each list in the 'vars' list comprise vectors of the values for each variable, labelled by the variable name
    # each of these lists is expanded factorially by expand.grid and placed into the below list of dataframes
    vars <- list( 
      fnames    = NA,
      fnamesB   = NA,
      pars      = NA,
      # list with same elements and names as pars but giving the fnames list name i.e. the process name to which each parameter belongs
      pars_proc = NA,
      # list with same elements and names as pars but each element is a code snippet as a string that once evaluated gives a vector or parameter values
      # allows different types of distributions to be specified for each parameter
      # this must be used for the Ye SA method
      pars_eval = NA,
      parsB     = NA,
      env       = NA
    )
    
    # input/output matrices and dataframes
    # with an associated length for input matrices
    dataf  <- list( 
      # variables matrices - created during runtime 
      fnames  = NULL,
      fnamesB = NULL,
      pars    = NULL,
      parsB   = NULL,
      env     = NULL,
      met     = NULL,          # a dataframe of sequential meteorological driving data, for running the analysis at a particular site for example 
      # row length of matrices
      lf      = NULL,
      lfA     = NULL,
      lfB     = NULL,
      lp      = NULL,
      lpB     = NULL,
      le      = NULL,
      lm      = NULL,
      # output dataframes
      # - maintained as dataframes as model output can be character strings
      obs     = NULL,          # a dataframe of observations against which to valiadate/ calculate likelihood of model
      obsse   = NULL,          # a dataframe of observation errors for the obs data, must exactly match the above dataframe
      out     = NULL,          # output dataframe
      out_saltelli = NULL     # saltelli output list 
    )
    
    # parameters specific to the wrapper object
    wpars <- list(
      multic   = F,           # multicore the simulation
      procs    = 6,           # number of processors to use if multic = T
      cverbose = F,           # write configuration output during runtime 
      UQ       = F,           # run a UQ analysis
      UQtype   = 'none',      # SA/UQ type - 'saltelli' and 'ye' available so far
      n        = numeric(1),  # parameter sample number
      nmult    = 1,           # parameter sample number multiplier for saltelli method
      eval_strings = F,       # switch tellin wrapper that vars$pars are to be evaluated from code string snippets in vars$pars_eval
      unit_testing = F
    )
    
    
    
    # Output processing functions
    ###########################################################################
    
    combine <- function(.,i,df) suppressWarnings(data.frame(.$dataf$met,df[i,]))
    
    output  <- function(.){
      # function that combines the "vars", "met", and "out" dataframes correctly for output 
      return(
        # if at least one of fnames, pars, and env are varied
        if(is.null(.$dataf$env)+is.null(.$dataf$pars)+is.null(.$dataf$fnames) < 3) {
#           vpars <- if(is.null(.$dataf$pars))  NULL  else if(.$wpars$UQ) .$vars$pars[1] else .$vars$pars
          # vpars    <- if(is.null(.$dataf$pars))    NULL else .$vars$pars[[1]]
          # vpars    <- if(is.null(.$dataf$pars))    NULL else .$vars$pars
          # venv     <- if(is.null(.$dataf$env))     NULL else .$vars$env
          # vfnames  <- if(is.null(.$dataf$fnames))  NULL else .$vars$fnames
          # vfnamesB <- if(is.null(.$dataf$fnamesB)) NULL else .$vars$fnamesB
          # vparsB  <- if(.$wpars$UQ) .$vars$parsB[[1]]  else NULL
          
          # run types if
          if(.$wpars$UQ) {
            if(.$wpars$UQtype=='ye') {
              
              # if Ye (i.e. process) UQ - need to write a function to do this
              # the number of rows in the resultant dataframe should be lf*lfB*n^2*le
              
              # convert dfs to matrices - this should probably be done at the beginning due to efficiency of matrices vs dataframes
              # vfnames  <- as.matrix(.$dataf$fnames)
              # vfnamesB <- as.matrix(.$dataf$fnamesB)
              # vpars    <- as.matrix(.$dataf$pars) 
              # vparsB   <- as.matrix(.$dataf$parsB)
              # venv     <- as.matrix(.$dataf$env)
              # vfnames  <- .$dataf$fnames
              # vfnamesB <- .$dataf$fnamesB
              # vpars    <- .$dataf$pars 
              # vparsB   <- .$dataf$parsB
              # venv     <- .$dataf$env
              
              # combine input into a single dataframe in order of output dataframe,
              #  - i.e. repeats lines in input dataframes/matrices to align with output
              vardf    <- cbind(
                # fnames of process A - length lf * lfB * le * n^2
                apply(.$dataf$fnames, 2,function(v) rep(v,each=.$dataf$lfB*.$dataf$le*.$wpars$n^2) ),
                # pars of process A
                apply(.$dataf$pars,   2,function(v) rep(rep(v,each=.$dataf$lfB*.$dataf$le*.$wpars$n),.$dataf$lfA) ),
                # fnames of process B
                apply(.$dataf$fnamesB,2,function(v) rep(rep(v,each=.$dataf$le*.$wpars$n),.$dataf$lfA*.$wpars$n) ),
                # parsB
                apply(.$dataf$parsB,  2,function(v) rep(v,each=.$dataf$le) ),
                # environment
                apply(.$dataf$env,    2,function(v) rep(v,.$dataf$lfA*.$dataf$lfB*.$wpars$n^2) )
              )
              
              # convert to dataframe
              vardf        <- as.data.frame(vardf)
              
            } else if(.$wpars$UQtype=='saltelli') { 
              # if Saltelli SA
              # this is not currently used by the run scripts as the Saltelli method expects a different organisation of output
              # can be called by user to inspect Saltelli run output from matrix AB in a user readable format

              # vpars   <- .$dataf$pars
              # vpars   <- as.data.frame(apply(vpars,2,function(v) rep(v,each=.$dataf$le) ))
              # vardf   <- suppressWarnings(cbind( 
              #   rep(unlist(vfnames),each=.$dataf$le*2*.$wpars$n),
              #   vpars,
              #   unlist(venv) 
              # ))
              
              vardf   <- cbind( 
                # fnames
                apply(.$dataf$fnames, 2,function(v) rep(v,each=.$dataf$le*2*.$wpars$n) ),
                # pars - I think this has the correct length specifications, need to check witha saltelli unit testing function
                apply(.$dataf$pars,   2,function(v) rep(rep(v,each=.$dataf$le),.$dataf$lf) ),
                # environment
                apply(.$dataf$env,    2,function(v) rep(v,.$dataf$lf*2*.$wpars$n) )
              )
              
              vardf        <- as.data.frame(vardf)
              # names(vardf) <- c(names(vfnames),names(.$dataf$pars),names(venv))
            }   
          } else { 
            # if factorial combination run
            vpars    <- if(is.null(.$dataf$pars))    NULL else .$vars$pars
            venv     <- if(is.null(.$dataf$env))     NULL else .$vars$env
            vfnames  <- if(is.null(.$dataf$fnames))  NULL else .$vars$fnames

            vardf <- expand.grid(c(venv,vpars,vfnames),stringsAsFactors=F)        
          }
          
          # if no met data
          if(is.null(.$dataf$met)) {
            # suppressWarnings(cbind(vardf,.$dataf$out))
            cbind(vardf,.$dataf$out)
            
          # if met data
          } else {
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
    
    output_saltelli <- function(.) {
      # creates output for a saltelli Sobol sensitivity analysis 
      # A and B matrices are stacked in a single matrix, which for each model and environment combination are then stored in an array 

      # AB output is an array
      # - dim 1 (rows)   model combination
      # - dim 2 (cols)   environment combination
      # - dim 3 (slices) sample
      # - dim 4          output variable (character variables are coerced to NAs)
      
      # ABi output is a nested list
      # - list 1 model combination
      # - list 2 environment combination
      # - each element of list 2 is a 3 dimesional arrays 
      #    - dim 1 (rows)   sample
      #    - dim 2 (cols)   output variable
      #    - dim 3 (slices) parameter that has used value from matrix B while all other par values are from matrix A
      
      # create AB output matrix array
      # - suppressWarnings on the as.numeric call as there can be character strings in .$dataf$out which cause warnings when converted to NAs
      ab_mat_out  <- array(suppressWarnings(as.numeric(unlist(.$dataf$out))),c(.$dataf$le,2*.$wpars$n,.$dataf$lf,length(.$dataf$out)))
      tab_mat_out <- aperm(ab_mat_out,c(3,1,2,4))
      
      list(AB=tab_mat_out,ABi=.$dataf$out_saltelli)
    }
    
    # Print functions
    ###########################################################################
    
    print_data <- function(.) {
      print('',quote=F)
      print("MAAT :: summary of data",quote=F)
      print('',quote=F)

      ens_n <- 
        if(.$wpars$UQ) {
          if(.$wpars$UQtype=='ye') .$dataf$lf*.$dataf$lfB*.$wpars$n^2*.$dataf$le
          else                     .$dataf$lf*.$wpars$n*(2+dim(.$dataf$pars)[2])*.$dataf$le
        } else .$dataf$lf*.$dataf$lp *.$dataf$le
      
      print(paste('ensemble number:',ens_n),quote=F)
      if(!is.null(.$dataf$met)) {
        print(paste('timesteps in met data:',.$dataf$lm),quote=F)                
        print(paste('total number of timesteps:',ens_n*.$dataf$lm),quote=F)                
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
    
    # general factorial test, with or without metdata
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
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening 
      
      ### Define meteorological and environment dataset
      ###############################
      # can load a met dataset here
      # below a trivial met dataset is created to be used as an example
#       metdata <- expand.grid(list(leaf.par = 500,leaf.ca_conc = seq(10,1200,10)))      
      metdata <- as.matrix(expand.grid(list(leaf.par = seq(0,1000,100),leaf.ca_conc = 400)))      
      if(metd) .$dataf$met <- metdata
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length,
      # if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
      
      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$static$fnames <- list(vcmax='f_vcmax_lin')
      .$vars$fnames <- list(
        leaf.etrans = c('f_j_farquhar1980','f_j_collatz1991'),
        leaf.rs     = c('f_r_zero','f_rs_medlyn2011')
      )
        
      .$vars$env <- list(
        leaf.vpd  = c(1,2),
        leaf.temp = c(5,20)
      )
        
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

    # test function for Ye method Sobol process sensitivity analysis
    .test_ye <- function(.,metd=F,mc=T,pr=4,oconf=F,n=3) {
      
      # source directory
      source('leaf_object.R')
      library(lattice)
      
      # clone the model object
      .$model <- as.proto(leaf_object$as.list(),parent=.)
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- oconf      
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc   # multicore the ensemble
      .$wpars$procs        <- pr   # number of cores to use if above is true
      .$wpars$UQ           <- T    # run a UQ/SA style ensemble 
      .$wpars$UQtype       <- 'ye' # Ye style SA ensemble 
      .$wpars$unit_testing <- T    # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions) 
      .$wpars$n            <- n    # number of parameter samples in each loop 
      .$wpars$eval_strings <- T    # parameters are passed as strings to be evaluated to allow for different sample numbers 
      .$coef_var           <- 0.1
      
      ### Define static variables 
      ###############################
      .$model$env$par <- 1000
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # "pars" lists must contain parameter vectors that are of equal length,

      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$static$fnames <- list(vcmax='f_vcmax_lin')
      .$vars$fnames <- list(
        leaf.Alim   = c('f_lim_farquhar1980','f_lim_collatz1991'),
        leaf.etrans = c('f_j_farquharwong1984','f_j_collatz1991','f_j_harley1992')
      )

      # .$vars$fnamesB <- list(
      #   leaf.etrans = c('f_j_farquhar1980','f_j_collatz1991')
      # )
      
      .$vars$pars <- list(
        leaf.avn_25   = NA,
        leaf.bvn_25   = NA,
        leaf.theta    = NA,
        leaf.e_ajv_25 = NA
      )

      .$vars$pars_eval <- list(
        leaf.avn_25   = ' 10 * rnorm(n,1,.$coef_var)',
        leaf.bvn_25   = '  5 * rnorm(n,1,.$coef_var)',
        leaf.theta    = '0.9 * rnorm(n,1,.$coef_var)',
        leaf.e_ajv_25 = '0.9 * rnorm(n,1,.$coef_var)'
      )
      
      .$vars$pars_proc <- list(
        leaf.avn_25   = 'leaf.Alim',
        leaf.bvn_25   = 'leaf.Alim',
        leaf.theta    = 'leaf.etrans',
        leaf.e_ajv_25 = 'leaf.etrans'
      ) 

      .$vars$env <- list(
        leaf.ca_conc  = c(400,600)
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
      p1 <- bwplot(leaf.ca_conc~A|leaf.Alim*leaf.etrans,df,auto.key=T,abline=0)
      list(df,p1)
    }    
    
    # test function for Saltelli method Sobol parametric sensitivity analysis
    .test_saltelli <- function(.,metd=F,mc=T,pr=4,oconf=F,n=3,eval_strings=T) {
      
      # source directory
      source('leaf_object.R')
      library(lattice)
      
      # clone the model object
      .$model <- as.proto(leaf_object$as.list(),parent=.)
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- oconf      
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc           # multicore the ensemble
      .$wpars$procs        <- pr           # number of cores to use if above is true
      .$wpars$UQ           <- T            # run a UQ/SA style ensemble 
      .$wpars$UQtype       <- 'saltelli'   # Saltelli style SA ensemble 
      .$wpars$unit_testing <- T            # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions) 
      .$wpars$n            <- n            # number of parameter samples in each loop 
      .$wpars$eval_strings <- eval_strings # parameters are passed as strings to be evaluated to allow for different sample numbers 
      .$coef_var           <- 0.1
      
      ### Define static variables 
      ###############################
      .$static$env    <- list(leaf.par=1000)
      .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # "pars" lists must contain parameter vectors that are of equal length,
      
      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$vars$fnames <- list(
        leaf.Alim   = c('f_lim_farquhar1980','f_lim_collatz1991'),
        leaf.etrans = c('f_j_farquharwong1984','f_j_collatz1991','f_j_harley1992')
      )

      if(eval_strings) {
        .$vars$pars_eval <- list(
          leaf.avn_25   = ' 10 * rnorm(n,1,.$coef_var)',
          leaf.bvn_25   = '  5 * rnorm(n,1,.$coef_var)',
          leaf.theta    = '0.9 * rnorm(n,1,.$coef_var)',
          leaf.e_ajv_25 = '0.9 * rnorm(n,1,.$coef_var)'
        )
      } else {
        n <- 2 * n
        .$vars$pars <- list(
          leaf.avn_25   =  10 * rnorm(n,1,.$coef_var),
          leaf.bvn_25   =   5 * rnorm(n,1,.$coef_var),
          leaf.theta    = 0.9 * rnorm(n,1,.$coef_var),
          leaf.e_ajv_25 = 0.9 * rnorm(n,1,.$coef_var)
        )
      }
      
      # .$vars$pars <- list(
      #   leaf.avn_25   = NA,
      #   leaf.bvn_25   = NA,
      #   leaf.theta    = NA,
      #   leaf.e_ajv_25 = NA
      # )
      # 
      # 
      # .$vars$pars_proc <- list(
      #   leaf.avn_25   = 'leaf.Alim',
      #   leaf.bvn_25   = 'leaf.Alim',
      #   leaf.theta    = 'leaf.etrans',
      #   leaf.e_ajv_25 = 'leaf.etrans'
      # )
      
      .$vars$env <- list(
        leaf.ca_conc  = c(400,600)
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
      out_list <- .$output_saltelli()
      # p1 <- xyplot(A~ci|leaf.etrans*leaf.rs,df,type='l',auto.key=T)
      # p1 <- bwplot(leaf.ca_conc~A|leaf.Alim*leaf.etrans,out_list[[1]],auto.key=T,abline=0)
      # c(out_list,p1)
      out_list
    }    
    
    ###########################################################################
    # end maat wrapper 
}) 











