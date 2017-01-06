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
      .$model$configure(func='write_fnames',df=data.frame(.$static$fnames,stringsAsFactors=F) )
      .$model$configure(func='write_pars',  df=data.frame(.$static$pars,stringsAsFactors=F)   )
      .$model$configure(func='write_env',   df=data.frame(.$static$env,stringsAsFactors=F)    )      
      
      # create dataframes of runtime variables  
      # expand the fnames and driving variables  
      .$dataf$fnames  <- if(!is.na(.$vars$fnames[1])) expand.grid(.$vars$fnames,stringsAsFactors=F) else NULL
      .$dataf$env     <- if(!is.na(.$vars$env[1]))    expand.grid(.$vars$env,stringsAsFactors=F   ) else NULL
      # - once Ye method has been generalised this can be revised and moved to the initial A/B selection loop
      if(.$wpars$UQ) {
        if(.$wpars$sobol) {
          .$dataf$fnames  <- if(!is.na(.$vars$fnames[1]))  expand.grid(c(.$vars$fnames,.$vars$fnamesB),stringsAsFactors=F) else NULL
        } else {
          .$dataf$fnamesB <- if(!is.na(.$vars$fnamesB[1])) expand.grid(.$vars$fnamesB,stringsAsFactors=F) else stop()
        }
      }
      # bind the parameter vectors OR expand pre-determined parameter vectors 
      # - this could be modified to allow a distribution function to be specifed
      if(.$wpars$UQ) { 
        # yet to solve is how to fix parameters to a relevant function (not 100% necessary, especially for MC runs). 
        # - however the ATS does this in a rigorous way allowing dynamic object construction.
        # - this will also be needed for generalising the Ye method to switch between process A and process(es) B
        .$dataf$pars    <- if(!is.na(.$vars$pars[1] ))   as.data.frame(do.call(cbind,.$vars$pars ))     else stop()
        # - once Ye method has been generalised this can be moved to the initial A/B selection loop
        .$dataf$parsB   <- if(!is.na(.$vars$parsB[1]))   as.data.frame(do.call(cbind,.$vars$parsB))     else stop()
        # this if for sobol is for sobol intended to be used in conjunction with a process rep analysis
        # - once Ye method has been generalised this can probably be removed
        if(.$wpars$sobol) .$dataf$pars  <- cbind(.$dataf$pars,.$dataf$parsB)
      } else {
        .$dataf$pars  <- if(!is.na(.$vars$pars[1]))   expand.grid(.$vars$pars,stringsAsFactors=F)  else NULL
      } 
      
      # add an extra column to the dataframes if they have only one column
      # - this prevents them being coerced to a vector in further functions
      if(!is.null(.$dataf$fnames)) if(dim(.$dataf$fnames)[2]==1) .$dataf$fnames <- data.frame(.$dataf$fnames,NA) 
      if(!is.null(.$dataf$pars))   if(dim(.$dataf$pars)[2]==1)   .$dataf$pars   <- data.frame(.$dataf$pars,NA) 
      if(!is.null(.$dataf$env))    if(dim(.$dataf$env)[2]==1)    .$dataf$env    <- data.frame(.$dataf$env,NA) 
      if(!is.null(.$dataf$met))    if(dim(.$dataf$met)[2]==1)    .$dataf$met    <- data.frame(.$dataf$met,NA)       
      
      # calculate dataframe lengths, if no dataframe return 1
      # - used to set the number of iterations in the run functions  
      # - parameter lengths 
      .$dataf$lf <- if(is.null(.$dataf$fnames)|(sum(dim(.$dataf$fnames)==0)==2)) 1 else length(.$dataf$fnames[,1]) 
      .$dataf$le <- if(is.null(.$dataf$env)|(sum(dim(.$dataf$env)==0)==2)      ) 1 else length(.$dataf$env[,1])    
      .$dataf$lm <- if(is.null(.$dataf$met)|(sum(dim(.$dataf$met)==0)==2)      ) 1 else length(.$dataf$met[,1])    
      if(.$wpars$UQ&!.$wpars$sobol) {
        .$dataf$lfB <- if(is.null(.$dataf$fnamesB)|(sum(dim(.$dataf$fnamesB)==0)==2)) 1 else length(.$dataf$fnamesB[,1]) 
        .$dataf$lp  <- .$wpars$n 
        .$dataf$lpB <- .$wpars$n^2 # not sure this line is needed anymore, the above line should be sufficient
      } else {
        .$dataf$lp  <- if(is.null(.$dataf$pars)|(sum(dim(.$dataf$pars)==0)==2) ) 1      else length(.$dataf$pars[,1])
      }
      
      # output summary of maat setup
      .$print_data()
      
      # run model
      .$print_run()
      
      .$dataf$out <-
        as.data.frame(
          do.call(
            rbind, {
              # multicore or not
              if(.$wpars$multic) {
                # process UQ run, or either general variable run or matrix A and B of Saltelli method
                if(.$wpars$UQ&!.$wpars$sobol) mclapply(1:.$dataf$lf, .$run_repA, mc.cores=.$wpars$procs)
                else                          mclapply(1:.$dataf$lf, .$runf,     mc.cores=.$wpars$procs)
              } else {
                # process UQ run, or either general variable run or matrix A and B of Saltelli method
                if(.$wpars$UQ&!.$wpars$sobol)   lapply(1:.$dataf$lf, .$run_repA)
                else                            lapply(1:.$dataf$lf, .$runf)
              }
            }
        ))
      
      # if Sobol analysis run saltellis output from parameter matrices ABi
      # output from parameter matrices A & B already generated by the above function
        if(.$wpars$sobol) {
          # prob should run output here and clear out df 
          print('Saltelli matrix AB completed', quote=F)
          print('',quote=F)      
          .$dataf$out_saltelli <-
            if(.$wpars$multic) mclapply(1:.$dataf$lf,.$runf_saltelli,mc.cores=.$wpars$procs)
            else                 lapply(1:.$dataf$lf,.$runf_saltelli)          
          print('Saltelli array Abi completed',quote=F)
          print('',quote=F)      
        } else NA             
      
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
#       data.frame(do.call(rbind,lapply(1:.$dataf$lp,.$runp)))
      data.frame(
        do.call(
          rbind, {
            if(.$wpars$multic) mclapply(1:.$dataf$lp,.$runp,mc.cores=floor(.$wpars$procs/.$dataf$lf))
            else                 lapply(1:.$dataf$lp,.$runp)
          }
      ))
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
    
    # Sobol al la Saltelli
    #################################
    
#     # Step 2 described in the WORD file to generate random parameter samples
#     # Initialise vectors
#     # Generate random parameter samples
#     X  <- rnorm(2*N) 
#     Y  <- rnorm(2*N)
#     # combine into a single matrix
#     pm <- cbind(X,Y)
#     
#     # generate function output vector & matrix
#     g  <- numeric(2*N) # this is the output vector from the two parameter matrices A & B
#     gp <- matrix(0,N,k) # where k is the number of parameters i.e. the same dims as matrices A, B, & AB
#     
#     for (i in 1:2*N ) {    
#       # Step 3 described in the WORD file to execute the model using the parameter samples
#       # Calculate model output g using parameters in row i
#       g[i] <- func(pm[i,])
#       
#       # Step 4 described in the WORD file to execute the model using the replaced parameter samples
#       if (i > N) {
#         # Caculate model output g using replaced parameters 
#         for(p in 1:k) { 
#           sub       <- rep(i-N,k)
#           sub[p]    <- i
#           smat      <- cbind(sub,1:k)
#           gp[i-N,k] <- func(pm[smat])
#         }
#       }
#     }
#     
#     list(g,gp)

    runf_saltelli <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnames[i,]),F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[i,])
      
      # call next run function
      if(.$wpars$multic) mclapply(1:.$dataf$le,.$rune_saltelli,mc.cores=floor(.$wpars$procs/.$dataf$lf))
      else                 lapply(1:.$dataf$le,.$rune_saltelli)
    }

    rune_saltelli <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
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
      # returns a numeric matrix
      do.call(rbind,lapply((.$wpars$n+1):.$dataf$lp,.$runp_saltelli,pk=p))      
    }
    
    runp_saltelli <- function(.,j,pk) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # create index matrix to create row on matrix ABi for the .$dataf$par matrix (which is matrix A stacked on top of matrix B)
      sub     <- rep(j-.$wpars$n,dim(.$dataf$pars)[2])
      sub[pk] <- j
      smat    <- cbind(sub,1:dim(.$dataf$pars)[2])      
      
      # create a dataframe from vector and add names
      psdf        <- data.frame(t(.$dataf$pars[smat]))
      names(psdf) <- names(.$dataf$pars)
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=psdf,F)
      if(.$wpars$cverbose) .$printc('pars',psdf)
      
#       # call environment run function
#       data.frame(do.call(rbind,lapply(1:.$dataf$le,.$rune)))      

#       print('par here')
      
      # run model
      # - this requires a numeric vector for the output
      # - warnings are suppressed because character strings are converted to NAs which comes with a warning   
      suppressWarnings(as.numeric( .$model$run() ))        
    }
    
    
    
    ###########################################################################
    # nested run functions for Process Sensitivity Analysis after Ye etal 201? - 2 process case 
    
    # This case will be specific to the core enzyme kinetic model of 6 equations and a solver as described in the documentation
    # In reality the two processes are the overall model vs electron transport models (including the Jmax~Vcmax relationship)

    # The Ye method is a nested system of loops, two loops for each process
    # - the first loop over each process representation of process A, the second the parameter sampling loop common to MC methods
    # - for this test case process A is the overall model, process B the electron transport model
    # - therefore the initial loop for process A consists of a single iteration so is not needed as a loop - function names will be set during init 
    
    # below is the loop structure:
    # Loop 1: parameter loop for process A             - use standard location for variable parameter values
    # Loop 2: process representation loop for proces B - use standard location for variable function values
    # Loop 3: parameter loop for proces B              - use an additional non-standard location for variable parameter values

    run_repA <- function(.,g) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnames[g,]),F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[g,])
      
      # call process A parameter run function
      # data.frame(do.call(rbind,lapply(1:.$dataf$lp,.$run_parA,offset=g)))      
      data.frame(
        do.call(rbind,
                if(.$wpars$multic) mclapply(1:.$dataf$lp,.$run_parA,offset=g,mc.cores=floor(.$wpars$procs/.$dataf$lf))
                else                 lapply(1:.$dataf$lp,.$run_parA,offset=g)
      ))
    }
        
    run_parA <- function(.,h,offset) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(func='write_pars',df=data.frame(.$dataf$pars[h,]),F)
      if(.$wpars$cverbose) .$printc('pars',.$dataf$pars[h,])
      
      # call process B process representation run function
      # os is the row in the factorial ftypes and pars matrix
      os <- .$wpars$n*(offset - 1) + h     
      data.frame(do.call(rbind,lapply(1:.$dataf$lfB,.$run_repB,offset=os)))      
    }

    run_repB <- function(.,i,offset) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure function names in the model
      if(!is.null(.$dataf$fnamesB)) .$model$configure(func='write_fnames',df=data.frame(.$dataf$fnamesB[i,]),F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnamesB[i,])
      
      # call process B parameter run function
      # os is the row in the parsB matrix
      os <- .$dataf$lfB*(offset - 1) + i    
      data.frame(do.call(rbind,lapply(os:(os+.$wpars$n-1),.$run_parB)))      
    }
    
    run_parB <- function(.,j) {
      # This wrapper function is called from an lapply or mclappy function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are independent and non-sequential
      
      # configure parameters in the model
      if(!is.null(.$dataf$parsB)) .$model$configure(func='write_pars',df=data.frame(.$dataf$parsB[j,]),F)
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
      fnames  = NA,
      fnamesB = NA,
      pars    = NA,
      parsB   = NA,
      env     = NA
    )
    
    # input/output data frames
    # currently dataframes with an associated length
    dataf  <- list( 
      fnames  = NULL,
      lf      = NULL,
      fnamesB = NULL,
      lfB     = NULL,
      pars    = NULL,
      lp      = NULL,
      parsB   = NULL,
      lpB     = NULL,
      env     = NULL,
      le      = NULL,
      met     = NULL,          # a dataframe of sequential meteorological driving data, for running the analysis at a particular site for example 
      lm      = NULL,
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
      n        = numeric(1),  # emsemble number in parameteric UQ
      sobol    = F,           # run a Sobol parameter sensitivity analysis
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
#           vpars <- if(is.null(.$dataf$pars))  NULL  else if(.$wpars$UQ) .$vars$pars[1] else .$vars$pars
          vpars    <- if(is.null(.$dataf$pars))    NULL else .$vars$pars[[1]]
          venv     <- if(is.null(.$dataf$env))     NULL else .$vars$env
          vfnames  <- if(is.null(.$dataf$fnames))  NULL else .$vars$fnames
          vfnamesB <- if(is.null(.$dataf$fnamesB)) NULL else .$vars$fnamesB
#           vparsB  <- if(.$wpars$UQ) .$vars$parsB[[1]]  else NULL
          
          # if Ye (i.e. process) UQ - need to write a function to do this
          # the number of rows in the resultant dataframe should be lf*lfB*n^2*le
          if(.$wpars$UQ&!.$wpars$sobol) {
            # pars matrix
            vpars    <- matrix(unlist(.$vars$pars),.$wpars$n) 
            # parsB matrix
            vparsB   <- matrix(unlist(.$vars$parsB),.$dataf$lf*.$dataf$lfB*.$wpars$n^2)
            
            # combine input into a single dataframe in order of output dataframe,
            #  - i.e. repeats lines in input dataframes/matrices to align with output
            # vardf    <- suppressWarnings(cbind(
            vardf    <- cbind(
              # fnames of process A - length lf * lfB * le * n^2
              rep(unlist(vfnames),each=.$dataf$lfB*.$dataf$le*.$wpars$n^2),
              # pars
              apply(vpars,2,function(v) rep(rep(v,each=.$dataf$lfB*.$dataf$le*.$wpars$n),.$dataf$lf) ),
              # fnames of process B
              # rep(rep(unlist(vfnames),each=.$dataf$le*.$wpars$n),.$wpars$n) ,
              rep(rep(unlist(vfnamesB),each=.$dataf$le*.$wpars$n),.$dataf$lf*.$wpars$n),
              # parsB - this will need to be broadened to a matrix when this function is expanded to > 2 processes
              apply(vparsB,2,function(v) rep(v,each=.$dataf$le) ),
              # environment
              unlist(venv) 
            )#)
            
            # add names to dataframe
            vardf        <- as.data.frame(vardf)
            names(vardf) <- c(names(vfnames),names(.$vars$pars),names(vfnamesB),names(.$vars$parsB),names(venv))
          
          } else if(.$wpars$sobol) { 
            vpars   <- .$dataf$pars
            vpars   <- as.data.frame(apply(vpars,2,function(v) rep(v,each=.$dataf$le) ))
            vardf   <- suppressWarnings(cbind( 
              rep(unlist(vfnames),each=.$dataf$le*2*.$wpars$n),
              vpars,
              unlist(venv) 
            ))
            vardf        <- as.data.frame(vardf)
            names(vardf) <- c(names(vfnames),names(.$dataf$pars),names(venv))
            
          } else 
            vardf <- expand.grid(c(venv,vpars,vfnames),stringsAsFactors=F)        
          
          # if no met data
          if(is.null(.$dataf$met)) {
            suppressWarnings(cbind(vardf,.$dataf$out))

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
    
    output_saltelli <- function(.) {
      # creates output for a saltelli Sobol sensitivity analysis 
      # A and B matrices are stacked in a single matrix, which for each model and environment combination are then stored in an array 

      # AB output is an array
      # - dim 1 (rows)   model combination
      # - dim 2 (cols)   environment combination
      # - dim 3 (slices) sample
      # - dim 4          output variable
      
      # ABi output is a nested list
      # - list 1 model combination
      # - list 2 environment combination
      # - list 2 is composed of 3 dimesional output arrays 
      #    - dim 1 (rows)   sample
      #    - dim 2 (cols)   output variable
      #    - dim 3 (slices) parameter that has used value from matrix B while all other par values are from matrix A
      
      # create AB output matrix array
      # - suppressWarnings on the as.numeric call as there can be character strings in .$dataf$out which cause warnings when converted to NAs
      ab_mat_out  <- array(suppressWarnings(as.numeric(unlist(.$dataf$out))),c(.$dataf$le,2*.$wpars$n,.$dataf$lf,length(.$dataf$out)))
      tab_mat_out <- aperm(ab_mat_out,c(3,1,2,4))
      
#       list(maat_out=.$output(),AB=tab_mat_out,ABi=.$dataf$out_saltelli)
      list(AB=tab_mat_out,ABi=.$dataf$out_saltelli)
    }
    
    print_data <- function(.) {
      print('',quote=F)
      print("MAAT :: summary of data",quote=F)
      print('',quote=F)

      ens_n <- 
        if(.$wpars$UQ) {
          if(.$wpars$sobol) .$dataf$lf*.$wpars$n*(2+dim(.$dataf$pars)[2])*.$dataf$le
          else              .$dataf$lf*.$dataf$lfB*.$wpars$n^2*.$dataf$le
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

    # test function for Ye method Sobol process sensitivity analysis
    .test_ye <- function(.,metd=F,mc=T,pr=4,oconf=F) {
      
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
      .$wpars$UQ           <- T   # run a UQ/SA style ensemble 
      .$wpars$sobol        <- F   # Ye style SA ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions) 
      .$wpars$n            <- 3   # number of parameter samples in each loop 
      
      ### Define static variables 
      ###############################
      .$model$env$par <- 1000
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # "pars" lists must contain parameter vectors that are of equal length,

      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$vars$fnames <- list(
        leaf.Alim   = c('f_lim_farquhar1980','f_lim_collatz1991')
      )

      .$vars$fnamesB <- list(
        leaf.etrans = c('f_j_farquhar1980','f_j_collatz1991')
      )
      
      .$model$fnames$vcmax <- 'f_vcmax_lin'
      .$vars$pars <- list(
        leaf.avn_25 = 9:11,
        leaf.bvn_25 = 4:6
      )

      coef_var <- 0.1
      n        <- 4 * 9 # lf * lfB * n^2
      .$vars$parsB <- list(
        theta       = 0.9 * rnorm(n,1,coef_var),
        e_ajv_25    = 0.9 * rnorm(n,1,coef_var)
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
    
    ###########################################################################
    # end maat wrapper 
}) 











