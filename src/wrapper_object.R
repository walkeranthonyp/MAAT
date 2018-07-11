################################
#
# Wrapper object for model objects 
# 
# AWalker December 2015
#
################################

library(proto)
library(parallel)

source('functions/general_functions.R')
source('functions/calc_functions.R')



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
    build <- function(., model, ... ) {
      .$model <- as.proto(get(model)$as.list())
      .$model$build(...)
    }
    
    # clean function to reset object
    clean <- function(.) {
      for(d in 1:length(.$dataf)) if(!is.null(.$dataf[[d]])) .$dataf[[d]] <- NA
      gc()
    }
    
    
    ###########################################################################
    # Main run function

    run   <- function(.,verbose=T) {
      
      # Initialise 
      if(!.$wpars$unit_testing) {
        if(.$wpars$UQtype=='ye'|.$wpars$UQtype=='mcmc') .$wpars$eval_strings <- T
        .$init()
      }

      # Initialisation checks
      # need to add a check for equal par vector lengths if this is a UQ run and not eval_strings
      # for Ye et al SA method
      # - due to different parameter sample numbers in process A and B loops,
      # - parameters samples must be generated from code snippets as strings
      if(.$wpars$eval_string&is.null(.$dynamic$pars_eval)) {
        stop(paste('wrapper: eval_strings = T but dynamic$pars_eval not set. \n
              vars$pars_eval:,\n',.$dynamic$pars_eval,'\n 
              NOTE: Ye method SA must draw parameter samples during runtime \n
              from code snippets written as strings in dynamic$pars_eval')) 
      }

      # configure initialisation lists 
      ########################################

      # initialise model with static variables
      if(!is.null(.$static$fnames)) .$model$configure(vlist='fnames', df=t(as.matrix(.$static$fnames,stringsAsFactors=F))[1,] ) else if(!.$wpars$unit_testing) stop('Static fnames not defined')
      if(!is.null(.$static$pars))   .$model$configure(vlist='pars',   df=t(as.matrix(.$static$pars,stringsAsFactors=F))[1,]   ) else if(!.$wpars$unit_testing) stop('Static pars not defined')
      if(!is.null(.$static$env))    .$model$configure(vlist='env',    df=t(as.matrix(.$static$env,stringsAsFactors=F))[1,]    ) else if(!.$wpars$unit_testing) stop('Static env not defined')      

      # create matrices of runtime variables  
      ######################################
      # expand the fnames and driving variables lists into matrices 
      .$dataf$fnames  <- if(!is.null(.$dynamic$fnames)) as.matrix(expand.grid(.$dynamic$fnames,stringsAsFactors=F)) else NULL
      .$dataf$env     <- if(!is.null(.$dynamic$env))    as.matrix(expand.grid(.$dynamic$env,stringsAsFactors=F   )) else NULL
      
      # if an SA/UQ run
      if(.$wpars$UQ) {
          
        # if Saltelli style Sobol
        if(.$wpars$UQtype=='saltelli') {
          
          if(is.null(.$dynamic$pars)) { 
            if(!is.null(.$dynamic$pars_eval)) { 
              # increase parameter sample number
              .$wpars$n <- .$wpars$n * .$wpars$nmult
              # sample parameters from character string code snippets to generate matrices A and B
              n <- 2 * .$wpars$n
              .$dynamic$pars <- lapply(.$dynamic$pars_eval,function(cs) eval(parse(text=cs)))
            } else  stop('wrapper: pars (or pars_eval) list in vars list is empty')
          }         

          # create pars matrix
          .$dataf$pars  <- do.call(cbind, .$dynamic$pars )
          
          # remove potentially large pars list 
          .$dynamic$pars   <- lapply(.$dynamic$pars, function(e) numeric(1) )        

        # Ye et al process SA method 
        } else if(.$wpars$UQtype=='ye') {
          
          # need a minimum of >1 processes
          if(dim(.$dataf$fnames)[2]<=1) stop('need more than one process for a process sesitivity analysis')
          
          # check input dynamic$pars* are same length
          test_in <- length(.$dynamic$pars_eval) - length(.$dynamic$pars_proc)
          if(test_in!=0) stop('wrapper: Parameter input vectors - pars_eval & pars_proc - are not the same length')
          
          # assign same list structure as dynamic$pars_eval to vars$pars 
          .$dynamic$pars   <- lapply(.$dynamic$pars_eval,function(e) numeric(1) )
          
          # check input vars$pars* elements have same names
          # - to be done
         
        # if MCMC 
        if(.$wpars$UQtype=='mcmc') {
          
          # sample parameters from character string code snippets to generate initial proposal from priors 
          n <- .$wpars$mcmc_chains
          .$dynamic$pars <- lapply(.$dynamic$pars_eval,function(cs) eval(parse(text=cs)))

          # create pars / proposal matrix 
          .$dataf$pars  <- do.call(cbind, .$dynamic$pars )
         
          # create accepted proposal array 
          .$dataf$pars_array    <- array(1, dim=c(dim(.$dataf$pars),.$wpars$mcmc_maxiter) )
 
          # create accepted proposal likelihood matrix 
          .$dataf$pars_lklihood <- matrix(1, .$wpars$mcmc_chains, .$wpars$mcmc_maxiter ) 
 
          # remove initialisation pars list 
          .$dynamic$pars        <- lapply(.$dynamic$pars, function(e) numeric(1) )        

        } else {
          stop(paste('wrapper: no method for SA/UQ type',.$wpars$UQtype))
        }
        
      } else {
        # not a formal SA/UQ run - factorial combination of variables specified in the vars lists  
        .$dataf$pars <- if(!is.null(.$dynamic$pars)) as.matrix(expand.grid(.$dynamic$pars,stringsAsFactors=F))   else NULL
      } 
      
      # calculate input matrix lengths 
      ################################
      # if no matrix return 1
      # - used to set the number of iterations in the run functions  
      if(.$wpars$UQ&.$wpars$UQtype=='ye') {
        # determine number of processes to be analaysed
        .$dataf$lf <- length(.$dynamic$fnames)
        
      } else {
        # any type of run other than Ye process sensitivity analysis 
        .$dataf$lf <- if(is.null(.$dataf$fnames)) 1 else length(.$dataf$fnames[,1]) 
        .$dataf$lp <- if(is.null(.$dataf$pars))   1 else length(.$dataf$pars[,1])
        
      }
      
      # enviroment matrix and met matrix
      .$dataf$le <- if(is.null(.$dataf$env)) 1 else length(.$dataf$env[,1])    
      .$dataf$lm <- if(is.null(.$dataf$met)) 1 else length(.$dataf$met[,1])    
      
      # store model output template (currently must be a vector)
      # - this will fail when output is a vector of variable length depending on parameter values
      .$dataf$mout <- .$model$output()
      
      # print summary of maat setup
      .$print_data()
      .$print_data(otype='run')


      # run model ensemble
      # call initial run function in the hierarchy of nested run functions
      ################################
      
      # process UQ run
      if(.$wpars$UQ&.$wpars$UQtype=='ye') {
        # - Ye process SA is not multicored at this stage as multicoring here messes with the processes A and B in the data structure  
        vapply(1:.$dataf$lf, .$run_general_process_SA, numeric(0) )
        
      # MCMC run
      } else if(.$wpars$UQ&.$wpars$UQtype=='mcmc') {
         
        # if an MCMC and no met dataset has been specified, stop
        if(is.null(.$dataf$met)&.$wpars$UQtype=='mcmc') stop('No current method to run MCMC without a met dataset')
        if(length(.$dataf$mout)!=1)                     stop('No current method to run MCMC with multiple model outputs')

        # initialise output matrix
        .$dataf$out <- matrix(0, .$dataf$lp, length(.$dataf$met) )

        # call run function
        if(.$wpars$multic) mclapply( 1:.$dataf$lf, .$runf_mcmc, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lp)), mc.preschedule=F )
        else                 lapply( 1:.$dataf$lf, .$runf_mcmc )
       
        # print summary of results
        .$print_output()
         
      # Saltelli SA matrix AB run or factorial run
      } else {
        
        # if a Saltelli SA and a met dataset has been specified, stop
        if(!is.null(.$dataf$met)&.$wpars$UQtype=='saltelli') stop('No current method to run Saltelli SA with a met dataset')

        # initialise output matrix
        .$dataf$out <- matrix(0, .$dataf$lm*.$dataf$le*.$dataf$lp*.$dataf$lf, length(.$dataf$mout) )
        colnames(.$dataf$out) <- names(.$dataf$mout)

        # call run function
        .$dataf$out[] <-
          do.call( 'rbind', {
              if(.$wpars$multic) mclapply( 1:.$dataf$lf, .$runf, mc.cores=min(.$dataf$lf,.$wpars$procs), mc.preschedule=F )
              else                 lapply( 1:.$dataf$lf, .$runf )
          })
       
        # print summary of results
        .$print_output()
      }
        
      # if Saltelli SA, write output from AB matrix iteration above, then run and write ABi matrix iteration
      if(.$wpars$UQtype=='saltelli') {
        
        # write AB output array
        if(.$wpars$unit_testing) { hd <- getwd(); setwd('~/tmp'); ofname <- 'Salt_test' } else setwd(odir)
        write_to_file(.$output_saltelli_AB(), paste(ofname,'salt','AB',sep='_'), type='rds' )  

        # write dataf matrices used in AB run
        write_to_file(list(fnames=.$dataf$fnames, pars=.$dataf$pars, env=.$dataf$env ), paste(ofname,'salt','AB','dataf',sep='_'), type='rds' )  
        
        # remove large out array
        if(!.$wpars$unit_testing) .$dataf$out <- NULL   
        .$print_saltelli()
        
        # initialise output array
        # - dim 1 (rows)      output variable
        # - dim 2 (columns)   sample
        # - dim 3 (slices)    parameter that has used value from matrix B while all other par values are from matrix A 
        # - dim 4 (cube rows) environment combination
        # - dim 5 (cube cols) model combination
        .$dataf$out_saltelli <- array(0, dim=c(length(.$dataf$mout), .$wpars$n, dim(.$dataf$pars)[2], .$dataf$le, .$dataf$lf ))
        dimnames(.$dataf$out_saltelli) <- list(names(.$dataf$mout), NULL, colnames(.$dataf$pars), NULL, apply(.$dataf$fnames, 1, toString) )
        
        # run over ABi matrices
        .$dataf$out_saltelli[] <- vapply(1:.$dataf$lf, .$runf_saltelli, .$dataf$out_saltelli[,,,,1] )          

        # write ABi output array 
        write_to_file(.$output_saltelli_ABi(), paste(ofname,'salt','ABi',sep='_'), type='rds' )  
        if(!.$wpars$unit_testing) .$dataf$out_saltelli <- matrix(1)   
        print(paste('Saltelli array ABi completed',Sys.time()),quote=F)
        print('',quote=F)
        
        if(.$wpars$unit_testing) setwd(hd)
      }             
    }
    
    
    ###########################################################################
    # nested run functions - 
    ###########################################################################

    printc <- function(.,r1,r2) {
      print(r1,quote=F)
      print(r2,quote=F)
    }
    
    # takes a >=3 D array and stacks it into a 2D matrix
    stack <- function(., a ) {
      apply(a, 2, function(v) v )
    }


    # for factorial runs
    ###########################################################################
    
    runf <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call runp

#      if(i == 1) {
#        print('runf',quote=F)
#        print(.$dataf$fnames,quote=F)
#        print(.$dataf$pars,quote=F)
#        print(.$dataf$env,quote=F)
#      }
      
      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[i,], F )
      if(.$wpars$cverbose)         .$printc('fnames', .$dataf$fnames[i,] )

      # call next run function
      do.call( 'rbind', {
          if(.$wpars$multic) mclapply(1:.$dataf$lp, .$runp, fi=i, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lf)), mc.preschedule=T  )
          else                 lapply(1:.$dataf$lp, .$runp, fi=i )
      })
    }
    
    runp <- function(.,j,fi) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
      # assumes that each row of the pars matrix are independent and non-sequential
      # call rune
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[j,], F )
      if(.$wpars$cverbose)       .$printc('pars', .$dataf$pars[j,] )
      
      # call next run function
      funv   <- if(is.null(.$dataf$met)) .$dataf$mout else array(0, dim=c(.$dataf$lm, length(.$dataf$mout) ) )
      out    <- vapply(1:.$dataf$le, .$rune, funv )

      # out has the potential to be a vector, matrix (needs transposed), or an array (needs stacking)
      # returns matrix
      if(class(out)=='array') .$stack(out) else t(out)
    }
    
    rune <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$env matrix to the model
      # assumes that each row of the env matrix are independent and non-sequential
      # call .$model$run or .$model$run_met if met data are provided
      
      # configure environment in the model
      if(!is.null(.$dataf$env)) .$model$configure(vlist='env', df=.$dataf$env[k,], F )
      if(.$wpars$cverbose)      .$printc('env', .$dataf$env[k,] )
      
      # call next run function
      if(is.null(.$dataf$met)) {
        .$model$run()        
      } else {
        t(vapply(1:.$dataf$lm, .$model$run_met, .$dataf$mout ))
      }  
    }

    
    # for MCMC runs
    ###########################################################################
    
    runf_mcmc <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call run_mcmc

      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[i,], F )
      if(.$wpars$cverbose)         .$printc('fnames', .$dataf$fnames[i,] )

      # evaluate model over initial proposals derived from prior
      .$dataf$out[]  <- 
        do.call( 'rbind', {
            if(.$wpars$multic) mclapply(1:.$dataf$lp, .$runp_mcmc, mc.cores=min(.$wpars$procs,.$dataf$lp), mc.preschedule=T  )
            else                 lapply(1:.$dataf$lp, .$runp_mcmc )
        })

      # calculate likelihood of initial proposal
      .$dataf$pars_lklihood[,1] <- .$proposal_lklihood()   

      # run MCMC 
      vapply(1:(.$wpars$mcmc_maxiter-1), .$run_mcmc, numeric(0) )
    }
    
    run_mcmc <- function(.,j) {
      # This wrapper function is called from a vapply function to iterate / steop chains in an MCMC
      # runs in serial as each step depends on the previous step
      # call runp_mcmc
    
      # generate proposal matrix
      .$dataf$pars[] <- .$gen_proposal_demc()   
    
      # evaluate model for proposal on each chain
      .$dataf$out[]  <- 
        do.call( 'rbind', {
            if(.$wpars$multic) mclapply(1:.$dataf$lp, .$runp_mcmc, mc.cores=min(.$wpars$procs,.$dataf$lp), mc.preschedule=F  )
            else                 lapply(1:.$dataf$lp, .$runp_mcmc )
        })
    
      # calculate likelihood of proposals on each chain
      lklihood <- .$proposal_lklihood()   
      
      # accept / reject proposals on each chain 
      accept <- .$proposal_accept()

      # update accepted proposal array (.$dataf$pars_array) and likelihood matrix (.$dataf$pars_lklihood) 
      .$dataf$pars_array[,,j+1]   <- if(accept) .$dataf$pars else .$dataf$pars_array[,,j]    
      .$dataf$pars_lklihood[,j+1] <- if(accept) lklihood     else .$dataf$pars_lklihood[,j]    

      # test for convergence every x iterations
 
      # return nothing - this is not part of the MCMC, allows use of the more stable vapply to call this function   
      numeric(0) 
    }
 
    runp_mcmc <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
      # runs each chain at each iteration in MCMC
      # assumes that each row of the pars matrix are independent and non-sequential
      # call run met
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[k,], F )
      if(.$wpars$cverbose)       .$printc('pars', .$dataf$pars[k,] )
      
      # call metdata run function
      #funv   <- if(is.null(.$dataf$met)) .$dataf$mout else array(0, dim=c(.$dataf$lm, length(.$dataf$mout) ) )
      vapply(1:.$dataf$lm, .$model$run_met, .$dataf$mout )
    }
   
    # generate proposal using DE-MC algorithm  
    gen_proposal_demc <- function(.) {
      
    }   
 
    # calculate proposal likelihood using ...  
    proposal_lklihood <- function(.) {
      
    }   
 
    # calculate proposal acceptance using ...  
    proposal_accept <- function(.) {
      
    }   
 
    # calculate convergence using ...  
    chain_convergence <- function(.) {
      
    }   
 
 
    # for ABi array for Sobol SA using Saltelli method
    ###########################################################################
    
    runf_saltelli <- function(.,i) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call rune_saltelli
      
      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames',df=.$dataf$fnames[i,],F)
      if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[i,])
      
      # call next run function
      vapply({
        if(.$wpars$multic) mclapply(1:.$dataf$le, .$rune_saltelli, mc.cores=min(.$dataf$le,.$wpars$procs), mc.preschedule=F ) 
        else                 lapply(1:.$dataf$le, .$rune_saltelli )
      }, function(a) a, .$dataf$out_saltelli[,,,1,1] )
    }

    rune_saltelli <- function(.,k) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$env matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call runpmat_saltelli
      
      # configure environment in the model
      if(!is.null(.$dataf$env)) .$model$configure(vlist='env', df=.$dataf$env[k,], F )
      if(.$wpars$cverbose)      .$printc('env', .$dataf$env[k,] )
      
      # call parameter matrix run function
      if(is.null(.$dataf$met)){
        
        # call next run function, wrapped within vapply to convert (mc)lapply list output to an array
        # returns a numeric array - model output variable (rows), sample (columns), parameter (slices)
        vapply({
          if(.$wpars$multic) mclapply(1:dim(.$dataf$pars)[2], .$runpmat_saltelli, mc.cores=max(1,floor(.$wpars$procs/.$dataf$le)), mc.preschedule=T ) 
          else                 lapply(1:dim(.$dataf$pars)[2], .$runpmat_saltelli )
        },function(a) a, .$dataf$out_saltelli[,,1,1,1] )
        
      } else {
        # met data run not yet supported with Sobol, but should be caught before getting here
        stop('Saltelli')
      }  
    }
    
    runpmat_saltelli <- function(.,p) {
      # This wrapper function is called from an lapply or mclappy function to be run once for each parameter (i.e. each column of the dataf$pars matrix)
      # call runp_saltelli

      # returns a numeric matrix - model output variable (rows), sample (columns)
      vapply((.$wpars$n+1):.$dataf$lp, .$runp_saltelli, .$dataf$mout, pk=p )
    }
    
    runp_saltelli <- function(.,j,pk) {
      # This wrapper function is called from an lapply or mclappy function 
      # wrapper subscripts parameter matrix AB to give the row on matrix ABi
      # assumes that each row of the matrix are independent and non-sequential
      # call .$model$run
      
      # create index matrix to create row on matrix ABi for the .$dataf$par matrix (which is matrix A stacked on top of matrix B)
      sub     <- rep(j-.$wpars$n, dim(.$dataf$pars)[2] )
      sub[pk] <- j
      smat    <- cbind(sub, 1:dim(.$dataf$pars)[2] )      
      
      # create a matrix from vector and add names
      psdf        <- t(.$dataf$pars[smat])
      names(psdf) <- colnames(.$dataf$pars)

      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=psdf, F )
      if(.$wpars$cverbose) .$printc('pars', psdf )
      
      # run model
      .$model$run()        
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
      # then creates the fnames and pars matrices for process A and process B,
      # calls the below run function,
      # outputs an .RDS for each process segregation

      # create the fnames matrices for process A and process B
      .$dataf$fnames  <- if(!is.na(.$dynamic$fnames[f]))       as.matrix(expand.grid(.$dynamic$fnames[f] ,stringsAsFactors=F)) else stop()
      .$dataf$fnamesB <- if(!any(is.na(.$dynamic$fnames[-f]))) as.matrix(expand.grid(.$dynamic$fnames[-f],stringsAsFactors=F)) else stop()
      
      # determine the number of the rows in the fnames process matrices
      .$dataf$lfA     <- if(is.null(.$dataf$fnames )) 1 else length(.$dataf$fnames[,1]) 
      .$dataf$lfB     <- if(is.null(.$dataf$fnamesB)) 1 else length(.$dataf$fnamesB[,1]) 
      
      # partition the parameters to process A and and process B
      .$procA_name    <- names(.$dynamic$fnames)[f]
      .$procA_subs    <- which(unlist(.$dynamic$pars_proc)==.$procA_name)

      # evaluate parameter strings to sample vectors
      # - this allows a distribution function to be specifed
      # - also allows the dynamic calcuation of n for process A and B parameter samples
      if(.$wpars$eval_strings) {
        n <- .$wpars$n
        .$dynamic$pars[.$procA_subs ] <- lapply(.$dynamic$pars_eval[.$procA_subs ],function(cs) eval(parse(text=cs)))
        n <- .$dataf$lfA * .$dataf$lfB * .$wpars$n^2
        .$dynamic$pars[-.$procA_subs] <- lapply(.$dynamic$pars_eval[-.$procA_subs],function(cs) eval(parse(text=cs)))
      }
      
      # bind the parameter vectors into run matrices 
      .$dataf$pars    <- if(!is.na(.$dynamic$pars[1])) do.call(cbind,.$dynamic$pars[.$procA_subs] ) else stop()
      .$dataf$parsB   <- if(!is.na(.$dynamic$pars[2])) do.call(cbind,.$dynamic$pars[-.$procA_subs]) else stop()
      .$dynamic$pars  <- lapply(.$dynamic$pars_eval,function(e) numeric(1) ) 

      # determine the number of the rows in parameter matrices
      .$dataf$lp      <- .$wpars$n # convert these to be the row number of the actual matrices
      .$dataf$lpB     <- .$dataf$lfA * .$dataf$lfB * .$wpars$n^2
      
      # initialise output array
      # - dim 1 (rows)        output variable
      # - dim 2 (columns)     environment combination
      # - dim 3 (slices)      process(es) B parameter sample 
      # - dim 4 (cube rows)   process(es) B representation(s)
      # - dim 5 (cube cols)   process A parameter sample
      # - dim 6 (cube slices) process A representation
      
      # if met data then ... .$dataf$out     <- array(0, c(length(.$dataf$mout), .$dataf$lm, .$dataf$le, .$wpars$n, .$dataf$lfB, .$wpars$n, .$dataf$lfA  ) )
      .$dataf$out           <- array(0, c(length(.$dataf$mout), .$dataf$le, .$wpars$n, .$dataf$lfB, .$wpars$n, .$dataf$lfA  ) )
      dimnames(.$dataf$out) <- list(names(.$dataf$mout), NULL, NULL, apply(.$dataf$fnamesB, 1, toString), NULL, .$dataf$fnames )
      
      print('',quote=F)
      print(paste('started process:', colnames(.$dataf$fnames), Sys.time()), quote=F )

      # call the below run function
      .$dataf$out[] <- vapply({
        if(.$wpars$multic) mclapply(1:.$dataf$lfA, .$run_repA, mc.cores=min(.$dataf$lfA,.$wpars$procs), mc.preschedule=F )
        else                 lapply(1:.$dataf$lfA, .$run_repA)
      }, function(a) a, .$dataf$out[,,,,,1] )
      
      # process & record output
      if(.$wpars$unit_testing) { hd <- getwd(); setwd('~/tmp'); ofname <- 'Ye_test' }
      else                     setwd(odir)
      write_to_file(.$dataf$out, paste(ofname, 'proc', f, sep='_' ), type='rds' )
      
      print(paste('completed process:', colnames(.$dataf$fnames), Sys.time() ), quote=F )
      if(.$wpars$unit_testing) setwd(hd)
      
      # return nothing
      numeric(0)
    }
    
    run_repA <- function(.,g) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
      # assumes that each row of the fnames matrix are independent and non-sequential
      # call run_parA
     
      print(paste('started representation:', .$dataf$fnames[g,], ', of process:', colnames(.$dataf$fnames)), quote=F )
 
      # configure function names in the model
      if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[g,] , F )
      if(.$wpars$cverbose) .$printc('fnames', .$dataf$fnames[g,] )

      # calculate offset to correctly subset parsB matrix
      osg <- .$wpars$n * .$dataf$lfB * (g-1) 
 
      # call process A parameter run function
      vapply({
        if(.$wpars$multic) mclapply(1:.$dataf$lp, .$run_parA, offset=osg, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lfA)), mc.preschedule=F  )
        else                 lapply(1:.$dataf$lp, .$run_parA, offset=osg)
      }, function(a) a, .$dataf$out[,,,,1,1] )
    }
        
    run_parA <- function(.,h,offset) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
      # assumes that each row of the pars matrix are independent and non-sequential
      # call run_repB
      
      # configure parameters in the model
      if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[h,], F )
      if(.$wpars$cverbose) .$printc('pars', .$dataf$pars[h,] )
      
      # calculate offset to correctly subset parsB matrix
      osh  <- offset + .$dataf$lfB * (h-1)        
      
      # call process B process representation run function
      vapply(1:.$dataf$lfB, .$run_repB, .$dataf$out[,,,1,1,1], offset=osh )
    }

    run_repB <- function(.,i,offset) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnamesB matrix to the model
      # assumes that each row of the fnamesB matrix are independent and non-sequential
      # call run_parB
      
      # configure function names in the model
      if(!is.null(.$dataf$fnamesB)) .$model$configure(vlist='fnames', df=.$dataf$fnamesB[i,], F )
      if(.$wpars$cverbose) .$printc('fnames', .$dataf$fnamesB[i,] )
      
      # calculate offset to correctly subset parsB matrix
      os  <- offset + i
      # oss is a vector of the row subscripts for the parsB matrix
      oss <- (.$wpars$n*(os-1) + 1):(.$wpars$n*(os))    
      
      # call process B parameter run function
      vapply(oss, .$run_parB, .$dataf$out[,,1,1,1,1] )
    }
    
    run_parB <- function(.,j) {
      # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$parsB matrix to the model
      # assumes that each row of the parsB matrix are independent and non-sequential
      # call rune
      
      # configure parameters in the model
      if(!is.null(.$dataf$parsB)) .$model$configure(vlist='pars', df=.$dataf$parsB[j,], F )
      if(.$wpars$cverbose)        .$printc('pars', .$dataf$parsB[j,] )
      
      # call the environment run function
      vapply(1:.$dataf$le, .$rune, .$dataf$out[,1,1,1,1,1] )
      
      # ignoring met data possibility for now 
      # # out has the potential to be a vector, matrix (needs transposed), or an array (needs stacking)
      # # returns matrix
      # if(class(out)=='array') .$stack(out) else t(out)
    }
        
    
    ###########################################################################
    # initialisation function

    # prefix name of model object to variable names
    # this flattens the object|variable hierarchy in the list structure
    # allowing single run matrices that contain variables for multiple model objects
    # each line of the matrix is passed to the configure function in the model object
    # the model object name in the variable name allows the configure function to correctly parse the variable.        
    init <- function(.) {

      # combines sublists in list l into l and names them '<sublistnameinl.variablenameinsublist>' 
      # i.e. removes list|sublist hierarchy and retains this information in the variable name
      # prefixes the names of the resulting list with '<modobj>.' 
      comb_init_list <- function(., l, modobj ) {
        if(sum(is.null(l))==length(l)|is.null(l)) NULL
        else {
          # check for any sublists in v 
          slss <- vapply(l, is.list, logical(1))
          if(any(slss)) {
            l1 <- l[which(!slss)]
  
            for(ss in which(slss)) {
              names(l[[ss]]) <- paste0(names(l[ss]),'.',names(l[[ss]]))
              l1 <- c(l1,l[[ss]])
            }
          } else l1 <- l

          names(l1) <- paste(modobj, names(l1), sep='.' ) 
          l1 
        }
      }
     
      # setup list names for assignment 
      mos    <- c(.$model$name, unlist(.$model$child_list) )
      type   <- c('static', 'dynamic')
      vlists <- c('fnames', 'pars', 'env' )
      
      for( mo in mos ) {
        for( t in type ) {
          for( vl in vlists ) {
            # input variables
            varlist <- .[[paste0('init_',t)]][[mo]][[vl]]
 
            # assign variables to wrapper, prefix variable names with the name of the model object that they belong to 
            .[[t]][[vl]] <- if(is.null(.[[t]][[vl]])) comb_init_list(l=varlist, modobj=mo ) 
                            else                     c(.[[t]][[vl]], comb_init_list(l=varlist, modobj=mo ) )
          }
        }
      
        # as above for pars code snippets (pars_eval input) and assigment of parameters to a process (pars_proc input)
        if(.$wpars$UQ) {
          if(is.null(.$init_dynamic[[mo]]$pars)&!is.null(.$init_dynamic[[mo]]$pars_eval)) .$wpars$eval_strings <- T
          t <- 'dynamic'

          if(.$wpars$eval_strings) {
            vl   <- 'pars_eval'
            varlist <- .[[paste0('init_',t)]][[mo]][[vl]]
            .[[t]][[vl]] <- if(is.null(.[[t]][[vl]])) comb_init_list(l=varlist, modobj=mo ) 
                            else                     c(.[[t]][[vl]], comb_init_list(l=varlist, modobj=mo ) )
          }
  
          if(.$wpars$UQtype=='ye') {
            vl   <- 'pars_proc'
            varlist <- .[[paste0('init_',t)]][[mo]][[vl]]
            .[[t]][[vl]] <- if(is.null(.[[t]][[vl]])) comb_init_list(l=varlist, modobj=mo ) 
                            else                     c(.[[t]][[vl]], comb_init_list(l=varlist, modobj=mo ) )
            for( vn in names(vars) ) if( !any(vn==names(.[['init_dynamic']][[mo]][['pars_eval']])) )
              stop(paste('\n Input variable:', vn, 'in pars_proc, not found in: pars_eval list.',
                         '\n The proc_pars input list must contain exactly the same parameter names as pars_eval input list.',
                         '\n The proc_pars is required to assign a parameter to a process as part of a process sensitivity analysis.'))
          }
        }
      }
    }
    
    
    ###########################################################################
    # variables that are to be varied in the SA/UQ
    
    # initialisation lists
    init_static  <- NULL
    init_dynamic <- NULL
    
    # static variables
    # all expected to be of class 'list'
    # each list in the below list comprise a single string or numeric value for each variable, labelled by the variable name
    static <- list( 
      fnames = NULL,
      pars   = NULL,
      env    = NULL
    )
    
    # dynamic variables
    # all expected to be of class 'list'
    # each list in the 'vars' list comprise vectors of the values for each variable, 
    # each element of the list is labelled by the variable name prefixed by the name of the model object that the variable belongs to
    # each of these lists is expanded, often factorially by expand.grid, and placed into the below list of dataframes
    dynamic <- list( 
      fnames        = NULL,
      fnamesB       = NULL,
      pars          = NULL,
      # list with same elements and names as pars but giving the fnames list name i.e. the process name to which each parameter belongs
      pars_proc     = NULL,
      # list with same elements and names as pars but each element is a code snippet as a string that once evaluated gives a vector or parameter values
      # allows different types of distributions to be specified for each parameter
      # this must be used for the Ye SA method
      pars_eval     = NULL,
      pars_lklihood = NULL,
      pars_array    = NULL,
      parsB         = NULL,
      env           = NULL
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
      met     = NULL,         # a dataframe of sequential meteorological driving data, for running the analysis at a particular site for example 
      # row length of above matrices
      lf      = NULL,
      lfA     = NULL,
      lfB     = NULL,
      lp      = NULL,
      lpB     = NULL,
      le      = NULL,
      lm      = NULL,
      # output matrices / arrays
      mout         = NULL,    # example model output vector, for setting up vapply functions  
      out          = NULL,    # output matrix
      out_saltelli = NULL,    # saltelli output list
      # observation matrices /dataframes
      obs     = NULL,         # a dataframe of observations against which to valiadate/ calculate likelihood of model
      obsse   = NULL          # a dataframe of observation errors for the obs data, must exactly match the above dataframe
      
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
      sobol_init   = T,       # initialise sobol sequence or not when calling rsobol. This should not be modified by the user. 
      mcmc_maxiter = 100,     # MCMC maximum number of iterations / steps in the chain 
      mcmc_chains  = 10,      # MCMC number of chains 
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

          # run types if - SA/UQ or not
          # if SA/UQ this output function is not used 
          # SA/UQ functionality is here to allow the user to call these to inspect SA/UQ output against inputs in a common dataframe 
          if(.$wpars$UQ) {
            # if Ye (i.e. process) SA 
            if(.$wpars$UQtype=='ye') {
              
              # combine input into a single dataframe in order of output dataframe,
              #  - i.e. repeats lines in input dataframes/matrices to align with output
              #  - the number of rows in the resultant dataframe is lf*lfB*n^2*le
 
              vardf    <- list(
                # fnames of process A - length lf * lfB * le * n^2
                fnames  = apply(.$dataf$fnames, 2,function(v) rep(v,each=.$dataf$lfB*.$dataf$le*.$wpars$n^2) ),
                # pars of process A
                pars    = apply(.$dataf$pars,   2,function(v) rep(rep(v,each=.$dataf$lfB*.$dataf$le*.$wpars$n),.$dataf$lfA) ),
                # fnames of process B
                fnamesB = apply(.$dataf$fnamesB,2,function(v) rep(rep(v,each=.$dataf$le*.$wpars$n),.$dataf$lfA*.$wpars$n) ),
                # parsB
                parsB   = apply(.$dataf$parsB,  2,function(v) rep(v,each=.$dataf$le) ),
                # environment
                env     = apply(.$dataf$env,    2,function(v) rep(v,.$dataf$lfA*.$dataf$lfB*.$wpars$n^2) )
              )
              
            # if Saltelli SA
            } else if(.$wpars$UQtype=='saltelli') { 

              vardf <- cbind( 
                # fnames
                apply(.$dataf$fnames, 2,function(v) rep(v,each=.$dataf$le*.$dataf$lp) ),
                # pars - I think this has the correct length specifications, need to check witha saltelli unit testing function
                apply(.$dataf$pars,   2,function(v) rep(rep(v,each=.$dataf$le),.$dataf$lf) ),
                # environment
                apply(.$dataf$env,    2,function(v) rep(v,.$dataf$lf*.$dataf$lp) )
              )
              
              vardf <- as.data.frame(vardf)
            }   

          # if factorial combination run
          } else { 
            vpars    <- if(is.null(.$dataf$pars))    NULL else .$dynamic$pars
            venv     <- if(is.null(.$dataf$env))     NULL else .$dynamic$env
            vfnames  <- if(is.null(.$dataf$fnames))  NULL else .$dynamic$fnames

            vardf <- expand.grid(c(venv,vpars,vfnames),stringsAsFactors=F)        
          }
          
          # if no met data
          if(is.null(.$dataf$met)) {
            # if(.$wpars$UQtype=='ye') {
            if(.$wpars$UQ) {
              print(paste('vardf:',length(vardf[[1]][,1]),', out:',length(.$dataf$out[,1])),quote=F)
              # return a list
              return(c(vardf, list(out=.$dataf$out) ))
              rm(vardf)            
            } else {
              print(paste('vardf:',length(vardf[,1]),', out:',length(.$dataf$out[,1])),quote=F)
              # return a dataframe
              return(cbind(vardf, .$dataf$out ))
              rm(vardf)            
            }

          # if met data
          # - so far will only work for factorial simulations
          } else {
            print(vardf)  
            print(.$dataf)  
            odf <- cbind(do.call(rbind , lapply(1:length(vardf[,1]) , .$combine, df=vardf ) ), .$dataf$out )            
            print(head(odf))
            if(dim(vardf)[2]==1) names(odf)[which(names(odf)=='df.i...')] <- names(vardf)
            rm(vardf)
            return(odf)
            rm(odf)            
          } 

        # if no vars  
        } else {
          # if met data
          if(!is.null(.$dataf$met)) cbind(.$dataf$met , .$dataf$out ) else .$dataf$out  
        }
      )
    }

    output_saltelli_AB <- function(.) {
      # creates output for a saltelli Sobol sensitivity analysis 
      # A and B matrices are stacked in a single matrix, which for each model and environment combination are then stored in an array 

      # AB output is an array
      # - dim 1 (rows)   model combination
      # - dim 2 (cols)   environment combination
      # - dim 3 (slices) sample
      # - dim 4          output variable (character variables are coerced to NAs)

      # create AB output matrix array
      AB  <- array(.$dataf$out, c(.$dataf$le, 2*.$wpars$n, .$dataf$lf, length(.$dataf$mout) ))
      dimnames(AB) <- list( NULL, NULL, apply(.$dataf$fnames, 1, toString), names(.$dataf$mout)  )
      
      # output a list composed of the AB matrix output array, the fnames that define each model combination, the parameter names
      aperm(AB, c(3,1,2,4) )
    }
  
    output_saltelli_ABi <- function(.) {
      # creates output for a saltelli Sobol sensitivity analysis 

      # ABi output is an array
      # - dim 1 (rows)         model combination
      # - dim 2 (cols)         environment combination
      # - dim 3 (slices)       sample
      # - dim 4 (cube rows)    output variable
      # - dim 5 (cube columns) parameter that has used value from matrix B while all other par values are from matrix A
      
      # .$dataf$out_saltelli needs permuting to acheive above array dim order
      # - dim 1 (rows)      output variable
      # - dim 2 (columns)   sample
      # - dim 3 (slices)    parameter that has used value from matrix B while all other par values are from matrix A 
      # - dim 4 (cube rows) environment combination
      # - dim 5 (cube cols) model combination
      
      # aperm(.$dataf$out_saltelli, c(5,4,1:3) )
      aperm(.$dataf$out_saltelli, c(5,4,2,1,3) )
    }
    
    
    # Print functions
    ###########################################################################
    
    print_data <- function(.,otype='data') {
      
      ens_n <- 
        if(.$wpars$UQ) {
          if(.$wpars$UQtype=='ye') .$dataf$lf * prod(unlist(lapply(.$dynamic$fnames,length))) * .$wpars$n^2 * .$dataf$le
          else                     .$dataf$lf * .$wpars$n * (2+dim(.$dataf$pars)[2]) * .$dataf$le
        } else                     .$dataf$lf * .$dataf$lp *.$dataf$le
     
      if(otype=='data') {         
        
        print('',quote=F)
        print('',quote=F)
        print('',quote=F)
        print('',quote=F)
        print("MAAT :: summary of data",quote=F)
        print('',quote=F)
        print('',quote=F)
  
        print("fnames ::",quote=F)
        print(summary(.$dataf$fnames),quote=F)
        print('',quote=F)
        print("pars ::",quote=F)
        if(!.$wpars$UQtype=='ye') print(summary(.$dataf$pars),quote=F)
        else {
          print(.$dynamic$pars_proc,quote=F)
          print(paste('sample n:',.$wpars$n),quote=F)
        }                      
        print('',quote=F)
        print("env ::",quote=F)
        print(summary(.$dataf$env),quote=F)
        print('',quote=F)
        print("met data ::",quote=F)
        print(summary(.$dataf$met),quote=F)
        print('',quote=F)
      
      } else if(otype=='run') {
        run_type <- if(!.$wpars$UQ) 'Factorial' else paste('SA/UQ, ',.$wpars$UQtype)
        
        print('',quote=F)
        print('',quote=F)
        print('',quote=F)
        print(paste("MAAT :: run model",Sys.time()),quote=F)
        print('',quote=F)
        print(paste(run_type,' ensemble'),quote=F)
        print(paste('ensemble number ::',ens_n),quote=F)
        if(!is.null(.$dataf$met)) {
          print(paste('timesteps in met data ::',.$dataf$lm),quote=F)                
          print(paste('total number of model calls ::',ens_n*.$dataf$lm),quote=F)                
        }
        print('',quote=F)
        
        if(.$wpars$multic) print(paste('parallel processing over ::',.$wpars$procs,'cores.'),quote=F)
        else               print(paste('serial processing.'),quote=F)
        print('',quote=F)
      } 
    }
    
    print_output <- function(.) {
      print("output ::",quote=F)
      print(paste('length ::', length(.$dataf$out[,1])), quote=F)
      print(head(.$dataf$out), quote=F)
      print('', quote=F)      
      print('', quote=F)      
      print(Sys.time(), quote=F)      
      print('', quote=F)      
    }
    
    print_saltelli <- function(.) {
      print('Saltelli matrix AB completed', quote=F)
      print('', quote=F)
      print('', quote=F)
      print('run Saltelli array ABi', quote=F)
      print('', quote=F)
    }
    
    
    
    ######################################################################################
    # Unit testing functions

    # simple test, run a single model instance with or without met data
    .test_simple <- function(., metd=F, mc=F, pr=4, oconf=F ) {
      
      # source directory
      setwd('system_models/leaf')
      source('leaf_object.R')
      setwd('../..')
      
      # clone the model object
      .$model <- as.proto(leaf_object$as.list(),parent=.)
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- oconf      
      print(.$model$output())
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc  # multicore the ensemble
      .$wpars$procs        <- pr  # number of cores to use if above is true
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening 
      
      # Define meteorological and environment dataset
      .$model$env$par     <- 1000
      .$model$env$ca_conc <- 400  
      .$model$env$vpd     <- 1  
      .$model$env$temp    <- 20  
      metdata <- as.matrix(expand.grid(list(leaf.par = seq(800,1000,100),leaf.ca_conc = 400)))      
      if(metd) .$dataf$met <- metdata
      
      # Define the static parameters and model functions  
      .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')
      .$dynamic$fnames <- list(
        leaf.etrans = c('f_j_farquhar1980')
      )
      
      .$dynamic$env <- list(
        leaf.vpd  = 1
      )
      
      .$dynamic$pars <- list(
        leaf.avn_25 = 10
      )
      
      # Run wrapper & model
      st <- system.time(
        .$run()
      )
      print('',quote=F)
      print('Run time:',quote=F)
      print(st)
      print('',quote=F)
      
      # process & record output
      .$output()
    }    
    
 
    # general factorial test, with or without metdata
    .test <- function(.,metd=T,mc=T,pr=4,oconf=F) {
      
      # source directory
      setwd('system_models/leaf')
      source('leaf_object.R')
      setwd('../..')
      
      library(lattice)
      
      # clone the model object
      .$model <- as.proto(leaf_object$as.list(),parent=.)
      .$model$cpars$verbose  <- F      
      .$model$cpars$cverbose <- oconf      
      print(.$model$output())
      
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc  # multicore the ensemble
      .$wpars$procs        <- pr  # number of cores to use if above is true
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening 
      
      ### Define meteorological and environment dataset
      ###############################
      # can load a met dataset here
      # below a trivial met dataset is created to be used as an example
      metdata <- as.matrix(expand.grid(list(leaf.par = seq(0,1000,100),leaf.ca_conc = 400)))      
      metdata <- as.matrix(expand.grid(list(leaf.par = seq(800,1000,100),leaf.ca_conc = 400)))      
      if(metd) .$dataf$met <- metdata
      else     { .$model$env$par <- 1000; .$model$env$ca_conc <- 400 } 

      ### Define the static parameters and model functions  
      ###############################
      .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')
            
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length,
      # if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
      
      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$dynamic$fnames <- list(
        leaf.etrans = c('f_j_farquhar1980','f_j_collatz1991'),
        leaf.rs     = c('f_r_zero','f_rs_medlyn2011')
      )
        
      .$dynamic$env <- list(
        leaf.vpd  = c(1,2),
        leaf.temp = c(5,20)
      )
        
      .$dynamic$pars <- list(
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
      
      # # process & record output
      # df <- .$output()
      .$output()
      # p1 <- xyplot(A~leaf.ca_conc|leaf.etrans*leaf.rs,df,groups=leaf.temp,type='l',auto.key=T,
      #              panel=function(...) { panel.abline(h=seq(0,20,2.5)) ; panel.xyplot(...) })
      # list(df,p1)
    }    


    # general factorial test, with or without metdata
    .test_init <- function(., 
                          metd=list(leaf.par = seq(800,1000,100), leaf.ca_conc = 400 ),
                          sfnames=list(fnames=list(vcmax='f_vcmax_lin')),
                          spars=list(pars=list(Ha=list(vcmax=7e4, jmax=4e4 ))),
                          senv=list(env=list(par=1000, ca_conc=400 )),
                          dfnames=NULL,
                          dpars=NULL,
                          denv=NULL
                          ) {
 
      #
      .$model <- NULL
      .$model$name <- 'leaf'      

      # Define the static model functions, parameters, and environment  
      .$init_static$leaf <- c(sfnames, spars, senv ) 
      print('init_static')     
      print(.$init_static)     
       
      # Define the dynamic model functions, parameters, and environment  
      .$init_dynamic$leaf <- c(dfnames, dpars, denv ) 
      print('init_dynamic')      
      print(.$init_dynamic)      

      # Run init function
      .$init()
      print('static')     
      print(.$static)     
      print('dynamic')     
      print(.$dynamic)      
      
    }    


    # general factorial test, with or without metdata
    .test_con <- function(., mc=T, pr=4, oconf=F,
                          metd=list(leaf.par = seq(800,1000,100),leaf.ca_conc = 400),
                          sfnames=list(leaf.vcmax='f_vcmax_lin'),
                          spars=NULL,
                          senv=list(par=1000, ca_conc=400 ),
                          dfnames=list(
                            leaf.etrans = c('f_j_farquhar1980','f_j_collatz1991'),
                            leaf.rs     = c('f_r_zero','f_rs_medlyn2011')
                          ),
                          dpars=list(
                            leaf.avn_25 = 9:11,
                            leaf.bvn_25 = 4:6
                          ),
                          denv=list(
                            leaf.vpd  = c(1,2),
                            leaf.temp = c(5,20)
                          )
                          ) {
      
      # source directory
      setwd('system_models/leaf')
      source('leaf_object.R')
      setwd('../..')
      
      # clone the model object
      .$model <- as.proto(leaf_object$as.list(),parent=.)
      .$model$cpars$verbose  <- F      
      .$model$cpars$cverbose <- oconf      
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc  # multicore the ensemble
      .$wpars$procs        <- pr  # number of cores to use if above is true
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening 
      
      # Define meteorological dataset
      if(!is.null(metd)) .$dataf$met  <- as.matrix(expand.grid(metd))

      # Define the static model functions, parameters, and environment  
      if(!is.null(sfnames)) .$static$fnames <- sfnames 
      if(!is.null(spars))   .$static$pars   <- spars 
      if(!is.null(senv))    .$static$env    <- senv 
            
      # Define the dynamic model functions, parameters, and environment  
      if(!is.null(dfnames)) .$dynamic$fnames <- dfnames 
      if(!is.null(dpars))   .$dynamic$pars   <- dpars 
      if(!is.null(denv))    .$dynamic$env    <- denv 
      
      # Run wrapper & model
      st <- system.time(
        .$run()
      )
      print('',quote=F)
      print('Run time:',quote=F)
      print(st)
      print('',quote=F)
      
      # process & record output
      .$output()
    }    


    # general factorial test with canopy object, with or without metdata
    .test_can <- function(., metd=T, mc=T, pr=4, verbose=F ) {
     
      mod_obj <- 'canopy'
       
      # source directory
      #setwd('system_models/canopy')
      #source('canopy_object.R')
      setwd(paste('system_models',mod_obj,sep='/'))
      source(paste(mod_obj,'object.R',sep='_'))
      # clone & build the model object
      init_default <- .$build(model=paste(mod_obj,'object',sep='_') )
      setwd('../..') 
  
      # load MAAT object(s) from source
      #setwd(paste('system_models',mod_obj,sep='/'))
      #source(paste(mod_obj,'object.R',sep='_'))
      # clone & build the model object
      #init_default <- .$build(model=paste(mod_obj,'object',sep='_'), mod_mimic=mod_mimic )
      #setwd('../..')      
  
      # clone the model object
      #.$model      <- as.proto(canopy_object$as.list(),parent=.)
      #.$model$leaf <- as.proto(leaf_object$as.list(),parent=.$model)      

      # define parameters for the model
      .$model$pars$verbose       <- verbose      
      .$model$leaf$pars$cverbose <- verbose      
      .$model$state$mass_a       <- 175
      .$model$state$C_to_N       <- 40
      
      # define parameters for the wrapper
      .$wpars$multic       <- mc  # multicore the ensemble
      .$wpars$procs        <- pr  # number of cores to use if above is true
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if faslse a fully factorial ensemble 
      .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening 
      
      # define meteorological and environment dataset
      metdata <- expand.grid(list(canopy.par_dir = 500,canopy.ca_conc = seq(10,1200,50)))      
      if(metd) .$dataf$met <- metdata
      
      # Define the parameters and model functions that are to be varied 
      .$dynamic$fnames <- list(
        canopy.can_scale_light = c('f_canlight_beerslaw_wrong','f_canlight_beerslaw'),
        leaf.etrans            = c('f_j_farquhar1980','f_j_collatz1991'),
        leaf.rs                = c('f_rs_medlyn2011','f_r_zero')
      )
      
      .$dynamic$env <- list(
        leaf.vpd  = c(1,2),
        leaf.temp = c(5,20)
      )
      
      .$dynamic$pars <- list(
        canopy.lai = seq(2,6,2),
        canopy.G   = seq(0.4,0.6,0.1)
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
      library(lattice)
      p1 <- xyplot(A~canopy.ca_conc|canopy.can_scale_light*leaf.rs,df,groups=canopy.lai,type='l',abline=5,auto.key=T)
      list(df,p1)
    }
      

    # simple test of init_function & model mimic set up 
    .test_mimic <- function(., mod_mimic='clm45_non_Tacclimation', mod_obj='leaf', metd=F, mc=F, pr=4, oconf=F ) {
      
      # load MAAT object(s) from source
      setwd(paste('system_models',mod_obj,sep='/'))
      source(paste(mod_obj,'object.R',sep='_'))
      # clone & build the model object
      init_default <- .$build(model=paste(mod_obj,'object',sep='_'), mod_mimic=mod_mimic )
      setwd('../..')      
      
      # define control parameters 
      .$model$cpars$verbose  <- F      
      .$model$cpars$cverbose <- oconf      
      .$wpars$multic       <- mc  # multicore the ensemble
      .$wpars$procs        <- pr  # number of cores to use if above is true
      .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble 
      .$wpars$unit_testing <- F   # tell the wrapper to run init function 
      
      ### Define the static parameters and model functions  
      ###############################
      init_default$leaf$env$ca_conc <- 400
      init_default$leaf$env$par     <- 2000
      init_default$leaf$env$vpd     <- 50 
      init_default$leaf$env$temp    <- 25
      .$init_static <- init_default
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length,
      # if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
      #maat$init_dynamic <- init_dynamic
      .$init_dynamic <- NULL 
      
      ### Define meteorological and environment dataset
      ###############################
      # can load a met dataset here
      # below a trivial met dataset is created to be used as an example

      metdata <- as.matrix(expand.grid(list(leaf.par = seq(0,1000,100), leaf.ca_conc = 400))  )      
      metdata <- as.matrix(expand.grid(list(leaf.par = 2000, leaf.ca_conc = seq(50,1500,50))) )      
      if(metd) .$dataf$met <- metdata

      # Run model
      st <- system.time(
        .$run()
      )
      print('',quote=F)
      print('Run time:',quote=F)
      print(st)
      print('',quote=F)
      
      # # process & record output
      .$output()
      #library(lattice)
      # p1 <- xyplot(A~leaf.ca_conc|leaf.etrans*leaf.rs,df,groups=leaf.temp,type='l',auto.key=T,
      #              panel=function(...) { panel.abline(h=seq(0,20,2.5)) ; panel.xyplot(...) })
      # list(df,p1)
    }    


    # test function for Ye method Sobol process sensitivity analysis
    .test_ye <- function(.,metd=F,mc=T,pr=4,oconf=F,n=3) {
      
      # source directory
      setwd('system_models/leaf')
      source('leaf_object.R')
      # clone & build the model object
      init_default <- .$build(model='leaf_object')
      setwd('../..')
      
      library(lattice)
      
      # define control parameters
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- oconf      
      .$wpars$multic       <- mc   # multicore the ensemble
      .$wpars$procs        <- pr   # number of cores to use if above is true
      .$wpars$UQ           <- T    # run a UQ/SA style ensemble 
      .$wpars$UQtype       <- 'ye' # Ye style SA ensemble 
      .$wpars$unit_testing <- T    # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions) 
      
      ### Define static variables 
      ###############################
      .$model$env$par <- 1000
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # "pars" lists must contain parameter vectors that are of equal length,

      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$wpars$n            <- n    # number of parameter samples in each loop 
      .$wpars$coef_var     <- 0.1  # coefficient of variation for prior parameter distribution 
      .$wpars$eval_strings <- T    # use evaluation strings to set parameter values

      .$static$fnames <- list(vcmax='f_vcmax_lin')
      .$dynamic$fnames <- list(
        leaf.Alim   = c('f_lim_farquhar1980','f_lim_collatz1991'),
        leaf.etrans = c('f_j_farquharwong1984','f_j_collatz1991','f_j_harley1992')
      )

      .$dynamic$pars <- list(
        leaf.avn_25   = NA,
        leaf.bvn_25   = NA,
        leaf.theta_j  = NA,
        leaf.e_ajv_25 = NA
      )

      .$dynamic$pars_eval <- list(
        leaf.avn_25   = ' 10 * rnorm(n,1,.$wpars$coef_var)',
        leaf.bvn_25   = '  5 * rnorm(n,1,.$wpars$coef_var)',
        leaf.theta_j  = '0.9 * rnorm(n,1,.$wpars$coef_var)',
        leaf.e_ajv_25 = '0.9 * rnorm(n,1,.$wpars$coef_var)'
      )
      
      .$dynamic$pars_proc <- list(
        leaf.avn_25   = 'leaf.Alim',
        leaf.bvn_25   = 'leaf.Alim',
        leaf.theta_j  = 'leaf.etrans',
        leaf.e_ajv_25 = 'leaf.etrans'
      ) 

      .$dynamic$env <- list(
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
      
    }    
    
    # test function for Saltelli method Sobol parametric sensitivity analysis
    .test_saltelli <- function(., metd=F, mc=T, pr=4, oconf=F, n=3, eval_strings=T ) {
      
      # source directory
      setwd('system_models/leaf')
      source('leaf_object.R')
      # clone & build the model object
      init_default <- .$build(model='leaf_object')
      setwd('../..')
      library(lattice)
      
      # define control parameters
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- oconf      
      .$wpars$multic       <- mc           # multicore the ensemble
      .$wpars$procs        <- pr           # number of cores to use if above is true
      .$wpars$UQ           <- T            # run a UQ/SA style ensemble 
      .$wpars$UQtype       <- 'saltelli'   # Saltelli style SA ensemble 
      .$wpars$unit_testing <- T            # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions) 
      
      ### Define static variables 
      ###############################
      .$static$env    <- list(leaf.par=1000)
      .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')
      
      ### Define the parameters and model functions that are to be varied 
      ###############################
      # "pars" lists must contain parameter vectors that are of equal length,
      
      # add the SA/UQ variables to the maat wrapper object
      # - the wrapper object takes care of combining these lists into the full ensemble      
      .$wpars$n            <- n            # number of parameter samples in each loop 
      .$wpars$eval_strings <- eval_strings # parameters are passed as strings to be evaluated to allow for different sample numbers 
      .$wpars$coef_var     <- 0.1
      
      .$dynamic$fnames <- list(
        leaf.Alim   = c('f_lim_farquhar1980','f_lim_collatz1991'),
        leaf.etrans = c('f_j_farquharwong1984','f_j_collatz1991','f_j_harley1992')
      )

      if(eval_strings) {
        .$dynamic$pars_eval <- list(
          leaf.avn_25   = ' 10 * rnorm(n,1,.$wpars$coef_var)',
          leaf.bvn_25   = '  5 * rnorm(n,1,.$wpars$coef_var)',
          leaf.theta_j  = '0.9 * rnorm(n,1,.$wpars$coef_var)',
          leaf.e_ajv_25 = '0.9 * rnorm(n,1,.$wpars$coef_var)'
        )
      } else {
        n <- 2 * n
        .$dynamic$pars <- list(
          leaf.avn_25   =  10 * rnorm(n,1,.$wpars$coef_var),
          leaf.bvn_25   =   5 * rnorm(n,1,.$wpars$coef_var),
          leaf.theta_j  = 0.9 * rnorm(n,1,.$wpars$coef_var),
          leaf.e_ajv_25 = 0.9 * rnorm(n,1,.$wpars$coef_var)
        )
      }
      
      .$dynamic$env <- list(
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
      list(AB=.$output_saltelli_AB(), ABi=.$output_saltelli_ABi() )
    }    
    
    # test function for Saltelli method Sobol parametric sensitivity analysis
    .test_mcmc_mixture <- function(., mc=T, pr=4, mcmc_chains=4 ) {
      
      library(lattice)
      
      # define control parameters
      .$model$pars$verbose  <- F      
      .$model$pars$cverbose <- F      
      .$wpars$multic       <- mc           # multicore the ensemble
      .$wpars$procs        <- pr           # number of cores to use if above is true
      .$wpars$UQ           <- T            # run a UQ/SA style ensemble 
      .$wpars$UQtype       <- 'mcmc'       # MCMC ensemble 
      .$wpars$unit_testing <- T            # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions) 
     
      # met dataf must be non-NULL

      # model output must be scalar
 
      # define the mixture model
      .$model <- function(.) {
        # write mixture model here
      } 

      # define priors - the mixture model needs to be a proto object with a configure, run, output functions etc ...
      

      # Run MCMC 
      st <- system.time(
        .$run()
      )
      print('',quote=F)
      print('Run time:',quote=F)
      print(st)
      print('',quote=F)

      # process & record output
      # output   
 
    }  
    
 
###########################################################################
# end maat wrapper 
}) 



### END ###
