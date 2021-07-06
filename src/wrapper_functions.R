################################
#
# Wrapper functions for wrapper object
#
# AWalker December 2015
#
################################

library(parallel)
library(lattice)

source('functions/general_functions.R')
source('functions/calc_functions.R')
source('wrapper_functions_mcmc.R')


#####################################

# expand the fnames, env, and pars input lists into matrices
# - each column is passed to the model sequentially by each run function in the run function cascade
generate_ensemble <- function(.) {

  .$dataf$fnames  <- if(!is.null(.$dynamic$fnames)) t(as.matrix(expand.grid(.$dynamic$fnames,stringsAsFactors=F))) else NULL
  .$dataf$env     <- if(!is.null(.$dynamic$env))    t(as.matrix(expand.grid(.$dynamic$env,stringsAsFactors=F   ))) else NULL
  .$generate_ensemble_pars()

  # remove potentially large data from pars list
  if(.$wpars$runtype!='factorial') .$dynamic$pars <- lapply(.$dynamic$pars, function(e) numeric(1) )

  # check names of ensemble matrices are character vectors
  if(!is.null(.$dataf$met)) .$model$configure_check(vlist='env', df=.$dataf$met[,1] )

  # calculate input matrix lengths
  # - used to set the number of iterations in the run functions
  # - if no matrix return 1
  if(.$wpars$UQ & .$wpars$runtype=='SAprocess_ye') {
    .$dataf$lf <- length(.$dynamic$fnames)

  # any type of run other than Ye process sensitivity analysis
  } else {
    .$dataf$lf <- if(is.null(.$dataf$fnames)) 1 else dim(.$dataf$fnames)[2]
    .$dataf$lp <- if(is.null(.$dataf$pars))   1 else dim(.$dataf$pars)[2]
  }

  # enviroment matrix and met matrix
  .$dataf$le <- if(is.null(.$dataf$env)) 1 else dim(.$dataf$env)[2]
  .$dataf$lm <- if(is.null(.$dataf$met)) 1 else dim(.$dataf$met)[2]
}


# parameter matrix for factorial run
generate_ensemble_pars_factorial <- function(.) {
  if(!.$wpars$parsinit_read) {
    .$dataf$pars <- if(!is.null(.$dynamic$pars)) t(as.matrix(expand.grid(.$dynamic$pars,stringsAsFactors=F))) else NULL
  }
}


# parameter matrix for Saltelli parameter SA
generate_ensemble_pars_SApar_saltelli <- function(.) {

  # sample parameters 
  if(!.$wpars$parsinit_read) {
    if(is.null(.$dynamic$pars)) {

      #from character string code snippets to generate matrices A and B
      if(!is.null(.$dynamic$pars_eval)) {

        .$wpars$n      <- .$wpars$n * .$wpars$nmult
        n              <- 2 * .$wpars$n
        .$dynamic$pars <- lapply(.$dynamic$pars_eval,function(cs) eval(parse(text=cs)))
        .$dataf$pars   <- t(do.call(cbind, .$dynamic$pars ))

      } else  stop('wrapper: pars (or pars_eval) list in dynamic list is empty')
    }

    ## remove potentially large pars list
    #.$dynamic$pars <- lapply(.$dynamic$pars, function(e) numeric(1) )
  }
}


# parameter matrix for Dai, Ye, process SA
# - the paramter matrices for the Ye method are generated in run function 1
# - here a list structure in .$dynamic$pars is created from .$dynamic$pars_eval
generate_ensemble_pars_SAprocess_ye <- function(.) {

  # need a minimum of >1 processes
  if(dim(.$dataf$fnames)[1]<=1) stop('need more than one process for a process sensitivity analysis')

  # check input dynamic$pars* are same length
  test_in <- length(.$dynamic$pars_eval) - length(.$dynamic$pars_proc)
  if(test_in!=0) stop('wrapper: Parameter input vectors - pars_eval & pars_proc - are not the same length')

  # assign same list structure as dynamic$pars_eval to dynamic$pars
  .$dynamic$pars <- lapply(.$dynamic$pars_eval,function(e) numeric(1) )

  # check input dynamic$pars* elements have same names
  # - to be done
}


# parameter matrix for initial proposal of MCMC
generate_ensemble_pars_mcmc <- function(.) {

  # determine boundary handling limits for parameter space
  .$boundary_handling_set()

  # read character string code snippets & sample initial proposal from prior
  if(!.$wpars$parsinit_read) {
    if(is.null(.$dynamic$pars)) {
    if(!is.null(.$dynamic$pars_eval)) {

      if(.$wpars$mcmc$mcmc_type=='dream') {
        n <- .$wpars$mcmc$chains
      } else if(.$wpars$mcmc$mcmc_type=='dreamzs') {
        n <- .$wpars$mcmc$prior_n 
      }
      .$dynamic$pars <- lapply(.$dynamic$pars_eval,function(cs) eval(parse(text=cs)))
      .$dataf$pars   <- t(do.call(cbind, .$dynamic$pars ))

      # boundary handling in case of bounded normal (or other bounded) distributions as prior 
      pm   <- .$dataf$pars
      bmax <- .$mcmc$boundary_max 
      bmin <- .$mcmc$boundary_min
      .$dataf$pars[] <- 
        t(sapply(1:dim(pm)[1], function(p, pm2=pm[p,], min2=bmin[p], max2=bmax[p], pt=.$dynamic$pars_eval[[p]]  )
           ifelse(pm2==min2 | pm2==max2, eval(parse(text=pt)), 
                  pm2) 
           ))

    } else  stop('wrapper: pars (or pars_eval) list in dynamic list is empty')
  }}

  # for DREAM-ZS reconfigure pars array 
  # - initialise past state array
  if(.$wpars$mcmc$mcmc_type=='dreamzs') {
    if(!.$wpars$parsinit_read) {
      lps     <- NULL
      lps_new <- dim(.$dataf$pars)[2] + .$wpars$mcmc$chains*floor(.$wpars$mcmc$maxiter/.$wpars$mcmc$iterappend)
      .$dataf$past_states    <- matrix(0, dim(.$dataf$pars)[1], lps_new )
      .$dataf$past_states[,1:.$wpars$mcmc$prior_n] <- .$dataf$pars
      .$dataf$pars           <- .$dataf$pars[,(.$wpars$mcmc$prior_n-.$wpars$mcmc$chains+1):.$wpars$mcmc$prior_n]
      .$dataf$lps            <- .$wpars$mcmc$prior_n

    # a restart  
    } else {
      lps     <- dim(.$dataf$mcmc_input$past_states)[2]
      lps_new <- lps + .$wpars$mcmc$chains*floor(.$wpars$mcmc$maxiter_restart/.$wpars$mcmc$iterappend)
      .$dataf$past_states         <- matrix(0, dim(.$dataf$pars)[1], lps_new )
      .$dataf$lps                 <- lps 
      .$dataf$past_states[,1:lps] <- .$dataf$mcmc_input$past_states
      print('')
      print('past states array')
      print(.$dataf$past_states[,(lps-10):lps])
    }
  }

  # observation subsampling
  # - evenly spaced or random subsampling, currently manual switch
  if(!.$wpars$parsinit_read & .$wpars$mcmc$thin_obs<1.0) {
    if(.$wpars$mcmc$thin_obs_random) {
      oss                   <- which(runif(dim(.$dataf$met)[2]) < .$wpars$mcmc$thin_obs)
    } else {
      .$wpars$mcmc$thin_obs <- min(0.5, .$wpars$mcmc$thin_obs )
      thin                  <- floor(1/.$wpars$mcmc$thin_obs)
      oss                   <- seq(1, dim(.$dataf$met)[2], thin )
    }
    .$dataf$met           <- .$dataf$met[,oss]
    .$dataf$obs           <- .$dataf$obs[oss]
    if(!is.null(.$dataf$obsse)) .$dataf$obsse <- .$dataf$obsse[oss]
  }
}



################################
# initialise output matrix/array
# - these functions initialise the output matrices used to store output from the ensemble

init_output_matrix_factorial <- function(.) {
  .$dataf$out <- matrix(0, .$dataf$lm*.$dataf$le*.$dataf$lp*.$dataf$lf, length(.$dataf$mout) )
  colnames(.$dataf$out) <- names(.$dataf$mout)
}

init_output_matrix_SApar_saltelli <- init_output_matrix_factorial


init_output_matrix_SApar_saltelli_ABi <- function(.) {
  # - dim 1 (rows)      output variable
  # - dim 2 (columns)   sample
  # - dim 3 (slices)    parameter that has used value from matrix B while all other par values are from matrix A
  # - dim 4 (cube rows) environment combination
  # - dim 5 (cube cols) model combination

  .$dataf$out_saltelli <- array(0, dim=c(length(.$dataf$mout), .$wpars$n, dim(.$dataf$pars)[1], .$dataf$le, .$dataf$lf ))
  dimnames(.$dataf$out_saltelli) <- list(names(.$dataf$mout), NULL, rownames(.$dataf$pars), NULL, apply(.$dataf$fnames, 2, toString) )
}


init_output_matrix_SAprocess_ye <- function(.) {
  # Ye method does not generate a single ensemble for whole simulation but rather a separate ensemble for each process
  # - dim 1 (rows)        output variable
  # - dim 2 (columns)     environment combination
  # - dim 3 (slices)      process(es) B parameter sample
  # - dim 4 (cube rows)   process(es) B representation(s)
  # - dim 5 (cube cols)   process A parameter sample
  # - dim 6 (cube slices) process A representation

  .$dataf$out           <- array(0, c(length(.$dataf$mout), .$dataf$le, .$wpars$n, .$dataf$lfB, .$wpars$n, .$dataf$lfA  ) )
  dimnames(.$dataf$out) <- list(names(.$dataf$mout), NULL, NULL, apply(.$dataf$fnamesB, 2, toString), NULL, .$dataf$fnames[,] )
}


init_output_matrix_mcmc <- function(.) {

  # number of parameters being estimated
  .$mcmc$pars_n   <- dim(.$dataf$pars)[1]
  .$mcmc$sd_state <- numeric(.$mcmc$pars_n)

  # preallocate memory space for algorithmic variables
  .$mcmc$current_state <- matrix(0, nrow=.$mcmc$pars_n, ncol=.$wpars$mcmc$chains )
  .$mcmc$jump          <- matrix(0, nrow=.$mcmc$pars_n, ncol=.$wpars$mcmc$chains )
  .$mcmc$lambda        <- numeric(.$wpars$mcmc$chains)

  # initialise accepted proposal & likelihood arrays
  .$dataf$pars_array    <- array(1, dim=c(dim(.$dataf$pars), .$wpars$mcmc$maxiter ))
  # APW: parameter names could be associated with dimnames of the array sooner 
  row.names(.$dataf$pars_array) <- row.names(.$dataf$pars)
  .$dataf$pars_lklihood <- matrix(1, .$wpars$mcmc$chains, .$wpars$mcmc$maxiter )

  # initialise output matricies
  .$dataf$out           <- matrix(0, .$dataf$lp, .$dataf$lm )
  .$dataf$out_mcmc      <- array(0, dim=c(.$dataf$lm, .$dataf$lp, .$wpars$mcmc$maxiter ))

  # initialise matricies for storing outlier & convergence diagnostics
  check_iter_n <-
    if(!.$wpars$parsinit_read) { 
      ceiling((.$wpars$mcmc$maxiter-.$wpars$mcmc$preburnin_iter+1) / .$wpars$mcmc$check_iter )
    } else {
      cin <- dim(.$dataf$mcmc_input$conv_check)[2] 
      cin <- cin + floor(.$wpars$mcmc$maxiter/.$wpars$mcmc$check_iter) - floor(.$wpars$mcmc$start_iter/.$wpars$mcmc$check_iter)
      if(.$wpars$mcmc$maxiter%%.$wpars$mcmc$check_iter!=0) cin <- cin + 1
      cin
    }
  .$dataf$omega         <- matrix(NA, .$wpars$mcmc$chains, check_iter_n )
  .$dataf$conv_check    <- matrix(0, ncol=check_iter_n, nrow=(dim(.$dataf$pars)[1]+3) )

  # initialise crossover variables if not a restart
  if(!.$wpars$parsinit_read) {
    .$mcmc$CR_counter      <- numeric(.$wpars$mcmc$n_CR)
    .$mcmc$jump_delta_norm <- numeric(.$wpars$mcmc$n_CR)
    .$mcmc$p_CR            <- numeric(.$wpars$mcmc$n_CR)
    .$mcmc$p_CR[]          <- 1/.$wpars$mcmc$n_CR
    .$mcmc$CR              <- numeric(.$wpars$mcmc$chains)
  
  # assign values from restart
  } else {

    .$dataf$out_mcmc[,,1:.$wpars$mcmc$start_iter-1]     <- .$dataf$mcmc_input$out_mcmc
    .$dataf$pars_array[,,1:.$wpars$mcmc$start_iter-1]   <- .$dataf$mcmc_input$pars_array
    .$dataf$pars_lklihood[,1:.$wpars$mcmc$start_iter-1] <- .$dataf$mcmc_input$pars_lklihood
    .$mcmc$check_ss                        <- dim(.$dataf$mcmc_input$conv_check)[2]
    .$dataf$conv_check[,1:.$mcmc$check_ss] <- .$dataf$mcmc_input$conv_check
    .$dataf$omega[,1:.$mcmc$check_ss]      <- .$dataf$mcmc_input$omega

    # remove potentialy large input data
    .$dataf$mcmc_input <- NULL
  }
}



###########################################################################
# run functions for the run cascasde
###########################################################################

################################
# for factorial runs or matrix A and B of Saltelli method for parametric sensitivity analysis (SA)
# run0
run0_factorial <- function(.) {

  # if a Saltelli SA and a met dataset has been specified, stop
  if(!is.null(.$dataf$met)&.$wpars$runtype=='SApar_saltelli') stop('No current method to run Saltelli SA with a met dataset')

  # generate output matrix
  .$init_output_matrix()

  # call run function
  .$dataf$out[] <-
    do.call( 'rbind', {
      if(.$wpars$multic) mclapply( 1:.$dataf$lf, .$run1, mc.cores=min(.$dataf$lf,.$wpars$procs), mc.preschedule=F )
      else                 lapply( 1:.$dataf$lf, .$run1 )
    })

  # call write function
  .$write_output()

  # if Saltelli run ABi ensemble
  if(.$wpars$runtype=='SApar_saltelli') .$run4()
  else                                  .$print_output()
}


# run1
run1_factorial <- function(.,i) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$fnames matrix to the model
  # assumes that each column of the fnames matrix are independent and non-sequential
  # call run2

  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[,i], F )
  if(.$wpars$cverbose)         .$printc('fnames', .$dataf$fnames[,i] )

  # call next run function
  do.call( 'rbind', {
      if(.$wpars$multic) mclapply(1:.$dataf$lp, .$run2, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lf)), mc.preschedule=T  )
      else                 lapply(1:.$dataf$lp, .$run2 )
  })
}


# run2
run2_factorial <- function(.,j) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$pars matrix to the model
  # assumes that each column of the pars matrix are independent and non-sequential
  # call run3

  # configure parameters in the model
  if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[,j], F )
  if(.$wpars$cverbose)       .$printc('pars', .$dataf$pars[,j] )

  # call next run function
  funv   <- if(is.null(.$dataf$met)) .$dataf$mout else array(0, dim=c(.$dataf$lm, length(.$dataf$mout) ) )
  out    <- vapply(1:.$dataf$le, .$run3, funv )

  # out has the potential to be a vector, matrix (needs transposed), or an array (needs stacking)
  # returns matrix
  #print('')
  #print(class(out))
  #print(out)
  if(class(out)[1]=='matrix') t(out) else if(class(out)[1]=='array') .$stack(out) else as.matrix(out)
}


# run3
run3_factorial <- function(.,k) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$env matrix to the model
  # assumes that each column of the env matrix are independent and non-sequential
  # call .$model$run or .$model$run_met if met data are provided

  # configure environment in the model
  if(!is.null(.$dataf$env)) .$model$configure(vlist='env', df=.$dataf$env[,k], F )
  if(.$wpars$cverbose)      .$printc('env', .$dataf$env[,k] )

  #print(.$init_static)
  #print(.$static)
  #print(.$model$pars)
  #print(.$model$state)
  print(.$model$fnames)
  #print(.$init_dynamic)
  #print(.$model$env)
  #print(.$model$fns$input())

  # call next run function
  if(is.null(.$dataf$met)) .$model$run() else .$model$run_met()
}


run4_factorial <- NULL
run5_factorial <- run4_factorial
run6_factorial <- run4_factorial
run7_factorial <- run4_factorial
run8_factorial <- run4_factorial



################################
# Run functions for Saltelli parameter SA method

# Copy above for Saltelli parameter SA matricies AB
run0_SApar_saltelli <- run0_factorial
run1_SApar_saltelli <- run1_factorial
run2_SApar_saltelli <- run2_factorial
run3_SApar_saltelli <- run3_factorial


# Additional run functions for ABi array of Saltelli parameter SA method
run4_SApar_saltelli <- function(.) {

  # print Saltelli
  .$print_saltelli()

  # load ABi specific functions
  .$init_output_matrix <- get('init_output_matrix_SApar_saltelli_ABi')
  .$write_output       <- get('write_output_SApar_saltelli_ABi')
  .$output             <- get('output_SApar_saltelli_ABi')

  # generate output matrix
  .$init_output_matrix()

  # run over ABi matrices
  .$dataf$out_saltelli[] <- vapply(1:.$dataf$lf, .$run5, .$dataf$out_saltelli[,,,,1] )

  # generate output
  .$write_output()

  # post-processing
  print(paste('Saltelli array ABi completed',Sys.time()),quote=F)
  print('',quote=F)
}


run5_SApar_saltelli <- function(.,i) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$fnames matrix to the model
  # assumes that each column of the fnames matrix are independent and non-sequential
  # call run6

  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames',df=.$dataf$fnames[,i],F)
  if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[,i])

  # call next run function
  vapply({
    if(.$wpars$multic) mclapply(1:.$dataf$le, .$run6, mc.cores=min(.$dataf$le,.$wpars$procs), mc.preschedule=F )
    else                 lapply(1:.$dataf$le, .$run6 )
  }, function(a) a, .$dataf$out_saltelli[,,,1,1] )
}


run6_SApar_saltelli <- function(.,k) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$env matrix to the model
  # assumes that each column of the fnames matrix are independent and non-sequential
  # call run7

  # configure environment in the model
  if(!is.null(.$dataf$env)) .$model$configure(vlist='env', df=.$dataf$env[,k], F )
  if(.$wpars$cverbose)      .$printc('env', .$dataf$env[,k] )

  # call parameter matrix run function
  if(is.null(.$dataf$met)){

    # call next run function, wrapped within vapply to convert (mc)lapply list output to an array
    # returns a numeric array - model output variable (rows), sample (columns), parameter (slices)
    vapply({
      #if(.$wpars$multic) mclapply(1:dim(.$dataf$pars)[2], .$run7, mc.cores=max(1,floor(.$wpars$procs/.$dataf$le)), mc.preschedule=T )
      #else                 lapply(1:dim(.$dataf$pars)[2], .$run7 )
      if(.$wpars$multic) mclapply(1:dim(.$dataf$pars)[1], .$run7, mc.cores=max(1,floor(.$wpars$procs/.$dataf$le)), mc.preschedule=T )
      else                 lapply(1:dim(.$dataf$pars)[1], .$run7 )
    },function(a) a, .$dataf$out_saltelli[,,1,1,1] )

  } else {
    # met data run not yet supported with Sobol, but should be caught before getting here
    stop('Saltelli cannot be run with met data')
  }
}


run7_SApar_saltelli <- function(.,p) {
  # This wrapper function is called from an lapply or mclappy function to be run once for each parameter (i.e. each column of the dataf$pars matrix)
  # call run8

  # returns a numeric matrix - model output variable (rows), sample (columns)
  vapply((.$wpars$n+1):.$dataf$lp, .$run8, .$dataf$mout, pk=p )
}


run8_SApar_saltelli <- function(.,j,pk) {
  # This wrapper function is called from an lapply or mclappy function
  # wrapper subscripts parameter matrix AB to give the column on matrix ABi
  # assumes that each column of the matrix are independent and non-sequential
  # call .$model$run

  # create index matrix to create column on matrix ABi for the .$dataf$par matrix (which is matrix A stacked on top of matrix B)
  #sub     <- rep(j-.$wpars$n, dim(.$dataf$pars)[2] )
  sub     <- rep(j-.$wpars$n, dim(.$dataf$pars)[1] )
  sub[pk] <- j
  #smat    <- cbind(sub, 1:dim(.$dataf$pars)[2] )
  smat    <- cbind(1:dim(.$dataf$pars)[1], sub )

  # create a matrix from vector and add names
  #psdf        <- t(.$dataf$pars[smat])
  #names(psdf) <- colnames(.$dataf$pars)
  # APW: perhaps an issue here
  psdf        <- .$dataf$pars[smat]
  names(psdf) <- rownames(.$dataf$pars)

  # configure parameters in the model
  if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=psdf, F )
  if(.$wpars$cverbose) .$printc('pars', psdf )

  # run model
  .$model$run()
}



################################
# for Process Sensitivity Analysis after Dai, Ye etal 2017

# The Ye method is a nested system of loops implemented by apply type functions
# Loop 1: switch and initialise process A process B loop
# Loop 2: process loop for process A               - use standard location for variable function values
# Loop 3: parameter loop for process A             - use standard location for variable parameter values
# Loop 4: process representation loop for proces B - use an additional non-standard location for variable function values
# Loop 5: parameter loop for proces B              - use an additional non-standard location for variable parameter values
# Loop 6: environment loop                         - use standard location for variable environment values

# run0
run0_SAprocess_ye <- function(.) {
  # Ye process SA is not multicored at this stage as multicoring here messes with the processes A and B in the data structure
  vapply(1:.$dataf$lf, .$run1, numeric(0) )
}


# run1
run1_SAprocess_ye <- function(.,f) {
  # this function is the overall wrapper function to run a generic process sensitivity analysis
  # The principle is to create a loop that runs the process SA nested loops once for each process to be analysed
  # This separates out the process in question - process A - from the other process(es) - process B (can be more than one process).
  # This function partitions the parameters to process A and and process B,
  # then creates the fnames and pars matrices for process A and process B,
  # calls the below run function,
  # outputs an .RDS for each process segregation

  # create the fnames matrices for process A and process B
  .$dataf$fnames  <- if(!is.na(.$dynamic$fnames[f]))       t(as.matrix(expand.grid(.$dynamic$fnames[f] ,stringsAsFactors=F))) else stop()
  .$dataf$fnamesB <- if(!any(is.na(.$dynamic$fnames[-f]))) t(as.matrix(expand.grid(.$dynamic$fnames[-f],stringsAsFactors=F))) else stop()

  # determine the number of the columns in the fnames process matrices
  #.$dataf$lfA     <- if(is.null(.$dataf$fnames )) 1 else length(.$dataf$fnames[,1])
  #.$dataf$lfB     <- if(is.null(.$dataf$fnamesB)) 1 else length(.$dataf$fnamesB[,1])
  .$dataf$lfA     <- if(is.null(.$dataf$fnames )) 1 else dim(.$dataf$fnames)[2]
  .$dataf$lfB     <- if(is.null(.$dataf$fnamesB)) 1 else dim(.$dataf$fnamesB)[2]

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
  .$dataf$pars    <- if(!is.na(.$dynamic$pars[1])) t(do.call(cbind,.$dynamic$pars[.$procA_subs] )) else stop()
  .$dataf$parsB   <- if(!is.na(.$dynamic$pars[2])) t(do.call(cbind,.$dynamic$pars[-.$procA_subs])) else stop()
  .$dynamic$pars  <- lapply(.$dynamic$pars_eval,function(e) numeric(1) )

  # determine the number of the columns in parameter matrices
  #.$dataf$lp      <- .$wpars$n # convert these to be the column number of the actual matrices
  #.$dataf$lpB     <- .$dataf$lfA * .$dataf$lfB * .$wpars$n^2
  .$dataf$lp      <- dim(.$dataf$pars)[2]
  .$dataf$lpB     <- dim(.$dataf$parsB)[2]

  # initialise output array
  .$init_output_matrix()

  # call run2 function
  print('',quote=F)
  print(paste('started process:', colnames(.$dataf$fnames), Sys.time()), quote=F )

  .$dataf$out[] <- vapply({
    if(.$wpars$multic) mclapply(1:.$dataf$lfA, .$run2, mc.cores=min(.$dataf$lfA,.$wpars$procs), mc.preschedule=F )
    else                 lapply(1:.$dataf$lfA, .$run2 )
  }, function(a) a, .$dataf$out[,,,,,1] )

  # process & record output
  if(.$wpars$unit_testing) { hd <- getwd(); setwd('~/tmp'); ofname <- 'Ye_test' }
  else                     setwd(odir)
  .$write_output(f=f)

  print(paste('completed process:', colnames(.$dataf$fnames), Sys.time() ), quote=F )

  # return nothing
  numeric(0)
}


# run2
run2_SAprocess_ye <- function(.,g) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$fnames matrix to the model
  # assumes that each column of the fnames matrix are independent and non-sequential
  # call run3

  print(paste('started representation:', .$dataf$fnames[,g], ', of process:', colnames(.$dataf$fnames)), quote=F )

  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[,g] , F )
  if(.$wpars$cverbose) .$printc('fnames', .$dataf$fnames[,g] )

  # calculate offset to correctly subset parsB matrix
  osg <- .$wpars$n * .$dataf$lfB * (g-1)

  # call process A parameter run function
  vapply({
    if(.$wpars$multic) mclapply(1:.$dataf$lp, .$run3, offset=osg, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lfA)), mc.preschedule=F  )
    else                 lapply(1:.$dataf$lp, .$run3, offset=osg )
  }, function(a) a, .$dataf$out[,,,,1,1] )
}


# run3
run3_SAprocess_ye <- function(., h, offset ) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$pars matrix to the model
  # assumes that each column of the pars matrix are independent and non-sequential
  # call run4

  # configure parameters in the model
  if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[,h], F )
  if(.$wpars$cverbose) .$printc('pars', .$dataf$pars[,h] )

  # calculate offset to correctly subset parsB matrix
  osh  <- offset + .$dataf$lfB * (h-1)

  # call process B process representation run function
  vapply(1:.$dataf$lfB, .$run4, .$dataf$out[,,,1,1,1], offset=osh )
}


# run4
run4_SAprocess_ye <- function(., i, offset ) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$fnamesB matrix to the model
  # assumes that each column of the fnamesB matrix are independent and non-sequential
  # call run5

  # configure function names in the model
  if(!is.null(.$dataf$fnamesB)) .$model$configure(vlist='fnames', df=.$dataf$fnamesB[,i], F )
  if(.$wpars$cverbose) .$printc('fnames', .$dataf$fnamesB[,i] )

  # calculate offset to correctly subset parsB matrix
  os  <- offset + i
  # oss is a vector of the column subscripts for the parsB matrix
  oss <- (.$wpars$n*(os-1) + 1):(.$wpars$n*(os))

  # call process B parameter run function
  vapply(oss, .$run5, .$dataf$out[,,1,1,1,1] )
}


# run5
run5_SAprocess_ye <- function(.,j) {
  # This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$parsB matrix to the model
  # assumes that each column of the parsB matrix are independent and non-sequential
  # call run6

  # configure parameters in the model
  if(!is.null(.$dataf$parsB)) .$model$configure(vlist='pars', df=.$dataf$parsB[,j], F )
  if(.$wpars$cverbose)        .$printc('pars', .$dataf$parsB[,j] )

  # call the environment run function
  vapply(1:.$dataf$le, .$run6, .$dataf$out[,1,1,1,1,1] )
}


run6_SAprocess_ye <- run3_factorial # i.e. over .$dataf$env
run7_SAprocess_ye <- run4_factorial # i.e. NULL
run8_SAprocess_ye <- run4_factorial # i.e. NULL



################################
# for MCMC runs

run0_mcmc <- function(.) {

  # MCMC run cannot run with more than one model
  if(.$dataf$lf>1) {
    print('MCMC cannot run with more than one process representation due to')
    print('current output file naming. Also to avoid the risk of specifying')
    print('parameters to estimate that do not belong to a specified process representation.')
    stop('More than one process representation requested, not compatible with MCMC run.')
  }

  # if more than one model output has been specified, stop
  if(length(.$dataf$mout)!=1) stop('No current method to run MCMC with multiple model outputs')

  # initialise output arrays & other MCMC data structures
  .$init_output_matrix()

  # call run function
  #if(.$wpars$multic) mclapply( 1:.$dataf$lf, .$run1, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lp)), mc.preschedule=F )
  #else                 lapply( 1:.$dataf$lf, .$run1 )
  vapply( 1:.$dataf$lf, .$run1, numeric(0) )

  # print summary of results
  .$print_output()
}


run1_mcmc <- function(.,i) {
  # assumes that each column of the fnames matrix are independent and non-sequential
  # call run2

  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[,i], F )
  if(.$wpars$cverbose)         .$printc('fnames', .$dataf$fnames[,i] )

  # evaluate model over initial proposals derived from prior
  # - if not a restart
  if(!.$wpars$parsinit_read) {
    .$dataf$out[]  <-
      do.call( 'rbind', {
          if(.$wpars$multic) mclapply(1:.$dataf$lp, .$run3, mc.cores=min(.$wpars$procs,.$dataf$lp), mc.preschedule=T  )
          else                 lapply(1:.$dataf$lp, .$run3 )
      })

    # add to pars array and calculate likelihood of initial proposal
    .$dataf$pars_array[,,1]   <- .$mcmc$prior_sample <- .$dataf$pars
    .$dataf$pars_lklihood[,1] <- .$proposal_lklihood()
    .$dataf$out_mcmc[,,1]     <- t(.$dataf$out[,]) 
  }

  # run MCMC
  vapply(.$wpars$mcmc$start_iter:.$wpars$mcmc$maxiter, .$run2, numeric(0) )

  # write output from MCMC
  .$write_output(i=i)

  numeric(0)
}


# This wrapper function is called from a vapply function to iterate / step chains in an MCMC
run2_mcmc <- function(.,j) {
  # runs in serial as each step depends on the previous step
  # call runp_mcmc

  # generate proposal matrix
  .$proposal_generate(j=j)
  .$boundary_handling()

  # evaluate model for proposal on each chain
  .$dataf$out[]  <-
    do.call( 'rbind', {
        if(.$wpars$multic) mclapply(1:.$dataf$lp, .$run3, mc.cores=min(.$wpars$procs,.$dataf$lp), mc.preschedule=F )
        else                 lapply(1:.$dataf$lp, .$run3 )
    })


  # calculate likelihood of proposals on each chain, accept/reject proposals
  prop_lklihood <- .$proposal_lklihood()
  curr_lklihood <- .$dataf$pars_lklihood[,j-1]
  accept        <- .$proposal_accept(prop_lklihood, curr_lklihood )

  # update chains
  for(ii in 1:.$dataf$lp) {
    if(accept[ii]) {
      .$dataf$pars_array[,ii,j]   <- .$dataf$pars[,ii]
      .$dataf$pars_lklihood[ii,j] <-  prop_lklihood[ii]
      .$dataf$out_mcmc[,ii,j]     <- .$dataf$out[ii,] 
    } else {
      .$dataf$pars_array[,ii,j]   <- .$dataf$pars_array[,ii,j-1]
      .$dataf$pars_lklihood[ii,j] <- curr_lklihood[ii]
      .$dataf$out_mcmc[,ii,j]     <- .$dataf$out_mcmc[,ii,j-1]
    }

    # compute squared normalized jumping distance
    # - count selected crosover values
    # - jump_delta_norm = Delta_m V2011
    # - CR_counter = L_m V2011
    # APW: potentially inefficient, does this need calculating for each chain or can it be simultaneous?, fix
    if(.$wpars$mcmc$adapt_pCR) {
      summation <- sum(((.$dataf$pars_array[,ii,j] - .$mcmc$current_state[,ii]) / .$mcmc$sd_state)^2 )
      .$mcmc$jump_delta_norm[.$mcmc$CR[ii]] <- .$mcmc$jump_delta_norm[.$mcmc$CR[ii]] + summation
      .$mcmc$CR_counter[.$mcmc$CR[ii]]      <- .$mcmc$CR_counter[.$mcmc$CR[ii]] + 1
    } 
  }

  # update the selection probability of crossover probabilities/values
  if(!.$wpars$parsinit_read)
    if(j==.$wpars$mcmc$preburnin_iter & .$wpars$mcmc$adapt_pCR ) .$mcmc$adapt_pCR <- T
  if(.$mcmc$adapt_pCR) .$mcmc_adapt_pCR() 

  
  # append history matrix
  # - for dreamzs
  if(.$wpars$mcmc$mcmc_type=='dreamzs') {
    if((.$mcmc$j_true%%.$wpars$mcmc$iterappend)==0) {
      new_lps                                       <- .$dataf$lps + .$wpars$mcmc$chains
      .$dataf$past_states[,(.$dataf$lps+1):new_lps] <- .$dataf$pars_array[,,j] 
      .$dataf$lps                                   <- new_lps 
    } 
  }


  # MCMC checks
  mcmc_check <- ((j>=.$wpars$mcmc$preburnin_iter)|.$wpars$parsinit_read) & (j%%.$wpars$mcmc$check_iter==0)
  if(mcmc_check | (j==.$wpars$mcmc$maxiter)) {

    # subscript for convergence etc arrays, and 50 % of current post-outlier samples in pars etc arrays
    .$mcmc$check_ss   <- .$mcmc$check_ss + 1 
    .$mcmc$j_burnin50 <-
      if(.$mcmc$outlier_detected | !.$wpars$parsinit_read) .$mcmc$j_start_burnin + ceiling((j-.$mcmc$j_start_burnin)/2)
      else                                                 j - ceiling((j+.$wpars$mcmc$start_iter-1)/2) + 1

    # test for and handle outlier chains 
    # - if mean of last conv_period convergence stat calculations for all parameters <1.2 no longer check outliers 
    if(.$mcmc$outlier_test) { 
      if(!is.null(.$wpars$mcmc$conv_period)) {
      if(.$mcmc$check_ss>.$wpars$mcmc$conv_period) {
        cc_ss <- (.$mcmc$check_ss-.$wpars$mcmc$conv_period):(.$mcmc$check_ss-1)
        if(!is.na(sum(.$dataf$conv_check[,cc_ss]))) {
          cc_mean <- apply(.$dataf$conv_check[,cc_ss], 1, mean, na.rm=T )
          if(sum(cc_mean<1.18) == dim(.$dataf$pars)[1]) {
            .$mcmc$outlier_test <- F 
            print('',quote=F); print('',quote=F)
            print(paste('Parameters converged over last',.$wpars$mcmc$check_iter*.$wpars$mcmc$conv_period,'iterations:'),quote=F); print(cc_mean,quote=F)
          }
        }
      }}
      outliers <- .$mcmc_outlier(j=j)
      if(length(outliers)>0) .$mcmc_outlier_handling(outliers, j ) 
    }

    # test for convergence
    R_hat <- .$mcmc_converge(j=j)
    .$dataf$conv_check[,.$mcmc$check_ss] <- c(.$mcmc$j_true, j, R_hat )
    if(j==.$wpars$mcmc$maxiter) .$mcmc_handle_iter_final(R_hat)
  }


  # update MCMC iteration counter
  .$mcmc$j_true[] <- .$mcmc$j_true + 1

  # return nothing - allows use of the stable vapply to call this function
  numeric(0)
}


# This wrapper function is called from an lapply or mclappy function to pass every column of the dataf$pars matrix to the model
run3_mcmc <- function(.,k) {
  # runs each chain at each iteration in MCMC

  # assumes that columns of the pars matrix (i.e. chains) are independent and non-sequential

  # configure parameters in the model
  if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[,k], F )
  if(.$wpars$cverbose)       .$printc('pars', .$dataf$pars[,k] )

  # call model/met run function
  if(is.null(.$dataf$met)) .$model$run() else .$model$run_met()
}

run4_mcmc <- run4_factorial # i.e. NULL
run5_mcmc <- run5_factorial # i.e. NULL
run6_mcmc <- run6_factorial # i.e. NULL
run7_mcmc <- run7_factorial # i.e. NULL
run8_mcmc <- run8_factorial # i.e. NULL



# Output processing functions
###########################################################################

################################
# write output

write_output_factorial <- function(.) {
  .$wpars$of_name <- .$wpars$of_name_stem
  .$write_to_file()
}


# write AB output array
write_output_SApar_saltelli <- function(.) {

  #write_to_file(.$output_saltelli_AB(), paste(ofname,'salt','AB',sep='_'), type='rds' )
  .$wpars$of_name <- paste(.$wpars$of_name_stem,'salt','AB',sep='_')
  .$write_to_file()

  # write dataf matrices used in AB run
  .$wpars$of_name <- paste(.$wpars$of_name_stem,'salt','AB','dataf',sep='_')
  .$write_to_file(df=list(fnames=.$dataf$fnames, pars=.$dataf$pars, env=.$dataf$env ))

  # remove large out array
  if(!.$wpars$unit_testing) .$dataf$out <- NULL

}


# write ABi output array
write_output_SApar_saltelli_ABi <- function(.) {
  .$wpars$of_name <- paste(.$wpars$of_name_stem,'salt','ABi',sep='_')
  .$write_to_file()

  # remove large out array
  if(!.$wpars$unit_testing) .$dataf$out_saltelli <- matrix(1)
}


# write Ye output array
write_output_SAprocess_ye <- function(.,f) {
  # similar to init matrix Ye method does not output for whole simulation but rather for each process
  .$wpars$of_name <- paste(.$wpars$of_name_stem, 'proc', f, sep='_' )
  .$write_to_file()
}


# write MCMC output list
write_output_mcmc <- function(.,i) {

  # write MCMC output
  #.$wpars$of_name <- paste(ofname, 'mcmc', 'f', i, sep='_' )
  .$write_to_file()

  # determine restart iteration number 
  hist_file_list <- list.files(pattern=paste0(.$wpars$of_name, '_history_' ))      
  if(length(hist_file_list)==0) {
    hn <- 1
  } else {
    hstart <- regexpr('history_', hist_file_list )
    hend   <- regexpr('.RDS', hist_file_list )
    hns    <- as.numeric(substr(hist_file_list, hstart[1]+attr(hstart,'match.length')[1], hend-1 ))
    hn     <- max(hns) + 1
  } 

  # make and write history file 
  df <- .$output(iter_out_start=1, iter_out_end=.$mcmc$j_burnin50-1 )
  saveRDS(df, paste0(.$wpars$of_name, '_history_', hn, '.RDS' ))
}



################################
# generate output

# function that combines the "vars", "met", and "out" dataframes correctly for output in a factorial simulation
# - not sure how much this is actually used, seems to be just for a factorial run these days
output_factorial  <- function(.){
  return(
    # if at least one of fnames, pars, and env are varied
    if(is.null(.$dataf$env)+is.null(.$dataf$pars)+is.null(.$dataf$fnames) < 3) {

      vpars   <- if(is.null(.$dataf$pars))    NULL else .$dynamic$pars
      venv    <- if(is.null(.$dataf$env))     NULL else .$dynamic$env
      vfnames <- if(is.null(.$dataf$fnames))  NULL else .$dynamic$fnames
      vardf   <- expand.grid(c(venv,vpars,vfnames), stringsAsFactors=F )

      #print('')
      #print('output processing')
      #print(.$dataf$pars)
      #print(vpars)

      # if no met data
      if(is.null(.$dataf$met)) {
        if(.$wpars$UQ) {
          # return a list
          return(c(vardf, list(out=.$dataf$out) ))
          rm(vardf)
        } else {
          # return a dataframe
          return(cbind(vardf, .$dataf$out ))
          rm(vardf)
        }

      # if met data
      # - so far will only work for factorial simulations
      } else {
        odf <- cbind(do.call(rbind, lapply(1:length(vardf[,1]), .$combine, df=vardf )), .$dataf$out )
        if(dim(vardf)[2]==1) names(odf)[which(names(odf)=='df.i...')] <- names(vardf)
        rm(vardf)
        return(odf)
        rm(odf)
      }

    # if no vars
    } else {
      # if met data
      if(!is.null(.$dataf$met)) cbind(t(.$dataf$met), .$dataf$out ) else .$dataf$out
    }
  )
}


# creates output for a saltelli Sobol sensitivity analysis
output_SApar_saltelli <- function(.) {
  # A and B matrices are stacked in a single matrix, which for each model and environment combination are then stored in an array

  # AB output is an array
  # - dim 1 (rows)   model combination
  # - dim 2 (cols)   environment combination
  # - dim 3 (slices) sample
  # - dim 4          output variable (character variables are coerced to NAs)

  # create AB output matrix array
  AB  <- array(.$dataf$out, c(.$dataf$le, 2*.$wpars$n, .$dataf$lf, length(.$dataf$mout) ))
  dimnames(AB) <- list( NULL, NULL, apply(.$dataf$fnames, 2, toString), names(.$dataf$mout)  )

  # output a list composed of the AB matrix output array, the fnames that define each model combination, the parameter names
  aperm(AB, c(3,1,2,4) )
}


# creates output for a Saltelli Sobol sensitivity analysis
output_SApar_saltelli_ABi <- function(.) {

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

  aperm(.$dataf$out_saltelli, c(5,4,2,1,3) )
}


# creates output for a Ye process sensitivity analysis
output_SAprocess_ye <- function(.) {
  .$dataf$out
}


# creates output for an MCMC simulation
output_mcmc <- function(., iter_out_start=.$mcmc$j_burnin50, 
                        iter_out_end=.$wpars$mcmc$maxiter ) {
  # print output  
  pa_dim <- dim(.$dataf$pars_array)
  cc_dim <- dim(.$dataf$conv_check)
  df1    <- data.frame(convstat=as.vector(t(.$dataf$conv_check[4:cc_dim[1],])), iter=.$dataf$conv_check[1,], pars=rep(letters[1:pa_dim[1]],each=cc_dim[2]) )
  df2    <- data.frame(omega=as.vector(t(.$dataf$omega[1:pa_dim[2],])), iter=.$dataf$conv_check[1,], chain=rep(letters[1:pa_dim[2]],each=dim(.$dataf$omega)[2]) )
  # iter_effective
  p1 <-
    xyplot(.$dataf$conv_check[3,]~.$dataf$conv_check[1,], df2, type='l', 
           abline=list(v=.$dataf$conv_check[1,is.na(.$dataf$conv_check[4,])],h=0), 
           ylab='Effective iteration', xlab='iteration (vert. line = outlier detected)' )
  # convergence check
  p2 <-
  xyplot(omega~iter, df2, groups=chain,type='l', 
         abline=list(v=.$dataf$conv_check[1,is.na(.$dataf$conv_check[4,])]), ylab='Outlier stat.', xlab='iteration'  )
  p3 <-
  xyplot(convstat~iter, df1, groups=pars, type='l', 
         abline=list(h=c(1.1,1.2)), ylab='Convergence stat.', xlab='iteration' )
  
  # parameter histograms
  mcmc_pa <- array(.$dataf$pars_array, dim=c(pa_dim[1],prod(pa_dim[2:3])) )
  mcmc_pa <- data.frame(par_val=as.vector(t(mcmc_pa)), par_name=rep(dimnames(.$dataf$pars_array)[[1]], each=prod(pa_dim[2:3]) ))
  p4      <- histogram(~par_val|par_name, mcmc_pa, breaks=3e1, scales='free', xlab='', as.table=T )
  pdf(paste0(.$wpars$of_name,'_stats.pdf'))
  print(p1, split=c(1,1,1,3), more=T )
  print(p2, split=c(1,2,1,3), more=T )
  print(p3, split=c(1,3,1,3), more=F )
  print(p4, split=c(1,1,1,1), more=F )
  dev.off()


  # combine output
  list(
    pars_array     = .$dataf$pars_array[,,iter_out_start:iter_out_end,drop=F],
    pars_lklihood  = .$dataf$pars_lklihood[,iter_out_start:iter_out_end,drop=F],
    out_mcmc       = .$dataf$out_mcmc[,,iter_out_start:iter_out_end,drop=F],
    past_states    = .$dataf$past_states,
    mod_out_final  = .$dataf$out,
    conv_check     = .$dataf$conv_check,
    omega          = .$dataf$omega,
    met            = .$dataf$met,
    obs            = .$dataf$obs,
    obsse          = .$dataf$obsse,

    # user defined MCMC parameters
    wpars = list(
      mcmc = .$wpars$mcmc
    ),

    # dynamic MCMC variables that need preserving across a restart
    mcmc = list(
      prior_sample    = .$mcmc$prior_sample,
      p_CR            = .$mcmc$p_CR,
      jump_delta_norm = .$mcmc$jump_delta_norm,
      CR_counter      = .$mcmc$CR_counter,
      j_true          = .$mcmc$j_true,
      adapt_pCR       = .$mcmc$adapt_pCR,
      outlier_test    = .$mcmc$outlier_test,
      obs_vars        = .$mcmc$obs_vars
    )
  )
}



### END ###
