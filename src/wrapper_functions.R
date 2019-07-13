################################
#
# Wrapper object for model objects 
# 
# AWalker December 2015
#
################################

library(parallel)

source('functions/general_functions.R')
source('functions/calc_functions.R')




#####################################

# expand the fnames, env, and pars input lists into matrices   
# - each row is passed to the model sequentially by each run function in the run function cascade
generate_ensemble <- function(.) {

  .$dataf$fnames  <- if(!is.null(.$dynamic$fnames)) as.matrix(expand.grid(.$dynamic$fnames,stringsAsFactors=F)) else NULL
  .$dataf$env     <- if(!is.null(.$dynamic$env))    as.matrix(expand.grid(.$dynamic$env,stringsAsFactors=F   )) else NULL
  .$generate_ensemble_pars()
  
  # calculate input matrix lengths - separate out as a function 
  # - used to set the number of iterations in the run functions  
  # - if no matrix return 1
  if(.$wpars$UQ&.$wpars$runtype=='SAprocess_ye') {
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
}


# parameter matrix for factorial run  
generate_ensemble_pars_factorial <- function(.) {
  .$dataf$pars <- if(!is.null(.$dynamic$pars)) as.matrix(expand.grid(.$dynamic$pars,stringsAsFactors=F)) else NULL
} 


# parameter matrix for Saltelli parameter SA  
generate_ensemble_pars_SApar_saltelli <- function(.) {
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
  .$dataf$pars   <- do.call(cbind, .$dynamic$pars )
  
  # remove potentially large pars list 
  .$dynamic$pars <- lapply(.$dynamic$pars, function(e) numeric(1) )        
}


# parameter matrix for Dai, Ye, process SA
# - the paramter matrices for the Ye method are generated in run function 1
# - here a list structure in .$dynamic$pars is created from .$dynamic$pars_eval
generate_ensemble_pars_SAprocess_ye <- function(.) {
  # need a minimum of >1 processes
  if(dim(.$dataf$fnames)[2]<=1) stop('need more than one process for a process sesitivity analysis')
  
  # check input dynamic$pars* are same length
  test_in <- length(.$dynamic$pars_eval) - length(.$dynamic$pars_proc)
  if(test_in!=0) stop('wrapper: Parameter input vectors - pars_eval & pars_proc - are not the same length')
  
  # assign same list structure as dynamic$pars_eval to vars$pars 
  .$dynamic$pars <- lapply(.$dynamic$pars_eval,function(e) numeric(1) )
  
  # check input vars$pars* elements have same names
  # - to be done
}

# debug: set seed functions (to reproduce sequences of quasi-random numbers)
################################

# debug: funciton (1), set seed for uniform_r generation
# set_seed1 <- function(.) {
#  set.seed(1703)
#  .$mcmc$uniform_r_seed <- matrix(data = 0, nrow = .$wpars$mcmc_maxiter, ncol = .$dataf$lp)
#  .$mcmc$uniform_r_seed[] <- runif((.$wpars$mcmc_maxiter * .$dataf$lp), min = -0.01, max = 0.01)
# }


# debug: function (2), set seed for R1 and R2 general_functions
# set_seed2 <- function(.) {
#   set.seed(4050)
# }


# debug: function (3), set seed for runif(1) value chosen in accept/reject step
# set_seed3 <- function(.) {
#   set.seed(1337)
#   .$mcmc$runif_seed <- matrix(data = 0, nrow = .$wpars$mcmc_maxiter, ncol = .$dataf$lp)
#   .$mcmc$runif_seed[] <- runif((.$wpars$mcmc_maxiter * .$dataf$lp), min = 0, max = 1)
# }

# MCMC functions
################################

# set parameter boundaries from prior distributions
boundary_handling_set <- function(.) {
  # number of samples to be used in boundary handling
  n <- 1000
  boundary_sample <- lapply(.$dynamic$pars_eval, function(cs) eval(parse(text=cs)) )
  .$mcmc$boundary_min <- unlist(lapply(boundary_sample,min))
  .$mcmc$boundary_max <- unlist(lapply(boundary_sample,max))
  rm(boundary_sample)
}


# restrict parameter proposals that are beyond the boundary
boundary_handling <- function(., ii, jj ) {
  if      (.$dataf$pars[ii,jj] < .$mcmc$boundary_min[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_min[jj]
  else if (.$dataf$pars[ii,jj] > .$mcmc$boundary_max[jj]) .$dataf$pars[ii,jj] <- .$mcmc$boundary_max[jj]
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
  # initialise output array
  # - dim 1 (rows)      output variable
  # - dim 2 (columns)   sample
  # - dim 3 (slices)    parameter that has used value from matrix B while all other par values are from matrix A 
  # - dim 4 (cube rows) environment combination
  # - dim 5 (cube cols) model combination
  .$dataf$out_saltelli <- array(0, dim=c(length(.$dataf$mout), .$wpars$n, dim(.$dataf$pars)[2], .$dataf$le, .$dataf$lf ))
  dimnames(.$dataf$out_saltelli) <- list(names(.$dataf$mout), NULL, colnames(.$dataf$pars), NULL, apply(.$dataf$fnames, 1, toString) )
}


init_output_matrix_SAprocess_ye <- function(.) {
  # Ye method does not generate a single for whole simulation but rather a separate ensemble for each process 
  # - dim 1 (rows)        output variable
  # - dim 2 (columns)     environment combination
  # - dim 3 (slices)      process(es) B parameter sample 
  # - dim 4 (cube rows)   process(es) B representation(s)
  # - dim 5 (cube cols)   process A parameter sample
  # - dim 6 (cube slices) process A representation
  # if met data then ... .$dataf$out     <- array(0, c(length(.$dataf$mout), .$dataf$lm, .$dataf$le, .$wpars$n, .$dataf$lfB, .$wpars$n, .$dataf$lfA  ) )
  .$dataf$out           <- array(0, c(length(.$dataf$mout), .$dataf$le, .$wpars$n, .$dataf$lfB, .$wpars$n, .$dataf$lfA  ) )
  dimnames(.$dataf$out) <- list(names(.$dataf$mout), NULL, NULL, apply(.$dataf$fnamesB, 1, toString), NULL, .$dataf$fnames )
  
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
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
  # assumes that each row of the fnames matrix are independent and non-sequential
  # call run2

  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[i,], F )
  if(.$wpars$cverbose)         .$printc('fnames', .$dataf$fnames[i,] )

  # call next run function
  do.call( 'rbind', {
      if(.$wpars$multic) mclapply(1:.$dataf$lp, .$run2, mc.cores=max(1,floor(.$wpars$procs/.$dataf$lf)), mc.preschedule=T  )
      else                 lapply(1:.$dataf$lp, .$run2 )
  })
}


# run2
run2_factorial <- function(.,j) {
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
  # assumes that each row of the pars matrix are independent and non-sequential
  # call run3
  
  # configure parameters in the model
  if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[j,], F )
  if(.$wpars$cverbose)       .$printc('pars', .$dataf$pars[j,] )
  
  # call next run function
  funv   <- if(is.null(.$dataf$met)) .$dataf$mout else array(0, dim=c(.$dataf$lm, length(.$dataf$mout) ) )
  out    <- vapply(1:.$dataf$le, .$run3, funv )

  # out has the potential to be a vector, matrix (needs transposed), or an array (needs stacking)
  # returns matrix
  if(class(out)=='matrix') t(out) else if(class(out)=='array') .$stack(out) else as.matrix(out)
}


# run3
run3_factorial <- function(.,k) {
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
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
  # assumes that each row of the fnames matrix are independent and non-sequential
  # call run6
  
  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames',df=.$dataf$fnames[i,],F)
  if(.$wpars$cverbose) .$printc('fnames',.$dataf$fnames[i,])
  
  # call next run function
  vapply({
    if(.$wpars$multic) mclapply(1:.$dataf$le, .$run6, mc.cores=min(.$dataf$le,.$wpars$procs), mc.preschedule=F ) 
    else                 lapply(1:.$dataf$le, .$run6 )
  }, function(a) a, .$dataf$out_saltelli[,,,1,1] )
}


run6_SApar_saltelli <- function(.,k) {
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$env matrix to the model
  # assumes that each row of the fnames matrix are independent and non-sequential
  # call run7
  
  # configure environment in the model
  if(!is.null(.$dataf$env)) .$model$configure(vlist='env', df=.$dataf$env[k,], F )
  if(.$wpars$cverbose)      .$printc('env', .$dataf$env[k,] )
  
  # call parameter matrix run function
  if(is.null(.$dataf$met)){
    
    # call next run function, wrapped within vapply to convert (mc)lapply list output to an array
    # returns a numeric array - model output variable (rows), sample (columns), parameter (slices)
    vapply({
      if(.$wpars$multic) mclapply(1:dim(.$dataf$pars)[2], .$run7, mc.cores=max(1,floor(.$wpars$procs/.$dataf$le)), mc.preschedule=T ) 
      else                 lapply(1:dim(.$dataf$pars)[2], .$run7 )
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
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnames matrix to the model
  # assumes that each row of the fnames matrix are independent and non-sequential
  # call run3
 
  print(paste('started representation:', .$dataf$fnames[g,], ', of process:', colnames(.$dataf$fnames)), quote=F )
 
  # configure function names in the model
  if(!is.null(.$dataf$fnames)) .$model$configure(vlist='fnames', df=.$dataf$fnames[g,] , F )
  if(.$wpars$cverbose) .$printc('fnames', .$dataf$fnames[g,] )

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
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$pars matrix to the model
  # assumes that each row of the pars matrix are independent and non-sequential
  # call run4
  
  # configure parameters in the model
  if(!is.null(.$dataf$pars)) .$model$configure(vlist='pars', df=.$dataf$pars[h,], F )
  if(.$wpars$cverbose) .$printc('pars', .$dataf$pars[h,] )
  
  # calculate offset to correctly subset parsB matrix
  osh  <- offset + .$dataf$lfB * (h-1)        
  
  # call process B process representation run function
  vapply(1:.$dataf$lfB, .$run4, .$dataf$out[,,,1,1,1], offset=osh )
}


# run4
run4_SAprocess_ye <- function(., i, offset ) {
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$fnamesB matrix to the model
  # assumes that each row of the fnamesB matrix are independent and non-sequential
  # call run5
  
  # configure function names in the model
  if(!is.null(.$dataf$fnamesB)) .$model$configure(vlist='fnames', df=.$dataf$fnamesB[i,], F )
  if(.$wpars$cverbose) .$printc('fnames', .$dataf$fnamesB[i,] )
  
  # calculate offset to correctly subset parsB matrix
  os  <- offset + i
  # oss is a vector of the row subscripts for the parsB matrix
  oss <- (.$wpars$n*(os-1) + 1):(.$wpars$n*(os))    
  
  # call process B parameter run function
  vapply(oss, .$run5, .$dataf$out[,,1,1,1,1] )
}


# run5
run5_SAprocess_ye <- function(.,j) {
  # This wrapper function is called from an lapply or mclappy function to pass every row of the dataf$parsB matrix to the model
  # assumes that each row of the parsB matrix are independent and non-sequential
  # call run6
  
  # configure parameters in the model
  if(!is.null(.$dataf$parsB)) .$model$configure(vlist='pars', df=.$dataf$parsB[j,], F )
  if(.$wpars$cverbose)        .$printc('pars', .$dataf$parsB[j,] )
  
  # call the environment run function
  vapply(1:.$dataf$le, .$run6, .$dataf$out[,1,1,1,1,1] )
}
    

run6_SAprocess_ye <- run3_factorial # i.e. over .$dataf$env
run7_SAprocess_ye <- run4_factorial # i.e. NULL 
run8_SAprocess_ye <- run4_factorial # i.e. NULL 



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


################################
# generate output

# function that combines the "vars", "met", and "out" dataframes correctly for output in a factorial simulation 
# - not sure how much this is actually used, seems to be just for a factorial run these days
output_factorial  <- function(.){
  return(
    # if at least one of fnames, pars, and env are varied
    if(is.null(.$dataf$env)+is.null(.$dataf$pars)+is.null(.$dataf$fnames) < 3) {

      vpars    <- if(is.null(.$dataf$pars))    NULL else .$dynamic$pars
      venv     <- if(is.null(.$dataf$env))     NULL else .$dynamic$env
      vfnames  <- if(is.null(.$dataf$fnames))  NULL else .$dynamic$fnames
      vardf <- expand.grid(c(venv,vpars,vfnames),stringsAsFactors=F)        
      
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
        odf <- cbind(do.call(rbind , lapply(1:length(vardf[,1]) , .$combine, df=vardf ) ), .$dataf$out )            
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
  dimnames(AB) <- list( NULL, NULL, apply(.$dataf$fnames, 1, toString), names(.$dataf$mout)  )
  
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

# DEMC functions
################################

# generate proposal using DE-MC algorithm
proposal_generate_demc <- function(., j ) {

  # Metropolis sampling

  # scaling factor
  d <- ncol(.$dataf$pars)
  gamma_star <- 2.38 / sqrt(d + d)

  # b-value should be small compared to width of target distribution
  # b_rand specifies range for drawn "randomization" value
  b_rand  <- 0.01

  # debug: temporarily hardcode uniform_r
  # uniform_r <- rep(-0.000378, d)

  uniform_r <- runif(1, min = (-b_rand), max = b_rand)

  # debug: (1) set the seed for uniform_r generation
  # uniform_r <- .$mcmc$uniform_r_seed[j, ii]
  # uniform_r <- .$mcmc$uniform_r_seed[j, 1:.$dataf$lp]
  # uniform_r <- rep(.$mcmc$uniform_r_seed[j, ii], d)

  # print("uniform_r = ")
  # print(uniform_r)

  # evaluate for each chain
  for (ii in 1:.$dataf$lp) {

    # debug: (1) set the seed for uniform_r generation
    # uniform_r <- rep(.$mcmc$uniform_r_seed[j, ii], d)

    # debug: R1 and R2 iniitalization to 0 inside the for-loop
    # randomly select two different numbers R1 and R2 unequal to j
    # from a uniform distribution without replacement
    R1 <- 0
    R2 <- 0

    while ((R1 == 0) | (R1 == ii))               R1 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)
    while ((R2 == 0) | (R2 == ii) | (R2 == R1))  R2 <- ceiling(runif(1,min=0,max=1)*.$dataf$lp)

    # debug: temporarily hardcode R1 and R2
    # R1 <- 6
    # R2 <- 5

    # print(paste0('<<<< iteration = ', j, ', chain = ', ii, ' <<<<'))

    # print('.$dataf$pars_array = ')
    # print(.$dataf$pars_array)

    # debug: check dimensions
    # print(paste0('dim of .$dataf$pars_array = ', dim(.$dataf$pars_array)))
    # print(paste0('dim of .$dataf$pars = ', dim(.$dataf$pars)))

    # print(paste0("R1 = ", R1, ", R2 = ", R2))

    # debug: store R1 and R2 values that were generated
    # .$dataf$R1_R2_storage[ii, 1, j] <- R1
    # .$dataf$R1_R2_storage[ii, 2, j] <- R2

    # evaluate for each parameter
    for (jj in 1:d) {

      # generate proposal via Differential Evolution

      # debug: restructure to make setting the seed easier
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[j]

      # debug: note that this one below was the original at beginning of bug-fixing process
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[jj]

      # debug: create "scalar" to measure jump differential
      # jump <- gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r[jj]
      # .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + jump

      # store uniform_r value to check that set.seed() code is functioning
      # .$dataf$uniform_r_storage[ii, jj, j] <- uniform_r[jj]

      # print(paste0('ii = ', ii))
      # print(paste0('jj = ', jj))
      # print(paste0('j = ', j))

      # print(paste0('gamma_star = ', gamma_star))

      # print(paste0('uniform_r[jj] = ', uniform_r[jj]))

      # print(paste0('jump = ', jump))

      # print(paste0('.$dataf$pars_array[R1, jj, j-1] = ', .$dataf$pars_array[R1, jj, j-1]))
      # print(paste0('.$dataf$pars_array[R2, jj, j-1] = ', .$dataf$pars_array[R2, jj, j-1]))
      # print(paste0('.$dataf$pars_array[ii, jj, j-1] = ', .$dataf$pars_array[ii, jj, j-1]))
      # print(paste0('.$dataf$pars[ii, jj] = ', .$dataf$pars[ii, jj]))

      .$dataf$pars[ii,jj] <- .$dataf$pars_array[ii,jj,j-1] + gamma_star * (.$dataf$pars_array[R1,jj,j-1] - .$dataf$pars_array[R2,jj,j-1]) + uniform_r

      # print(paste0('uniform_r = ', uniform_r[jj]))

      # debug: temporarily comment out boundary handling
      # call boundary handling function
      boundary_handling(., ii, jj )
    }

    # print('proposal generated = ')
    # print(.$dataf$pars[ii, ])

    # store proposals being generated (regardless of whether or not they are being accepted)
    # .$dataf$prop_storage[ii, 1:d, j] <- .$dataf$pars[ii, 1:d]

  }

  # debug: print statements at the end of each run
  # if (j == .$wpars$mcmc_maxiter) {
    # print('array of proposals that were generated = ')
    # print(.$dataf$prop_storage)
    # print('uniform_r values (randomly generated) = ')
    # print(.$dataf$uniform_r_storage)
    # print('R1 (1st col.) and R2 (2nd col.) values (randomly generated) = ')
    # print(.$dataf$R1_R2_storage)
  # }

}


# calculate proposal acceptance using the Metropolis ratio (for DE-MC algorithm)
proposal_accept_demc <- function(., j, lklihood ) {

  # Metropolis ratio
  metrop_ratio <- exp(lklihood - .$dataf$pars_lklihood[ ,j-1])

  .$dataf$metrop_ratio_storage[1:.$wpars$mcmc_chains, 1, j] <- t(metrop_ratio)

  # print(paste0('likelihood of proposal = ', lklihood))

  # print(paste0('likelihood of current chain = ', .$dataf$pars_lklihood[ ,j-1]))

  # print(paste0('dim of .$dataf$pars_lklihood ' = dim(.$dataf$pars_lklihood)))

  # print(paste0('Metropolis ratio = ', metrop_ratio))
  # print(length(metrop_ratio))

  alpha        <- pmin(1, metrop_ratio)

  .$dataf$alpha_storage[1:.$wpars$mcmc_chains, 1, j] <- t(alpha)

  # debug: maybe try computing alpha this way (more numerically comprehensive/rigorous)
  # alpha <- numeric(.$dataf$lp)
  # for (q in 1:.$dataf$lp) {
  #   if (.$dataf$pars_lklihood[q, j-1] > 0) {
  #     alpha[q] <- min(1, metrop_ratio[q])
  #   } else {
  #     alpha[q] <- 1
  #   }
  # }

  # print('alpha = ')
  # print(alpha)

  # debug: (3) set the seed for runif(1) value chosen in accept/reject step
  # runif_val <- .$mcmc$runif_seed[j, ii]
  # runif_val <- .$mcmc$runif_seed[j, 1:.$dataf$lp]

  # print('runif_val = ')
  # print(runif_val)

  # debug: store runif(1) value
  # .$dataf$runif_val_storage[1:.$wpars$mcmc_chains, 1, j] <- t(runif_val)

  # print(paste0('length of alpha = ', length(alpha)))
  # print(paste0('type of alpha = ', typeof(alpha)))

  # print(paste0('length of runif_val = ', length(runif_val)))
  # print(paste0('type of runif_val = ', typeof(runif_val)))

  # evaluate for each chain
  for(ii in 1:.$dataf$lp) {

    # print(paste0('<<<< iteration = ', j, ', chain = ', ii, ' <<<<'))

    # accept if Metropolis ratio > random number from uniform distribution on interval (0,1)
    accept <- log(alpha[ii]) > log(runif(1, min = 0, max = 1))

    # debug: temporarily hard-code runif-value for debugging
    # accept <- log(alpha[ii]) > log(0.843332)

    # debug:(3) set the seed for runif(1) value chosen in accept/reject step
    accept <- log(alpha[ii]) > log(runif_val[ii])

    # debug: store accept value
    # .$dataf$accept_storage[ii, 1, j] <- accept

    # print(paste0('runif_val = ', runif_val[ii]))

    # print(paste0('length of accept = ', length(accept)))

    # print('accept = ')
    # print(accept)

    .$dataf$pars_array[ii,,j]   <- if(accept) .$dataf$pars[ii,] else .$dataf$pars_array[ii,,j-1]
    .$dataf$pars_lklihood[ii,j] <- if(accept) lklihood[ii]      else .$dataf$pars_lklihood[ii,j-1]


    # debug: store every model evaluation debugging purposes
    # out_n <- .$wpars$mcmc_maxiter / 2
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
    # APW: not a huge deal, but this is not quite right in the case where j==out_n+1 but not accepted as it adds the output from the rejected proposal
    # if ((j+1) > out_n)
    #   .$dataf$out_mcmc[ii,,(j+1-out_n)] <- if(accept | j==out_n+1) .$dataf$out[ii,] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }

  # debug: print statements
  # if (j == .$wpars$mcmc_maxiter) {
    # print('Metropolis ratio = ')
    # print(.$dataf$metrop_ratio_storage)
    # print('alpha = ')
    # print(.$dataf$alpha_storage)
    # print('randomly generated runif_val = ')
    # print(.$dataf$runif_val_storage)
    # print('accept/reject values (1 = TRUE and 0 = FALSE) = ')
    # print(.$dataf$accept_storage)
  # }

}



# DREAM functions
################################

# static part of DREAM algorithm
static_dream <- function(.) {

  # number of parameters being estimated
  .$mcmc$d <- ncol(.$dataf$pars)

  # preallocate memory space for algorithmic variables
  .$mcmc$J             <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$n_id          <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$CR            <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$p_CR          <- numeric(.$wpars$mcmc_n_CR)
  .$mcmc$R             <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$dataf$lp - 1)
  .$mcmc$current_state <- matrix(data = 0, nrow = .$dataf$lp, ncol = .$mcmc$d)
  .$mcmc$p_state       <- numeric(.$dataf$lp)
  .$mcmc$sd_state      <- numeric(.$mcmc$d)
  .$mcmc$jump          <- matrix(data=0, nrow=.$dataf$lp,   ncol=.$mcmc$d )
  .$mcmc$draw          <- matrix(data=0, nrow=.$dataf$lp-1, ncol=.$dataf$lp )
  .$mcmc$lambda        <- matrix(data=0, nrow=.$dataf$lp,   ncol=1 )

  # index of chains for Differential Evolution
  for (ii in 1:.$dataf$lp) .$mcmc$R[ii, ] <- setdiff(1:.$dataf$lp, ii)

  # crossover values
  .$mcmc$CR[] <- 1:.$wpars$mcmc_n_CR / .$wpars$mcmc_n_CR

  # selection probability of crossover values
  .$mcmc$p_CR[] <- 1 / .$wpars$mcmc_n_CR

  # vector that stores how many times crossover value indices are used
  # ALJ: initialized to 1's in order to avoid numeric issues
  # debug: maybe noteworthy that this was originally initialized to 0's in Vrugt's algorithm
  # debug: play around with whether or not this makes a difference
  .$mcmc$n_id[] <- 1
}


# generate proposal using DREAM algorithm
proposal_generate_dream <- function(., j ) {

  # reset matrix of jump vectors to zero
  .$mcmc$jump[] <- 0

  # current state, 'mcmc_chains' number of samples of a d-variate distribution
  # APW: I think you can just use the pars_array here and later
  # ALJ: you could, but I am hesitant to change it
  #      becasue using the current_state and jump matrices circumvents "function call order" issue we ran into with the DE-MC algorithm
  #      i.e., I think these "placeholder" matrices  are relevant to parallelization and/or other forms of the DREAM alg
  #      I would like to leave them until we are absolutely certain that the DREAM algorithm is up and running perfectly in MAAT
  # future work: play around with whether the current_state and jump matrices are absolutely essential
  .$mcmc$current_state[] <- matrix(.$dataf$pars_array[ , , j-1], nrow = .$dataf$lp, ncol = .$mcmc$d)

  # permute [1,2,...,mcmc_chains-1] mcmc_chains number of times
  .$mcmc$draw[]          <- apply(matrix(runif((.$dataf$lp - 1) * .$dataf$lp), .$dataf$lp - 1, .$dataf$lp), 2, function(v) sort(v, index.return = T)$ix)

  # create a .$dataf$lp x 1 matrix of continuous uniform random values between -c_rand and c_rand
  .$mcmc$lambda[]        <- matrix(runif(.$dataf$lp * 1, -.$wpars$mcmc_c_rand, .$wpars$mcmc_c_rand), .$dataf$lp)

  # compute standard deviation of each dimension (ie,compute standard deviation of each column of current_state matrix)
  .$mcmc$sd_state[]     <- apply(.$mcmc$current_state, 2, sd)

  # create proposals
  # future work: vectorize this for-loop to improve computational efficiency, but this is non-trivial
  for (ii in 1:.$dataf$lp) {

    # select delta (equal selection probability) (ie, choose 1 value from the vector [1:delta] with replacement)
    D <- sample(1:.$wpars$mcmc_delta, 1, replace = T)

    # extract vectors a and b not equal to ii
    a <- .$mcmc$R[ii, .$mcmc$draw[1:D, ii]]
    b <- .$mcmc$R[ii, .$mcmc$draw[(D + 1):(2 * D), ii]]

    # select index of crossover value (weighted sample with replacement)
    .$mcmc$id <- sample(1:.$wpars$mcmc_n_CR, 1, replace = T, prob = .$mcmc$p_CR)

    # draw d values from uniform distribution between 0 and 1
    zz <- runif(.$mcmc$d)

    # derive subset A of selected dimensions
    A  <- which(zz < .$mcmc$CR[.$mcmc$id])

    #  how many dimensions are sampled
    d_star <- length(A)

    # make sure that A contains at least one value
    if (d_star == 0) {
      A <- which.min(zz)
      d_star <- 1
    }

    # calculate jump rate
    gamma_d <- 2.38 / sqrt(2 * D * d_star)

    # select gamma
    # future work: figure out a way to make this chunk of code prettier
    temp1 <- c(gamma_d, 1)
    temp2 <- c(1 - .$wpars$mcmc_p_gamma, .$wpars$mcmc_p_gamma)
    gamma <- sample(temp1, 1, replace = T, prob = temp2)

    # compute jump differential evolution of ii-th chain
    .$mcmc$jump[ii,A] <- .$wpars$mcmc_c_ergod * rnorm(d_star) + (1 + .$mcmc$lambda[ii]) * gamma * sum((.$mcmc$current_state[a,A] - .$mcmc$current_state[b,A]), dim = 1)

    # compute proposal of ii-th chain
    .$dataf$pars[ii,1:.$mcmc$d] <- .$mcmc$current_state[ii,1:.$mcmc$d] + .$mcmc$jump[ii,1:.$mcmc$d]

    # debug: temporarily comment out boundary handling
    # call boundary handling function
    for (jj in 1:.$mcmc$d) boundary_handling(., ii, jj )
  }
}


# proposal acceptance function for the DREAM algorithm
proposal_accept_dream <- function(., j, lklihood) {

  # likelihood of current state
  # APW: is this assignment totally necessary? Could we just use the pars_lklihood matrix?
  # ALJ: same reasoning as above
  # future work: play around with whether p_state is absolutely essential
  .$mcmc$p_state[] <- .$dataf$pars_lklihood[ ,j-1]

  for (ii in 1:.$dataf$lp) {

    # compute Metropolis acceptance probability
    alpha <- min(1, exp(lklihood[ii] - .$mcmc$p_state[ii]))

    # future work: figure out how to make this clunky block of code prettier

    # determine if p_acc is larger than random number drawn from uniform distribution on interval [0,1]
    if (alpha > runif(1, min = 0, max = 1)) {

      # accept the proposal
      accept <- TRUE

      # APW: are these assignments totally necessary? Could we just use the pars & lklihood matrix / vector?
      # ALJ: no, I don't think all of these assignments are totally necessary
      #      but I left them all in to remain as true as possible to Vrugt's original pseudocode
      #      and because they make debugging easier
      #      once DREAM is integrated in the new MAAT (and we're confident that it's working), I can alter some of the unnecessary assignments
      .$mcmc$current_state[ii, 1:.$mcmc$d]  <- .$dataf$pars[ii,1:.$mcmc$d]
      .$mcmc$p_state[ii]                    <- lklihood[ii]

      # append accepted current_state and probability density to storage data frames
      .$dataf$pars_array[ii,1:.$mcmc$d,j]   <- .$mcmc$current_state[ii,1:.$mcmc$d]
      .$dataf$pars_lklihood[ii,j]           <- .$mcmc$p_state[ii]

      # debug: store generated proposal (regardless of whether accepted or not)
      # .$dataf$prop_storage[ii,1:.$mcmc$d,j] <- .$dataf$pars[ii,1:.$mcmc$d]

    } else {

      # reject the proposal
      accept <- FALSE

      # set jump back to zero for p_CR
      .$mcmc$jump[ii, 1:.$mcmc$d] <- 0

      # repeat previous current_state and probability density in storage data frames
      .$dataf$pars_array[ii, 1:.$mcmc$d,j]   <- .$dataf$pars_array[ii,1:.$mcmc$d,j-1]
      .$dataf$pars_lklihood[ii,j]            <- .$dataf$pars_lklihood[ii,j-1]

      # debug: store generated proposal (regardless of whether accepted or not)
      # .$dataf$prop_storage[ii, 1:.$mcmc$d,j] <- .$dataf$pars[ii, 1:.$mcmc$d]
    }

    # update jump distance crossover index
    # APW: this step is introducing NaNs when a value in sd_state is 0
    # ALJ: I think this issue is now resolved?
    .$mcmc$J[.$mcmc$id]    <- .$mcmc$J[.$mcmc$id] + sum((.$mcmc$jump[ii, 1:.$mcmc$d] / .$mcmc$sd_state)^2)

    # number of times index crossover is used
    .$mcmc$n_id[.$mcmc$id] <- .$mcmc$n_id[.$mcmc$id] + 1

    # debug: see above comments in analagous DE-MC function
    # out_n <- .$wpars$mcmc_maxiter / 2
    out_n <- 0
    if (j > out_n)
      .$dataf$out_mcmc[ii,,(j-out_n)] <- if(accept | j == out_n + 1) .$dataf$out[ii, ] else .$dataf$out_mcmc[ii,,(j-out_n-1)]
  }

  # print(.$mcmc$id)

  # print(.$mcmc$n_id)

  # print(sum(.$mcmc$jump[ii,]))

  # print(.$mcmc$sd_state)

  # print(.$mcmc$J)

  # print(.$mcmc$p_CR)

  # update selection probability of crossover
  # ALJ: altered original algorithm here to account for numerical issues
  # APW: still something odd going on here.
  #      with the linear example mcmc$J has two zeros in vector elements 2 & 3
  #      and then at iteration 54 and after the first element is getting a NaN
  # APW: Actually the behaviour is weirder than that, I've seen NaNs come in anywhere from iteration 12 and up
  #      this only happens when .$mcmc$id is 1,
  #      and for any full MCMC run using linear test, .$mcmc$id is always the same value
  # APW: this is probably my fault from moving the variables from .$mcmc to .$wpars but I can't see where
  #      seems to be because mcmc$p_CR always contains just 1 and 0's
  # APW: update, p_CR coverges to a 1 and 0's acfter the second proposal, perhaps this is expected behaviour but worth more investigation
  # debug: test less than versus greater than
  # debug: testing below - to prevent immediate transition to 1/0 states for p_CR, seems to work
  if ((j > (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
  # if ((j < (.$wpars$mcmc_maxiter / 10)) & (sum(.$mcmc$J) > 0)) {
    .$mcmc$p_CR <- .$mcmc$J    / .$mcmc$n_id
    .$mcmc$p_CR <- .$mcmc$p_CR / sum(.$mcmc$p_CR)
  }
}

# future work: function for detection and correction of outlier chains
# outlier_check <- function(.) {
# }

# future work: test for convergence using the R-statistic convergence diagnostic of Gelman and Rubin
# chain_converge <- function(.) {
# }

# Likelihood functions
################################

# expects model output to be probability - as in the output from the mixture model
f_proposal_lklihood_log <- function(.) {
  # derive log density   
  log(.$dataf$out)

  # print(paste0('model likelihood = ', log(.$dataf$out)))

  return(log(.$dataf$out))
}


# standard error probability density function with iid error residuals
f_proposal_lklihood_ssquared <- function(.) {

  # number of measured data points
  obs_n <- length(.$dataf$obs)

  # calculate error residual
  # - each chain is on rows of dataf$out take transpose
  error_residual_matrix <- t(.$dataf$out) - .$dataf$obs

  print(.$dataf$obs)

  # calculate sum of squared error
  # - error_residual_matrix now has chains on columns
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2) )

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2)*log(SSR)
}


# standard error probability density function with iid error residuals
# this function incorporates measurement errors
f_proposal_lklihood_ssquared_se <- function(.) {

  # read in measurement error and remove zeros from measurement error  
  sspos <- which(.$dataf$obsse > 1e-9)

  # number of measured data points (that do not have zero uncertainty)
  obs_n <- length(sspos)

  # observed error
  obsse <- if(.$wpars$mcmc_homosced)   rep(mean(.$dataf$obsse[sspos]), obs_n)
           else                .$dataf$obsse[sspos]

  # calculate error residual (each chain is on rows of dataf$out, take transpose)
  error_residual_matrix <- ( t(.$dataf$out)[sspos,] - .$dataf$obs[sspos] ) / obsse

  # calculate sum of squared error
  SSR <- apply(error_residual_matrix, 2, function(v) sum(v^2))

  # return log-likelihood vector corresponding to each chain/row in .$dataf$pars matrix
  -(obs_n/2)*log(2*pi) - sum(log(obsse)) - 0.5*SSR
}



### END ###
