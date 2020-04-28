###############################
#
# Wrapper object for model objects
#
# AWalker December 2015, June 2019
#
################################

library(proto)
source('wrapper_functions.R')
source('functions/general_functions.R')
source('functions/packagemod_functions.R')
source('functions/calc_functions.R')
source('wrapper_functions_mcmc.R')



### HIGH LEVEL WRAPPER FUNCTION
#####################################
wrapper_object <- proto(expr={})

# Object name, expected child objects & build function
wrapper_object$name <- 'maat'

# model child object
# - a proto object that is the model in which to run the SA/UQ
wrapper_object$model <- NULL

# function to configure wrapper and build model object
wrapper_object$build <- function(., ... ) {

  # build wrapper
  # - a set of if/else statements to initialise the wrapper functions with specific functions
  # - could be built around get paste type setup based on runtype argument
  .$generate_ensemble_pars <- get(paste0('generate_ensemble_pars_',.$wpars$runtype))
  .$init_output_matrix     <- get(paste0('init_output_matrix_',.$wpars$runtype))
  .$write_output           <- get(paste0('write_output_',.$wpars$runtype))
  .$output                 <- get(paste0('output_',.$wpars$runtype))
  .$run0                   <- get(paste0('run0_',.$wpars$runtype))
  .$run1                   <- get(paste0('run1_',.$wpars$runtype))
  .$run2                   <- get(paste0('run2_',.$wpars$runtype))
  .$run3                   <- get(paste0('run3_',.$wpars$runtype))
  .$run4                   <- get(paste0('run4_',.$wpars$runtype))
  .$run5                   <- get(paste0('run5_',.$wpars$runtype))
  .$run6                   <- get(paste0('run6_',.$wpars$runtype))
  .$run7                   <- get(paste0('run7_',.$wpars$runtype))
  .$run8                   <- get(paste0('run8_',.$wpars$runtype))

  # MCMC specific functions
  if(grepl('mcmc',.$wpars$runtype)) {
    .$proposal_generate     <- get(paste0('proposal_generate_',.$wpars$runtype))
    .$proposal_accept       <- get(paste0('proposal_accept_',.$wpars$runtype))
    .$proposal_lklihood     <- get(paste0('f_proposal_lklihood_',.$wpars$mcmc$lklihood))
    .$init_mcmc             <- get(paste0('init_',.$wpars$runtype))
    .$mcmc_outlier          <- get(paste0('mcmc_outlier_', .$wpars$mcmc$outlier))
    .$mcmc_converge         <- get(paste0('mcmc_converge_', .$wpars$mcmc$converge))
    .$mcmc_bdry_handling    <- get(paste0('mcmc_bdry_handling_', .$wpars$mcmc$bdry_handling))
    .$mcmc_prior            <- get(paste0('mcmc_prior_', .$wpars$mcmc$prior))
    .$boundary_handling_set <- boundary_handling_set
  }

  # build model
  setwd(paste('system_models', .$wpars$mod_obj, sep='/' ))
  mod_obj <- paste(.$wpars$mod_obj, 'object', sep='_' )
  source(paste0(mod_obj, '.R' ))
  .$model <- as.proto(get(mod_obj)$as.list(), parent=. )
  .$model$build(...)
  rm(mod_obj)
  setwd('../..')
}

# clean function to reset object
wrapper_object$clean <- function(.) {
  for(d in 1:length(.$dataf)) if(!is.null(.$dataf[[d]])) .$dataf[[d]] <- NA
  gc()
}



# Main run function
###########################################################################
wrapper_object$run   <- function(.,verbose=T) {

  # Initialise
  if(!.$wpars$unit_testing) {
    # if(.$wpars$runtype=='SAprocess_ye')  .$wpars$eval_strings <- T
    if(.$wpars$runtype=='SAprocess_ye' | .$wpars$mcmc)  .$wpars$eval_strings <- T
    .$init()
  } else {
    .$wpars$of_dir       <- '~/tmp'
    .$wpars$of_type      <- 'csv'
    .$wpars$of_name_stem <- 'unit_test'
    hd <- getwd()
  }

  # Initialisation checks
  # need to add a check for equal par vector lengths if this is a UQ run and not eval_strings
  # for Ye et al SA method
  # - due to different parameter sample numbers in process A and B loops,
  # - parameters samples must be generated from code snippets as strings
  # if(.$wpars$eval_string&is.null(.$dynamic$pars_eval)) {
  if(.$wpars$eval_strings & is.null(.$dynamic$pars_eval)) {
    stop(paste('wrapper: eval_strings = T but dynamic$pars_eval not set. \n
          vars$pars_eval:,\n',.$dynamic$pars_eval,'\n
          NOTE: Ye method SA must draw parameter samples during runtime \n
          from code snippets written as strings in dynamic$pars_eval'))
  }

  # initialise model with static variables
  if(!is.null(.$static$pars))   .$model$configure(vlist='pars',   df=.$static$pars   )
  if(!is.null(.$static$env))    .$model$configure(vlist='env',    df=.$static$env    )
  if(!is.null(.$static$fnames)) .$model$configure(vlist='fnames', df=.$static$fnames )

  # create matrices of runtime variables
  .$generate_ensemble()

  # print summary of maat setup
  .$print_data()
  .$print_data(otype='run')

  # store model output template (currently must be a vector)
  # - code will fail when output is a vector of variable length depending on parameter values
  .$dataf$mout <- .$model$output()

  # run model ensemble
  # - call initial run function in the run function cascade
  .$run0()

  # reset wd if unit testing
  if(.$wpars$unit_testing) setwd(hd)
}


# run cascade
###########################################################################
wrapper_object$run0 <- function(.,i) {}
wrapper_object$run1 <- function(.,i) {}
wrapper_object$run2 <- function(.,j) {}
wrapper_object$run3 <- function(.,k) {}
wrapper_object$run4 <- function(.,l) {}
wrapper_object$run5 <- function(.,m) {}
wrapper_object$run6 <- function(.,n) {}
wrapper_object$run7 <- function(.,o) {}
wrapper_object$run8 <- function(.,p) {}


# print function for run cascade
wrapper_object$printc <- function(.,r1,r2) {
  print(r1,quote=F)
  print(r2,quote=F)
}


# takes a >=3 D array and stacks it into a 2D matrix
wrapper_object$stack <- function(., a ) {
  apply(a, 2, function(v) v )
}


# initialisation functions
###########################################################################

# this flattens the object|variable hierarchy in the list structure
# allowing single run matrices that contain variables for multiple model objects
# each line of the matrix is passed to the configure function in the model object
wrapper_object$init <- function(.) {

  # setup list names for assignment
  type   <- c('static', 'dynamic')
  vlists <- c('pars', 'env', 'fnames' )

  # assign standard input lists to wrapper data structure
  for( t in type ) {
    for( vl in vlists ) {
      # input variables
      .[[t]][[vl]] <-
        if(t == 'static') unlist(.[[paste0('init_',t)]][[vl]])
        else if(t == 'dynamic' & !is.null(unlist(.[[paste0('init_',t)]][[vl]])) )
          lapply(rapply(.[[paste0('init_',t)]][[vl]], enquote, how="unlist" ), eval )
    }
  }

  if(.$wpars$UQ|.$wpars$mcmc) .$init_uq()
}


# as above for pars code snippets (pars_eval input) and assigment of parameters to a process (pars_proc input)
wrapper_object$init_uq <- function(.) {

  if(is.null(unlist(.$init_dynamic$pars))&!is.null(unlist(.$init_dynamic$pars_eval))) .$wpars$eval_strings <- T
  t <- 'dynamic'

  if(.$wpars$eval_strings) {
    vl   <- 'pars_eval'
    if(!is.null(unlist(.[[paste0('init_',t)]][[vl]])))
      .[[t]][[vl]] <- lapply(rapply(.[[paste0('init_',t)]][[vl]], enquote, how="unlist" ), eval )
  }

  if(.$wpars$runtype=='SAprocess_ye') {
    vl   <- 'pars_proc'
    if(!is.null(unlist(.[[paste0('init_',t)]][[vl]])))
      .[[t]][[vl]] <- lapply(rapply(.[[paste0('init_',t)]][[vl]], enquote, how="unlist" ), eval )
    for( vn in names(.[[t]][[vl]]) ) if( !any(vn==names(.[[t]][['pars_eval']])) )
      stop(paste('\n Variable:', vn, 'in pars_proc, not found in: pars_eval list.',
                 '\n The proc_pars input list must contain exactly the same parameter names as the pars_eval input list.',
                 '\n The proc_pars is required to assign a parameter to a process as part of a process sensitivity analysis.'))
  }
}



# Wrapper object data struture
###########################################################################

# initialisation lists
wrapper_object$init_static  <- NULL
wrapper_object$init_dynamic <- NULL


# static variables
# each element in the below list is a character or numeric vector to overwrite default initialisation values
wrapper_object$static <- list(
  fnames = NULL,
  pars   = NULL,
  env    = NULL
)


# dynamic variables
# all elements expected to be of class 'list'
# each list in the 'dynamic' list comprise vectors of the values for each variable,
# each element of the list is labelled by the variable name prefixed by the name of the model object that the variable belongs to
# each of these lists is expanded, often factorially by expand.grid, and placed into the below list of dataframes
wrapper_object$dynamic <- list(
  fnames    = NULL,
  fnamesB   = NULL,
  pars      = NULL,
  # list with same elements and names as pars but giving the fnames list name i.e. the process name to which each parameter belongs
  pars_proc = NULL,
  # list with same elements and names as pars but each element is a code snippet as a string that once evaluated gives a vector or parameter values
  # allows different types of distributions to be specified for each parameter
  # this must be used for the Ye SA method
  pars_eval = NULL,
  parsB     = NULL,
  env       = NULL
)


# input/output matrices and dataframes
# with an associated length for input matrices
wrapper_object$dataf  <- list(
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
  mcmc_input   = NULL,    # list of output matrices and arrays from a previous MCMC run
  mout         = NULL,    # example model output vector, for setting up vapply functions
  out          = NULL,    # output matrix
  out_saltelli = NULL,    # saltelli output list
  # observation matrices /dataframes
  obs          = NULL,    # a dataframe of observations against which to valiadate/ calculate likelihood of model
  obsse        = NULL     # a dataframe of observation errors for the obs data, must exactly match the above dataframe
)


# parameters specific to the wrapper object
wrapper_object$wpars <- list(
  multic          = F,             # multicore the simulation
  procs           = 6,             # number of processors to use if multic = T
  cverbose        = F,             # write configuration output during runtime
  UQ              = F,             # run a UQ analysis
  runtype         = 'none',        # ensemble type - 'factorial', 'SApar_saltelli', and 'SAprocess_ye' available so far
  of_dir          = '~/tmp',       # output directory
  of_type         = 'csv',         # output file type - 'csv' or 'rds'
  of_name         = '',            # output file name - excluding file extension
  of_name_history = '',            # history output file name - excluding file extension
  of_name_stem    = 'MAAT_output', # output file name stem - all output file in an ensemble will begin with this
  n               = numeric(1),    # parameter sample number
  parsinit_read   = F,             # parameter samples have been read from a file
  nmult           = 1,             # parameter sample number multiplier for saltelli method
  eval_strings    = F,             # switch telling wrapper that vars$pars are to be evaluated from code string snippets in vars$pars_eval
  sobol_init      = T,             # initialise sobol sequence or not when calling rsobol. This should not be modified by the user.
  unit_testing    = F,

  mcmc = list(
    run_type      = 'dream',
    lklihood      = 'ssquared',
    outlier       = 'iqr',
    mcmc_converge = 'Gelman_Rubin',
    bdry_handling = 'bound',
    init_prior    = 'uniform',
    chains        = 7,
    maxiter       = 1000,
    start_iter    = 2,
    thin          = 0.1,
    thin_obs      = 1,
    homosced      = F,
    chain_delta   = 3,
    c_rand        = 0.01,
    c_ergod       = 1e-12,
    p_gamma       = 0.2,
    n_CR          = 3,
    adapt_pCR     = T,
    CR_burnin     = 1e4,
    check_ss      = numeric(1),
    check_iter    = 10
  )
)


# MCMC specific data, size depends on MCMC set up
wrapper_object$mcmc <- list(
  outlier_detected = F,
  j_start_burnin   = 1,
  j_burnin50       = numeric(1),
  d                = numeric(1),
  CR               = numeric(1),
  p_CR             = numeric(1),
  R                = matrix(),
  current_state    = matrix(),
  p_state          = numeric(1),
  sd_state         = numeric(1),
  jump             = matrix(),
  draw             = matrix(),
  lambda           = matrix(),
  boundary_min     = numeric(1),
  boundary_max     = numeric(1),
  del              = numeric(1),
  L                = numeric(1),
  t                = numeric(1),
  m                = numeric(1),
  CR_burnin        = T,
  d_star           = numeric(1)
)



# Output processing functions
###########################################################################

# function to combine factorial ensemble with met data
# - for each ensemble member all columns of met matirx are run
# - this is called from an lapply to expand each each ensemble member values of fnames, pars, and env with every column of the met matrix
wrapper_object$combine <- function(.,i,df) suppressWarnings(data.frame(t(.$dataf$met),df[i,]))


# function to write ensemble output data to file
wrapper_object$write_to_file <- function(., df=.$output(), app=F ) {

  setwd(.$wpars$of_dir)
  if(.$wpars$of_type=='csv')      write.table(format(df,width=12), paste(.$wpars$of_name,'.csv',sep=''),
                                              quote=F, row.names=F, col.names=!app, sep=',', append=app )
  else if(.$wpars$of_type=='rds') saveRDS(df, paste(.$wpars$of_name,'RDS',sep='.') )
  else print(paste('No methods for output file format:',.$wpars$of_type))
}



# Print functions
###########################################################################
wrapper_object$print_data <- function(.,otype='data') {

  ens_n <-
    if(.$wpars$UQ) {
      if(.$wpars$runtype=='SAprocess_ye') .$dataf$lf * prod(unlist(lapply(.$dynamic$fnames,length))) * .$wpars$n^2 * .$dataf$le
      else                                .$dataf$lf * .$wpars$n * (2+dim(.$dataf$pars)[2]) * .$dataf$le
    } else                                .$dataf$lf * .$dataf$lp *.$dataf$le

  if(otype=='data') {

    print('',quote=F)
    print('',quote=F)
    print('',quote=F)
    print('',quote=F)
    print("MAAT :: summary of data",quote=F)
    print('',quote=F)
    print('',quote=F)
    print("fnames ::",quote=F)
    if(!is.null(.$dataf$fnames)) print(summary(t(.$dataf$fnames)), quote=F )
    else                         print(NULL, quote=F )

    print('',quote=F)
    print("pars ::",quote=F)
    if(!.$wpars$runtype=='SAprocess_ye')
      if(!is.null(.$dataf$pars)) print(summary(t(.$dataf$pars)), quote=F )
      else                       print(NULL, quote=F )
    else {
      print(.$dynamic$pars_proc,quote=F)
      print(paste('sample n:',.$wpars$n),quote=F)
    }

    print('',quote=F)
    print("env ::",quote=F)
    if(!is.null(.$dataf$env)) print(summary(t(.$dataf$env)), quote=F )
    else                      print(NULL, quote=F )

    print('',quote=F)
    print("met data ::",quote=F)
    if(!is.null(.$dataf$met)) print(summary(t(.$dataf$met)), quote=F )
    else                      print(NULL, quote=F )
    print('',quote=F)

  } else if(otype=='run') {

    print('',quote=F)
    print('',quote=F)
    print('',quote=F)
    print(paste("MAAT :: run model",Sys.time()), quote=F )
    print('',quote=F)
    print(paste(.$wpars$runtype,' ensemble'), quote=F )
    print(paste('ensemble number ::',ens_n), quote=F )
    if(!is.null(.$dataf$met)) {
      print(paste('timesteps in met data ::',.$dataf$lm), quote=F )
      print(paste('total number of model calls ::',ens_n*.$dataf$lm), quote=F )
    }
    print('',quote=F)

    if(.$wpars$multic) print(paste('parallel processing over ::',.$wpars$procs,'cores.'), quote=F )
    else               print(paste('serial processing.'), quote=F )
    print('',quote=F)
  }
}

wrapper_object$print_output <- function(.) {
  print("output ::",quote=F)
  print(paste('length ::', length(.$dataf$out[,1])), quote=F)
  print(head(.$dataf$out), quote=F)
  print('', quote=F)
  print('', quote=F)
  print(Sys.time(), quote=F)
  print('', quote=F)
}

wrapper_object$print_saltelli <- function(.) {
  print('Saltelli matrix AB completed', quote=F)
  print('', quote=F)
  print('', quote=F)
  print('run Saltelli array ABi', quote=F)
  print('', quote=F)
}



######################################################################################
# Unit testing functions

# simple test, run a single model instance with or without met data
wrapper_object$.test_simple <- function(., gen_metd=F, mc=F, pr=4, oconf=F ) {

  # source directory
  .$wpars$UQ      <- F           # run a fully factorial ensemble
  .$wpars$runtype <- 'factorial'
  .$wpars$mod_obj <- 'leaf'
  .$build()

  # verbose params
  .$model$pars$verbose  <- F
  .$model$pars$cverbose <- oconf

  # define parameters for the wrapper
  .$wpars$multic       <- mc  # multicore the ensemble
  .$wpars$procs        <- pr  # number of cores to use if above is true
  .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening

  # Define meteorological and environment dataset
  .$model$env$par     <- 1000
  .$model$env$ca_conc <- 400
  .$model$env$vpd     <- 1
  .$model$env$temp    <- 20
  metdata <- t(as.matrix(expand.grid(list(leaf.par = seq(800,1000,100),leaf.ca_conc = 400))))
  if(gen_metd) .$dataf$met <- metdata

  # Define the static parameters and model functions
  .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')
  .$dynamic$fnames <- list(
    leaf.etrans = c('f_etrans_farquhar1980')
  )

  .$dynamic$env <- list(
    leaf.vpd  = 1
  )

  .$dynamic$pars <- list(
    leaf.avn_25 = 10
  )

  # Run wrapper & model
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

  # process & record output
  .$output() # this should be runtype specific, set up in build
}


# general factorial test, with or without metdata
wrapper_object$.test <- function(.,gen_metd=T,mc=T,pr=4,oconf=F) {

  library(lattice)

  # build wrapper and the model object
  .$wpars$UQ      <- F           # run a fully factorial ensemble
  .$wpars$runtype <- 'factorial'
  .$wpars$mod_obj <- 'leaf'
  .$build()

  # define parameters for the wrapper
  .$wpars$multic       <- mc  # multicore the ensemble
  .$wpars$procs        <- pr  # number of cores to use if above is true
  .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble
  .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening

  ### Define meteorological and environment dataset
  ###############################
  # can load a met dataset here
  # below a trivial met dataset is created to be used as an example
  metdata <- t(as.matrix(expand.grid(list(leaf.par = seq(0,1000,100),leaf.ca_conc = 400))))
  #metdata <- as.matrix(expand.grid(list(leaf.par = seq(800,1000,100),leaf.ca_conc = 400)))
  if(gen_metd) .$dataf$met <- metdata
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
    leaf.etrans = c('f_etrans_farquhar1980','f_etrans_collatz1991'),
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
  st <- system.time(.$run())
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
wrapper_object$.test_init <- function(.,
                      metd=list(leaf.par = seq(800,1000,100), leaf.ca_conc = 400 ),
                      sfnames=list(fnames=list(leaf=list(vcmax='f_vcmax_lin'))),
                      spars=list(pars=list(leaf=list(Ha=list(vcmax=7e4, jmax=4e4 )))),
                      senv=list(env=list(leaf=list(par=1000, ca_conc=400 ))),
                      dfnames=NULL,
                      dpars=NULL,
                      denv=NULL
                      ) {

  # Define the static model functions, parameters, and environment
  .$init_static <- c(sfnames, spars, senv )
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
wrapper_object$.test_con <- function(., mc=T, pr=4, oconf=F,
                                     metd=list(leaf.par = seq(800,1000,100),leaf.ca_conc = 400),
                                     sfnames=list(leaf.vcmax='f_vcmax_lin'),
                                     spars=NULL,
                                     senv=list(par=1000, ca_conc=400 ),
                                     dfnames=list(
                                       leaf.etrans = c('f_etrans_farquhar1980','f_etrans_collatz1991'),
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

  # build wrapper and the model object
  .$wpars$UQ      <- F           # run a fully factorial ensemble
  .$wpars$runtype <- 'factorial'
  .$wpars$mod_obj <- 'leaf'
  .$build()

  # verbose parameters
  .$model$cpars$verbose  <- F
  .$model$cpars$cverbose <- oconf

  # define parameters for the wrapper
  .$wpars$multic       <- mc  # multicore the ensemble
  .$wpars$procs        <- pr  # number of cores to use if above is true
  .$wpars$UQ           <- F   # run a UQ style ensemble, or if false a fully factorial ensemble
  .$wpars$unit_testing <- T   # tell the wrapper unit testing is happening

  # Define meteorological dataset
  if(!is.null(metd))    .$dataf$met  <- t(as.matrix(expand.grid(metd)))

  # Define the static model functions, parameters, and environment
  if(!is.null(sfnames)) .$static$fnames <- sfnames
  if(!is.null(spars))   .$static$pars   <- spars
  if(!is.null(senv))    .$static$env    <- senv

  # Define the dynamic model functions, parameters, and environment
  if(!is.null(dfnames)) .$dynamic$fnames <- dfnames
  if(!is.null(dpars))   .$dynamic$pars   <- dpars
  if(!is.null(denv))    .$dynamic$env    <- denv

  # Run wrapper & model
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

  # process & record output
  .$output()
}


# general factorial test with canopy object, with or without metdata
wrapper_object$.test_can <- function(., metd=T, mc=T, pr=4, verbose=F ) {

  # build wrapper and the model object
  .$wpars$UQ      <- F           # run a fully factorial ensemble
  .$wpars$runtype <- 'factorial'
  .$wpars$mod_obj <- 'canopy'
  .$build()

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
  metdata <- t(as.matrix(expand.grid(list(canopy.par_dir = 500,canopy.ca_conc = seq(10,1200,50)))))
  if(metd) .$dataf$met <- metdata

  # Define the parameters and model functions that are to be varied
  .$dynamic$fnames <- list(
    canopy.can_scale_light = c('f_canlight_beerslaw_wrong','f_canlight_beerslaw'),
    leaf.etrans            = c('f_etrans_farquhar1980','f_etrans_collatz1991'),
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
  st <- system.time(.$run())
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
wrapper_object$.test_mimic <- function(., mod_mimic='clm45_non_Tacclimation', mod_obj='leaf', metd=F, mc=F, pr=4, oconf=F ) {

  # build wrapper and the model object
  .$wpars$UQ      <- F
  .$wpars$runtype <- 'factorial'
  .$wpars$mod_obj <- mod_obj
  .$build()

  # define control parameters
  .$model$cpars$verbose  <- F
  .$model$cpars$cverbose <- oconf
  .$wpars$multic         <- mc  # multicore the ensemble
  .$wpars$procs          <- pr  # number of cores to use if above is true
  .$wpars$unit_testing   <- T   # tell the wrapper unit testing is happening

  # Define the static parameters and model functions
  init_static = list( env = list( leaf = list()))
  init_static$env$leaf$ca_conc <- 400
  init_static$env$leaf$par     <- 2000
  init_static$env$leaf$vpd     <- 50
  init_static$env$leaf$temp    <- 25
  .$init_static <- init_static

  # Define the parameters and model functions that are to be varied
  # if this is a UQ analysis, the "pars" list must contain parameter vectors that are of equal length,
  # if not a UQ analysis the parameter vectors in the "pars" list can be of different lengths
  .$init_dynamic <- NULL

  # run wrapper init to test mimic
  .$init()

  # Define meteorological and environment dataset
  # can load a met dataset here
  # below a trivial met dataset is created to be used as an example
  metdata <- t(as.matrix(expand.grid(list(leaf.par = seq(0,1000,100), leaf.ca_conc = 400))  ))
  metdata <- t(as.matrix(expand.grid(list(leaf.par = 2000, leaf.ca_conc = seq(50,1500,50))) ))
  if(metd) .$dataf$met <- metdata

  # Run model
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

  # # process & record output
  .$output()
}


# test function for Ye method Sobol process sensitivity analysis
wrapper_object$.test_ye <- function(.,metd=F,mc=T,pr=4,oconf=F,n=3) {

  # build wrapper and the model object
  .$wpars$UQ      <- T
  .$wpars$runtype <- 'SAprocess_ye'
  .$wpars$mod_obj <- 'leaf'
  .$build()

  # define control parameters
  .$model$pars$verbose  <- F
  .$model$pars$cverbose <- oconf
  .$wpars$multic        <- mc   # multicore the ensemble
  .$wpars$procs         <- pr   # number of cores to use if above is true
  .$wpars$unit_testing  <- T    # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions)

  ### Define static variables
  ###############################
  .$model$env$par <- 1000

  ### Define the parameters and model functions that are to be varied
  ###############################
  # "pars" lists must contain parameter vectors that are of equal length,

  # add the SA/UQ variables to the maat wrapper object
  .$wpars$n            <- n    # number of parameter samples in each loop
  .$wpars$coef_var     <- 0.1  # coefficient of variation for prior parameter distribution
  .$wpars$eval_strings <- T    # use evaluation strings to set parameter values


  # - the wrapper object takes care of combining these lists into the full ensemble
  .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')
  .$dynamic$fnames <- list(
    leaf.Alim   = c('f_Alim_farquhar1980','f_Alim_collatz1991'),
    leaf.etrans = c('f_etrans_farquharwong1984','f_etrans_collatz1991','f_etrans_harley1992')
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
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

}

# test function for Saltelli method Sobol parametric sensitivity analysis
wrapper_object$.test_saltelli <- function(., metd=F, mc=T, pr=4, oconf=F, n=3, eval_strings=T ) {

  # build wrapper and the model object
  .$wpars$UQ      <- T
  .$wpars$runtype <- 'SApar_saltelli'
  .$wpars$mod_obj <- 'leaf'
  .$build()

  # define control parameters
  .$model$pars$verbose  <- F
  .$model$pars$cverbose <- oconf
  .$wpars$multic       <- mc           # multicore the ensemble
  .$wpars$procs        <- pr           # number of cores to use if above is true
  .$wpars$unit_testing <- T            # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions)

  ### Define static variables
  .$static$env    <- list(leaf.par=1000)
  .$static$fnames <- list(leaf.vcmax='f_vcmax_lin')

  ### Define the parameters and model functions that are to be varied
  # "pars" lists must contain parameter vectors that are of equal length,

  # add the SA/UQ variables to the maat wrapper object
  # - the wrapper object takes care of combining these lists into the full ensemble
  .$wpars$n            <- n            # number of parameter samples in each loop
  .$wpars$eval_strings <- eval_strings # parameters are passed as strings to be evaluated to allow for different sample numbers
  .$wpars$coef_var     <- 0.1

  .$dynamic$fnames <- list(
    leaf.Alim   = c('f_Alim_farquhar1980','f_Alim_collatz1991'),
    leaf.etrans = c('f_etrans_farquharwong1984','f_etrans_collatz1991','f_etrans_harley1992')
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
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

  # process & record output
  list(AB=.$dataf$out, ABi=.$dataf$out_saltelli)
}

# test function for MCMC parameter estimation using mixture model with tri-modal distribution
wrapper_object$.test_mcmc_mixture <- function(., mc=F, pr=4, mcmc_type='dream',
                                              mcmc_chains=8, mcmc_maxiter=100,
                                              mu_vector=c(-8,0,8),
                                              sd_vector=c(1,1,1),
                                              height_vector=c(0.2,0.5,0.3),
                                              mixture_scale=1e12,
                                              verbose=F, cverbose=F, diag=F
                                              ) {


  library(lattice)
  # build wrapper and the model object
  .$wpars$UQ      <- T
  .$wpars$runtype <- paste0('mcmc_',mcmc_type)
  .$wpars$mod_obj <- 'mcmc_test'

  # define control parameters
  .$wpars$unit_testing  <- T            # tell the wrapper unit testing is happening - bypasses the model init function (need to write a separate unite test to test just the init functions)
  .$wpars$verbose       <- verbose
  .$wpars$cverbose      <- cverbose
  .$wpars$multic        <- mc           # multicore the ensemble
  .$wpars$procs         <- pr           # number of cores to use if above is true
  .$wpars$UQ            <- T            # run a UQ/SA style ensemble
  .$wpars$UQtype        <- 'mcmc'       # MCMC ensemble
  .$wpars$mcmc$type     <- mcmc_type    # MCMC type, 'demc' or 'dream'
  .$wpars$mcmc$chains   <- mcmc_chains  # MCMC number of chains
  .$wpars$mcmc$maxiter  <- mcmc_maxiter # MCMC max number of steps / iterations on each chain
  .$wpars$mcmc$lklihood <- 'log'        # MCMC likelihood function
  .$build(mod_out='mixture', switches=c(diag,verbose,cverbose) )

  # set model system function
  .$model$fnames$sys    <- 'f_sys_mixture'

  # Define static variables
  .$static$fnames <- list(mcmc_test.sys='f_sys_mixture')

  # set problem specific parameters
  .$model$pars$mu1      <- mu_vector[1]
  .$model$pars$mu2      <- mu_vector[2]
  .$model$pars$mu3      <- mu_vector[3]
  .$model$pars$sd1      <- sd_vector[1]
  .$model$pars$sd2      <- sd_vector[2]
  .$model$pars$sd3      <- sd_vector[3]
  .$model$pars$height1  <- height_vector[1]
  .$model$pars$height2  <- height_vector[2]
  .$model$pars$height3  <- height_vector[3]
  .$model$pars$mixture_scale  <- mixture_scale

  # define priors
  .$dynamic$pars_eval <- list(
    mcmc_test.proposal1  = 'runif(n,-20,20)',
    mcmc_test.proposal2  = 'runif(n,-20,20)',
    mcmc_test.proposal3  = 'runif(n,-20,20)',
    mcmc_test.proposal4  = 'runif(n,-20,20)'
  )

  # define ofname
  .$ofname <- 'mcmc_mixture_test'

  # Run MCMC
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

  # process & record output
  mcmc_pars_hist <-
    hist(.$dataf$pars_array, breaks=200,
         col='darkmagenta', border='darkmagenta',
         xlab='Mixture Model Parameters', main='Posterior (Target) Parameter Distributions for Mixture Model')

  dims <- dim(.$dataf$pars_lklihood)
  df1  <- data.frame(lklihood=as.vector(t(.$dataf$pars_lklihood)), chain=rep(1:dims[1],each=dims[2]) )
  lklihood_plot <-
    xyplot(lklihood ~ rep(1:dims[2],dims[1]), df1, groups=chain, auto.key=T, type='l' )
  lklihood_plot2 <-
    xyplot(lklihood ~ rep(1:dims[2],dims[1])|chain, df1, auto.key=T, type='l' )

  # output
  list(pars_array=.$dataf$pars_array, pars_lklihood=.$dataf$pars_lklihood,
       hist=mcmc_pars_hist, lklihood_plot=lklihood_plot, lklihood_plot2=lklihood_plot2 )
}


# test function for MCMC parameter estimation in a linear regression
wrapper_object$.test_mcmc_linreg <- function(., mc=F, mcmc_chains=7, pr=mcmc_chains,
                                             mcmc_type='dream', mcmc_lklihood='ssquared',
                                             mcmc_homosced=T, mcmc_maxiter=3,
                                             x=1:10, a_mu=-5, b_mu=15, a_sd=1, b_sd=1, standard_err=0.5,
                                             mcmc_test.a  = 'runif(n,-30,30)',
                                             mcmc_test.b  = 'runif(n,-30,30)',
                                             verbose=F, cverbose=F, diag=F
                                             ) {

  ### currently does not work with multicoring,
  ### probably due to assignment to . datastructure during the forked processes
  ### perhaps setting parent in build function would avoid problem, should be shared memory, but perhaps not
  # redefine mc if mc=T
  #if(mc) {mc <- F; print('mc specified as T but does not work in unit testing with current MAAT config, redefining mc as F') }

  library(lattice)
  # build wrapper and the model object
  .$wpars$UQ      <- T
  .$wpars$runtype <- paste0('mcmc_',mcmc_type)
  .$wpars$mod_obj <- 'mcmc_test'

  # define control parameters
  .$wpars$unit_testing  <- T
  .$wpars$verbose       <- verbose
  .$wpars$cverbose      <- cverbose
  .$wpars$mcmc$lklihood <- mcmc_lklihood # MCMC likelihood function
  .$wpars$multic        <- mc            # multicore the ensemble
  .$wpars$procs         <- pr            # number of cores to use if above is true
  .$wpars$UQ            <- T             # run a UQ/SA style ensemble
  .$wpars$UQtype        <- 'mcmc'        # MCMC ensemble
  .$wpars$mcmc$type     <- mcmc_type     # MCMC type, 'demc' or 'dream'
  .$wpars$mcmc$chains   <- mcmc_chains   # MCMC number of chains
  .$wpars$mcmc$homosced <- mcmc_homosced # MCMC homoscedastic error
  .$wpars$mcmc$maxiter  <- mcmc_maxiter  # MCMC max number of steps / iterations on each chain
  .$build(mod_out='regression', switches=c(diag,verbose,cverbose) )

  # Define static variables
  .$static$fnames <- list(
    mcmc_test.sys      = 'f_sys_regression',
    mcmc_test.reg_func = 'f_reg_func_linear'
  )

  # set problem specific parameters
  .$model$pars$obs_error <- standard_err
  .$model$pars$syn_a_mu  <- a_mu
  .$model$pars$syn_b_mu  <- b_mu
  .$model$pars$syn_a_sd  <- a_sd
  .$model$pars$syn_b_sd  <- b_sd

  # met data
  .$dataf$met            <- t(matrix(x, length(x), 1 ))
  rownames(.$dataf$met)  <- 'mcmc_test.linreg_x'

  # define priors
  .$dynamic$pars_eval <- list(
    mcmc_test.a  = mcmc_test.a,
    mcmc_test.b  = mcmc_test.b
  )

  # define ofname
  .$ofname <- 'lin_reg_test'

  # Run MCMC
  st <- system.time(.$run())
  print('',quote=F)
  print('Run time:',quote=F)
  print(st)
  print('',quote=F)

  # process & record output
  mcmc_pars_hist <-
    hist(.$dataf$pars_array, breaks=200,
         col='darkmagenta', border='darkmagenta',
         xlab='Regression Model Parameters', main='Posterior (Target) Parameter Distributions for Regression Model')

  dims <- dim(.$dataf$pars_lklihood)
  df1  <- data.frame(lklihood=as.vector(t(.$dataf$pars_lklihood)), chain=rep(1:dims[1],each=dims[2]) )
  lklihood_plot <-
    xyplot(lklihood ~ rep(1:dims[2],dims[1]), df1, groups=chain, auto.key=T, type='l' )
  lklihood_plot2 <-
    xyplot(lklihood ~ rep(1:dims[2],dims[1])|chain, df1, auto.key=T, type='l' )

  # output
  list(pars_array=.$dataf$pars_array, pars_lklihood=.$dataf$pars_lklihood,
       hist=mcmc_pars_hist, lklihood_plot=lklihood_plot, lklihood_plot2=lklihood_plot2 )
}



### END ###
