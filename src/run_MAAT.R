################################
#
# MAAT Model - run script
#
# AWalker (walkerap@ornl.gov)
# December 2015
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script is called to run MAAT

###
# This script initialises MAAT and runs it in the following steps:
# - set default arguments
# - parse command line arguments
# - set arguments that depend on other arguments
# - load MAAT objects from source
# - load init scripts
# - Configure and initialise MAAT
# - Run MAAT

###################################################################

# any one of the below objects to line 150 or so can be specified as a single string command line argument to this script
# the sub-arguments in the string (separated by a space) are interpreted individually as R code.

#       Rscript run_MAAT.R "object1<-value1 object2<-value2"
#  e.g. Rscript run_MAAT.R "dir<-'/home/alp/' multic<-T"

##################################
# command line options and defaults

# model object to use, any directory name in the src/system_models directory
mod_obj <- NULL

# directory paths - set these on the command line or modify these here
# source directory (full path)
# - must be added (here or as a commandline option to this script) before model will run
# - can be modified to target source copied to a static directory rather that a reposiotory with version control
srcdir  <- NULL
# project directory (full path)
pdir    <- NULL
# meteorological data directory (full path)
mdir    <- NULL
# evalutaion data directory (full path)
edir    <- NULL
# output data directory (full path)
odir    <- NULL

# wrapper object options
# multicore the ensemble
multic  <- F

# number of cores to use if above is true
procs   <- 4

# run an emsemble that combines variables in factorial
# - if set to TRUE this will over-ride a UQ analysis
factorial  <- T

# run an SA/UQ style ensemble, or if -uq- is false a fully factorial ensemble
uq         <- F

# types of SA/UQ run
# process SA
procSA     <- T
# Saltelli Sobol SA
salt       <- F
# MCMC parameter estimation
mcmc       <- F
# type of MCMC
mcmc_type  <- 'dream'

# initialise pars with output from an MCMC run
# - used either to 'restart' an MCMC run or run a predictive 'ensemble' by sampling posterior distributions
# - if 'ensemble' factorial must be TRUE
parsinit_mcmc <- NULL
# directory with MCMC data - relative to pdir
mcmcdir       <- NULL
# output filename for MCMC data
mcmcout       <- NULL
# number of samples from MCMC posterior
parsinit_n    <- NULL


# run options
# meteorological data file name
metdata       <- NULL

# evaluation data file name
evaldata      <- NULL 

# load standard error from evaluation data file (T/F)
evalse        <- F

# pass code snippets as strings to generate parameter samples, see init_MAAT.R and wrapper for use
eval_strings  <- F

# initialise in the configuration of the below specified model
# options are: clm40,
mod_mimic     <- NULL

# static and dynamic initialisation files are XMLs, if false init file is the R script named below
xml           <- F

# initialisation data file name if not an XML
init          <- 'init_MAAT'

# run i.d. - used as suffix/prefix for in/out files
runid         <- NULL

# basic output file name
of_main       <- 'out'

# model output switch
mod_out       <- 'run'

# output file format.  supported: rds, csv (default)
of_format     <- 'csv'

# verbose - ouput various things during runtime for diagnostics
verbose       <- F
cverbose      <- F
diag          <- F

# parameters for SA run
# ensemble number for an SA/UQ style ensemble, not used if -uq- is false
psa_n      <- 10
# coefficient of variation if used in parameter sampling in process sensitivity analysis
coef_var   <- 0.1
# multiplier on process ensemble n for Saltelli ensemble n
salt_nmult <- 100

# MCMC parameters
# likelihood function (options: log, ssquared, ssquared_se, ...)
mcmc_lklihood       <- 'ssquared'
# outlier handling (options: none, iqr)
mcmc_outlier        <- 'iqr'
# MCMC convergence testing (options: none, Gelman_Rubin)
mcmc_converge       <- 'Gelman_Rubin'
# MCMC option for parameter treatment in bounded search spaces (options: none, bound, reflect, fold)
boundary_handling   <- 'fold'
# number of chains to run (minumum = 2 * mcmc_chain_delta + 1)
mcmc_chains         <- 7
# number of samples from prior distribution to initialise past_states, recommended 10 * number of estimated parameters  
mcmc_prior_n        <- 40
# number of iterations / steps in chain
mcmc_maxiter        <- 1000
# iterations before convergence checking & adapt pCR  starts, not used when a restart 
mcmc_preburnin_iter <- 100 
# thinning for posterior, as a proportion
mcmc_thin           <- 0.1
# thinning for observations, as a proportion
mcmc_thin_obs       <- 1
# thinning for observations, random (T) or evenly spaced (F)
mcmc_thin_obs_random<- F
# option to assume homoscedastic error in measured observations (else, heteroscedastic)
mcmc_homosced       <- F
# DREAM number chain pairs in proposal
# - max number of chain pairs used to calculate the jump for each chain
mcmc_chain_delta    <- 3
# DREAM randomization, scalar noise min/max = 1 +- c_rand 
mcmc_c_rand         <- 0.01
# DREAM ergodicicty, additive noise, c_ergod is sd of normal distribution - should be small relative to width of posterior distribution
mcmc_c_ergod        <- 1e-12
# DREAM probability of unit jump rate (i.e. probability gamma = 1)
mcmc_p_gamma        <- 0.2
# DREAM-ZS probability of orthogonal 'snooker' update 
mcmc_psnooker       <- 0.1
# DREAM number of crossover values 
mcmc_n_CR           <- 3
# adapt probability of selecting crossover values
mcmc_adapt_pCR      <- T
# if true, number of iteration to adapt crossover selection probabilities, "burnin" for crossover prob adjustement
mcmc_CR_burnin      <- 1e4
# checking for convergence and outlier chains every N iterations
mcmc_check_iter     <- 10
# conv_period * check_iter iterations is the averaging period over which to assess convergence, once converged stop outlier detection (DREAM only) 
mcmc_conv_period    <- 40
# every N iterations append past_states matrix with updated state on all chains (DREAMZS only)
mcmc_iterappend     <- 10



##################################
# parse command line arguments

print('',quote=F); print('',quote=F)
print('Read command line arguments:',quote=F)
print('',quote=F)
print(commandArgs(T),quote=F)
if(length(commandArgs(T))>=1) {
  for( ca in 1:length(commandArgs(T)) ) {
    eval(parse(text=commandArgs(T)[ca]))
  }
}

print('',quote=F); print('',quote=F)
if(is.null(srcdir))  stop('srcdir argument not specified but required, check command line arguments to run_MAAT.R')
print(paste('Source directory:',srcdir) ,quote=F)
if(is.null(pdir))    stop('pdir argument not specified but required, check command line arguments to run_MAAT.R')
print(paste('Project directory:',pdir) ,quote=F)
if(is.null(mod_obj)) stop('mod_obj argument not specified but required, check command line arguments to run_MAAT.R')
print(paste('Model:',mod_obj) ,quote=F)
if(!is.null(runid)) if(runid=='') runid <- NULL
print(paste('Run ID:',runid) ,quote=F)

# set default values if not specified on command line
# - these are set after parsing command line arguments as they depend on other arguments that could be set on the command line
if(is.null(of_main)) of_main <- proj

# create output directory
setwd(pdir)
if(is.null(odir)) {
  date   <- Sys.Date()
  odir1  <- paste(pdir,'results',sep='/')
  odir   <- paste(odir1,date,sep='/')
  if(!file.exists(odir1)) dir.create(odir1)
  if(!file.exists(odir))  dir.create(odir)
}

# create input/output filenames
# initialisation file if not xml
initf       <- if(is.null(runid))     paste(init, 'R', sep='.' ) else paste(init, '_', runid, '.R', sep='' )
# prefix for output files
ofname      <- if(is.null(runid))     of_main                    else paste(of_main, runid, sep='_' )
ofname      <- if(is.null(mod_mimic)) ofname                     else paste(mod_mimic, ofname, sep='_' )

# factorial analysis over-rides UQ analysis
#if(!uq) factorial <- T
if(factorial&(uq|mcmc)) {
 uq   <- F
 mcmc <- F
 print('',quote=F); print('',quote=F)
 print(paste('Both factorial and UQ or MCMC run specified: Factorial ensemble will be run'),quote=F)
}

# select run/ensemble type and output
runtype <-
  if(factorial) 'factorial'      else
  if(procSA)    'SAprocess_ye'   else
  if(salt)      'SApar_saltelli' else
  if(mcmc)      'mcmc'

print('',quote=F); print('',quote=F)
print(paste('Run type selected:',runtype) ,quote=F)

if((uq|mcmc)&of_format!='rds') {
  of_format <- 'rds'
  print('',quote=F); print('',quote=F)
  print(paste('of_format changed to rds due to high output volume with SA/UQ ensembles'),quote=F)
}



##################################
# Clone and assign arguments to the maat wrapper

setwd(srcdir)
source('wrapper_object.R')
maat <- as.proto(wrapper_object$as.list())
rm(wrapper_object)

# define run parameters
maat$wpars$mod_obj       <- mod_obj
maat$wpars$runtype       <- runtype
maat$wpars$UQ            <- uq
maat$wpars$multic        <- multic
maat$wpars$procs         <- procs
maat$wpars$n             <- psa_n
maat$wpars$coef_var      <- coef_var
maat$wpars$nmult         <- salt_nmult
maat$wpars$eval_strings  <- eval_strings
maat$wpars$of_name_stem  <- ofname
maat$wpars$of_type       <- of_format
maat$wpars$of_dir        <- odir
maat$wpars$parsinit_read <- !is.null(parsinit_mcmc)

# define MCMC run parameters (sublist of wpars)
# APW: if statement doesn't work for all MCMC functions as maat builds immediately after this before the restart can be read, fix 
maat$wpars$mcmc$mcmc_type         <- mcmc_type
maat$wpars$mcmc$lklihood          <- mcmc_lklihood
maat$wpars$mcmc$outlier           <- mcmc_outlier
maat$wpars$mcmc$boundary_handling <- boundary_handling
maat$wpars$mcmc$converge          <- mcmc_converge
if((runtype=='mcmc') & is.null(parsinit_mcmc)) { 
  maat$wpars$mcmc$chains        <- mcmc_chains
  maat$wpars$mcmc$prior_n       <- mcmc_prior_n
  maat$wpars$mcmc$thin          <- mcmc_thin
  maat$wpars$mcmc$thin_obs      <- mcmc_thin_obs
  maat$wpars$mcmc$thin_obs_random <- mcmc_thin_obs_random
  maat$wpars$mcmc$homosced      <- mcmc_homosced
  maat$wpars$mcmc$chain_delta   <- mcmc_chain_delta
  maat$wpars$mcmc$c_rand        <- mcmc_c_rand
  maat$wpars$mcmc$c_ergod       <- mcmc_c_ergod
  maat$wpars$mcmc$p_gamma       <- mcmc_p_gamma
  maat$wpars$mcmc$psnooker      <- mcmc_psnooker
  maat$wpars$mcmc$n_CR          <- mcmc_n_CR
  maat$wpars$mcmc$adapt_pCR     <- mcmc_adapt_pCR
  maat$wpars$mcmc$CR_burnin     <- mcmc_CR_burnin
  # APW: unnecessary, fix 
  if((mcmc_maxiter/mcmc_check_iter)<3) mcmc_check_iter <- mcmc_check_iter / 2
  if(mcmc_check_iter<5) {
    mcmc_check_iter <- 4                    #ALJ: prevents numerical error in indexing ??
    if(mcmc_maxiter<12) mcmc_maxiter <- 12
  } 
  maat$wpars$mcmc$check_iter     <- mcmc_check_iter
  maat$wpars$mcmc$conv_period    <- mcmc_conv_period
  maat$wpars$mcmc$iterappend     <- mcmc_iterappend
  maat$wpars$mcmc$maxiter        <- mcmc_maxiter
  maat$wpars$mcmc$preburnin_iter <- mcmc_preburnin_iter
}

# build maat and model objects
maat$build(mod_mimic=mod_mimic, mod_out=mod_out )

# set debugging flags
# APW: these could now be set as an argument to build 
maat$model$cpars$verbose  <- verbose
maat$model$cpars$cverbose <- cverbose
maat$model$cpars$diag     <- diag



##################################
# Initialise the MAAT wrapper

# load init xml's or list from init R script
setwd(pdir)

# read user defined values of static variables
if(xml) {

  # read user defined XMLs of static variables
  staticxml   <- paste(mod_obj,'user','static.xml',sep='_')
  init_static <- if(file.exists(staticxml)) readXML(staticxml) else list(NULL)

  # read user defined XMLs of dynamic variables
  dynamicxml   <- paste(mod_obj,'user','dynamic.xml',sep='_')
  init_dynamic <- if(file.exists(dynamicxml)) readXML(dynamicxml) else list(NULL)

  # convert NAs to NULLs
  init_static  <- rapply(init_static,  function(x) if(is.na(x)) NULL else x, how='replace' )
  init_dynamic <- rapply(init_dynamic, function(x) if(is.na(x)) NULL else x, how='replace' )

  # otherwise read init list R script
} else source(initf)


# add init lists to wrapper
maat$init_static  <- init_static
maat$init_dynamic <- init_dynamic


# output static & dynamic values used in simulation 
# - need to move to build to get complete record of model setup (perhaps just need to move static output)
# - could move static config to build function 
print('',quote=F)
print('Write record of static & dynamic run variables:',quote=F)
setwd(odir)
listtoXML(paste(ofname,'setup_static.xml',sep='_'),  'static',  sublist=init_static)
listtoXML(paste(ofname,'setup_dynamic.xml',sep='_'), 'dynamic', sublist=init_dynamic)



##################################
# check process representation functions specified in input exist

search_fnames <-  function(v, ln ) {
  for( c1 in v ) if(!(c1 %in% ls(pos=1)))
    stop(paste('The function: ',c1,' , specified in init', ln ,'does not exist')) else 'exists'
}
print('',quote=F); print('',quote=F)
print('Check static fnames requested exist:', quote=F )
out <- lapply(init_static$fnames,
         function(l) lapply(l, function(l1) if(is.list(l1)) lapply(l1, search_fnames, ln='static') else search_fnames(l1, ln='static' )) )
print('  all static fnames requested exist.', quote=F )

print('', quote=F )
print('Check dynamic fnames requested exist:', quote=F )
if(!is.null(unlist(init_dynamic$fnames))) {
  out <- lapply(init_dynamic$fnames,
           function(l) lapply(l, function(l1) if(is.list(l1)) lapply(l1, search_fnames, ln='dynamic') else search_fnames(l1, ln='dynamic' )) )
}
print('  all dynamic fnames requested exist.', quote=F )


# get mcmc filename, delete history files if not a restart and they exist
if(runtype=='mcmc') {
  maat$wpars$of_name <- paste(ofname, 'mcmc', sep='_' )
  if(is.null(parsinit_mcmc)) {
    setwd(odir)
    hist_file_list <- list.files(pattern=paste0(maat$wpars$of_name, '_history_' ))      
    if(length(hist_file_list)>0) system(paste('rm', paste(hist_file_list,collapse=' ') ))
    setwd(pdir)
  }
}


# write directly to dataf$pars if parsinit specified
mcmc_restart <- F 
if(!is.null(parsinit_mcmc)) {
  if(parsinit_mcmc=='restart') mcmcout <- paste0(maat$wpars$of_name,'.RDS')
  if(is.null(mcmcdir))         mcmcdir <- paste0('results/',date)

  print('',quote=F); print('',quote=F)
  print('Read input from MCMC run output.',quote=F)
  print('directory:',quote=F)
  print(paste(' ',mcmcdir),quote=F)
  print('filename:',quote=F)
  print(paste(' ',mcmcout),quote=F)

  # read and write pars
  setwd(mcmcdir)
  maat$dataf$mcmc_input <- readRDS(mcmcout)
  mcmc_restart_pars_dim <- dim(maat$dataf$mcmc_input$pars_array)

  if(parsinit_mcmc=='restart') {
    print('',quote=F)
    print('MCMC restart:',quote=F)
    print('  no MCMC parameters will be read from init file,',quote=F)
    print('  all set from restarted MCMC run,',quote=F)
    print('  all output will be saved in the directory (see above) of restarted MCMC run.',quote=F)
    mcmc_restart <- T
    
    # debugging
#    print('')
#    print('MCMCinput:')
#    print(maat$dataf$mcmc_input)
#    print('pars arry dim:')
#    print(dim(maat$dataf$mcmc_input$pars_array))
#    print(mcmc_restart_pars_dim)
#    print('')

    # assign starting parameter values
    parsinit <- maat$dataf$mcmc_input$pars_array[,,mcmc_restart_pars_dim[3]]

    # update user-defined MCMC parameters passed from restart
    maat$wpars$mcmc[names(maat$dataf$mcmc_input$wpars$mcmc)] <- maat$dataf$mcmc_input$wpars$mcmc[names(maat$dataf$mcmc_input$wpars$mcmc)]

    # update beginning and end iteration counters
    maat$wpars$mcmc$start_iter      <- mcmc_restart_pars_dim[3] + 1
    #APW can delete: maat$mcmc$start_iter_thin       <- mcmc_restart_pars_dim[3]%%maat$wpars$mcmc$check_iter
    maat$wpars$mcmc$maxiter         <- mcmc_maxiter + mcmc_restart_pars_dim[3]
    maat$wpars$mcmc$maxiter_restart <- mcmc_maxiter
    print('',quote=F)
    print('Start iterations at:',quote=F)
    print(maat$wpars$mcmc$start_iter,quote=F)
    print('End iterations at:',quote=F)
    print(maat$wpars$mcmc$maxiter,quote=F)

    # update MCMC variables passed from restart
    maat$mcmc[names(maat$dataf$mcmc_input$mcmc)] <- maat$dataf$mcmc_input$mcmc[names(maat$dataf$mcmc_input$mcmc)]

    # reload met and obs data
    maat$dataf$met     <- maat$dataf$mcmc_input$met
    maat$dataf$obs     <- maat$dataf$mcmc_input$obs
    maat$dataf$obsse   <- maat$dataf$mcmc_input$obsse
    maat$mcmc$obs_vars <- maat$dataf$mcmc_input$mcmc$obs_vars
   
    # set output to mcmc dir 
    maat$wpars$of_dir <- mcmcdir
   

  } else if(parsinit_mcmc=='ensemble') {
   
    print('forward run from MCMC calibrated parameters',quote=F)
    if(runtype!='factorial') {
      print(paste0('parsinit requested but can only work with factorial runtype, requested runtype:', runtype ),quote=F)
      stop('parsinit requested with incompatible runtype.',quote=F)
    }

    # drop to a matrix
    parsinit <- as.matrix(maat$dataf$mcmc_input$pars_array, nrow=mcmc_restart_pars_dim[1] )

    # sub-sample
    if(parsinit_n>=dim(parsinit)[2]) print('parsinit_n greater than samples in mcmcout')
    parsinit <- parsinit[, sample(dim(parsinit)[2], parsinit_n )]
  }

  # add parsinit to wrapper data structure
  # APW: this could be moved to the generate_pars_ensemble code
  #      would probably help readability
  maat$dataf$pars <- parsinit
}


# output static variables used in simulation
# APW: this will not report pars that have been read by the above if statement
print('',quote=F); print('',quote=F)
print('Write record of static run variables:',quote=F)
setwd(odir)
listtoXML(paste(ofname,'setup_static.xml',sep='_'),  'static',  sublist=init_static)
listtoXML(paste(ofname,'setup_dynamic.xml',sep='_'), 'dynamic', sublist=init_dynamic)



##################################
# Load meteorological / environment dataset
# - each ensemble member is run over this entire dataset
# - not used unless specified

kill <- F
if(!mcmc_restart & !is.null(metdata)) {

  # read user defined met data translator
  setwd(pdir)
  if(file.exists(paste0(mod_obj,'_user_met.xml'))) {
    met_trans <- readXML(paste0(mod_obj,'_user_met.xml'))

    if(any(names(met_trans$env)==mod_obj)) {
      met_trans <- met_trans$env[[which(names(met_trans$env)==mod_obj)]]

    } else {
      print('',quote=F); print('',quote=F)
      print('Met translator file:',quote=F)
      print(paste0(mod_obj,'_user_met.xml'),quote=F)
      print('does not contain list for:',quote=F)
      print(mod_obj,quote=F)
      stop()
    }

  } else {
    print('',quote=F); print('',quote=F)
    print('Met translator file:',quote=F)
    print(paste0(mod_obj,'_user_met.xml'),quote=F)
    print('does not exist in:',quote=F)
    print(pdir,quote=F)
    stop()
  }

  # read metdata file
  print('',quote=F); print('',quote=F)
  print('Met data directory, filename, & translation:' , quote=F )
  print(mdir, quote=F )
  setwd(mdir)
  if(file.exists(metdata)&!kill) {
    print(metdata, quote=F )
    metdffull <- read.csv(metdata,strip.white=T)
    print(head(metdffull), quote=F )
    # ALJ: if doing a Sphagnum simulation
    #      screen out night-time values from met data file
    #      subset it to remove 0's and negative values
    # APW: thinking of a flexible way to do this, it's not that easy. Long term might need to do this in the met data itself
    # sub_idx   <- which(metdffull$EM_PAR_8100_x > 0)
    # metdffull <- metdffull[sub_idx, ]
    #print(head(metdffull), quote=F )

    # order met data in metfile according to that specified in the <mod_obj>_user_met.XML
    # - need to add a trap to catch met data files that do not contain all the data specified in <mod_obj>_user_met.XML
    cols      <- match(unlist(sapply(met_trans,function(l) l)), names(metdffull) )
    metdf     <- metdffull[,cols, drop=F ]
    print('', quote=F )
    print(met_trans, quote=F )
    print(head(metdf), quote=F )
    # if time variable specified put it first and don't rename it
    if('time'%in%names(met_trans)) {
      tcol  <- which(names(met_trans)=='time')
      metdf <- metdf[, c(tcol,c(1:length(metdf))[-tcol]), drop=F ]

      # rename to maat variables as defined in <mod_obj>_user_met.XML and prefix with the model object for compatibility with the configure function
      names(metdf)[1] <- 'time'
      names(metdf)[2:length(metdf)] <- paste(mod_obj,names(met_trans)[-tcol],sep='.')

      # due to matrix data structure a character time vector will mess the numeric matrix up
      # - could do something with posix convention
      if(!is.numeric(metdf[,1])) metdf <- metdf[,-1,drop=F ]

    } else {
      names(metdf) <- paste(mod_obj,names(met_trans),sep='.')
    }

    # add to MAAT object
    print(head(metdf), quote=F )
    maat$dataf$met <- t(as.matrix(metdf))

    # remove met data file
    rm(metdf)
    if(is.null(evaldata)) rm(metdffull) else if(evaldata!='metdata') rm(metdffull)

  } else {
    print('',quote=F); print('',quote=F)
    print('Met data file:',quote=F)
    print(metdata,quote=F)
    print('does not exist in:',quote=F)
    print(mdir,quote=F)
    stop('ERROR: metdata file does not exist')
  }
  setwd(pdir)
}



##################################
# Load evaluation dataset
# - not used unless specified

kill <- F
if(!mcmc_restart & !is.null(evaldata)) {

  # read user defined eval data translator
  setwd(pdir)
  if(file.exists(paste0(mod_obj,'_user_eval.xml'))) {
    eval_trans <- readXML(paste0(mod_obj,'_user_eval.xml'))

    if(any(names(eval_trans$state)==mod_obj)) {
      eval_trans <- eval_trans$state[[which(names(eval_trans$state)==mod_obj)]]

    } else {
      print('',quote=F); print('',quote=F)
      print('Evaluation data translator file:',quote=F)
      print(paste0(mod_obj,'_user_eval.xml'),quote=F)
      print('does not contain list for:',quote=F)
      print(mod_obj,quote=F)
      stop()
    }

  } else {
    print('',quote=F); print('',quote=F)
    print('Evaluation data translator file:',quote=F)
    print(paste0(mod_obj,'_user_eval.xml'),quote=F)
    print('does not exist in:',quote=F)
    print(pdir,quote=F)
    stop()
  }

  # read eval data file
  if(is.null(edir)) edir <- mdir
  print('',quote=F); print('',quote=F)
  print('Eval data directory, filename, & translation:' , quote=F )
  print(edir, quote=F )

  setwd(edir)
  if(evaldata=='metdata') {
    evaldfobs <- metdffull
    rm(metdffull)

  } else if(file.exists(evaldata)&!kill) {
    print(evaldata, quote=F )
    evaldfobs   <- read.csv(evaldata,strip.white=T)
    print(head(evaldfobs), quote=F )

    # check met and eval data sets are of same length
    if(dim(maat$dataf$met)[2]!=dim(evaldfobs)[1])
      stop(paste('Eval data not same length as met data: check equal length of files:', evaldata, metdata ))

  } else {
    print('',quote=F); print('',quote=F)
    print('Eval data file:',quote=F)
    print(evaldata,quote=F)
    print('does not exist in:',quote=F)
    print(evaldir,quote=F)
    stop('ERROR: eval file does not exist')
  }

  # order eval data in evalfile according to that specified in the <mod_obj>_user_eval.XML
  # - need to add a trap to catch eval data files that do not contain all the data specified in <mod_obj>_user_eval.XML
  cols   <- match(unlist(sapply(eval_trans,function(l) l)), names(evaldfobs) )
  evaldf <- evaldfobs[,cols] # if a single column will be dropped to a vector

  # add to MAAT object
  maat$dataf$obs     <- evaldf
  maat$mcmc$obs_vars <- names(eval_trans)
  print('', quote=F )
  print(eval_trans, quote=F )

  # Read standard error for eval data
  if(grepl('_se',mcmc_lklihood) & !evalse) {
    evalse <- T
    print('MCMC likelihood function that requires standard error selected, "evalse" changed to TRUE.', quote=F )
  }
  if(evalse) {
    print('Standard error selected for eval data, variables should be named same as eval data appended by ".se"', quote=F )

    # extract se data
    cols             <- match(paste(unlist(sapply(eval_trans,function(l) l)),'se',sep='.'), names(evaldfobs))
    evaldfse         <- evaldfobs[,cols]
    maat$dataf$obsse <- evaldfse
  }

  # remove eval data file
  rm(evaldfobs)
}



##################################
###  Run MAAT

setwd(pdir)
for(i in 1:5) print('',quote=F)
print(paste('Run MAAT',runtype,':'),quote=F)
st <- system.time(maat$run())

for(i in 1:3) print('',quote=F)
print('MAAT runtime:',quote=F)
print(st,quote=F)

maat$clean()
print('',quote=F)
print('MAAT system memory used:',quote=F)
gc()



### END ###
