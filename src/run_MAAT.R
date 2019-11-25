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

# parameters for SA run
# ensemble number for an SA/UQ style ensemble, not used if -uq- is false
psa_n      <- 10
# coefficient of variation if used in parameter sampling in process sensitivity analysis
coef_var   <- 0.1
# multiplier on process ensemble n for Saltelli ensemble n
salt_nmult <- 100

# parameters for MCMC run
# MCMC likelihood function (options: log, ssquared, ssquared_se, ...)
mcmc_lklihood      <- 'ssquared'
# MCMC outlier handling (options: none, iqr)
mcmc_outlier       <- 'iqr'
# MCMC convergence testing (options: none, Gelman_Rubin)
mcmc_converge      <- 'Gelman_Rubin'
# MCMC option for parameter treatment in bounded search spaces (options: none, bound, reflect, fold)
mcmc_bdry_handling <- 'bound'
# MCMC option for initializing Markov chains with chosen prior distribution (options: uniform, normal, none)
mcmc_prior         <- 'uniform'
# number of MCMC chains to run (minumum = 2 * mcmc_delta + 1)
mcmc_chains        <- 7
# number of iterations / steps in MCMC chain
mcmc_maxiter       <- 1000
# MCMC thinning for posterior, as a proportion
mcmc_thin          <- 0.1
# MCMC thinning for observations, as a proportion
mcmc_thin_obs      <- 1
# MCMC option to assume homoscedastic error in measured observations (else, heteroscedastic)
mcmc_homosced      <- F
# MCMC DREAM number chain pair proposal
mcmc_delta         <- 3
# MCMC DREAM randomization (default value)
mcmc_c_rand        <- 0.01
# MCMC DREAM ergodicicty (default value)
mcmc_c_ergod       <- 1e-12
# MCMC DREAM probability of unit jump rate (probability gamma = 1) (default value)
mcmc_p_gamma       <- 0.2
# MCMC DREAM number of crossover values (default value)
mcmc_n_CR          <- 3
# MCMC option whether or not to adapt probability of selecting crossover values
mcmc_adapt_pCR     <- T
# MCMC option determining how long to adapt crossover selection probabilities, as a proportion
mcmc_CR_burnin     <- 0.1
# MCMC option for checking for convergence and outlier chains every N iterations
mcmc_check_iter    <- 10

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



##################################
# parse command line arguments

print('',quote=F)
print('Read command line arguments',quote=F)
print(commandArgs(T),quote=F)
if(length(commandArgs(T))>=1) {
  for( ca in 1:length(commandArgs(T)) ) {
    eval(parse(text=commandArgs(T)[ca]))
  }
}

print('',quote=F)
if(is.null(srcdir))  stop('srcdir argument not specified but required, check command line arguments to run_MAAT.R')
print(paste('Source directory:',srcdir) ,quote=F)
if(is.null(pdir))    stop('pdir argument not specified but required, check command line arguments to run_MAAT.R')
print(paste('Project directory:',pdir) ,quote=F)
if(is.null(mod_obj)) stop('mod_obj argument not specified but required, check command line arguments to run_MAAT.R')
print(paste('Model:',mod_obj) ,quote=F)
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
initf       <- if(is.null(runid))     paste(init,'R',sep='.') else paste(init,'_',runid,'.R',sep='')
# prefix for output files
ofname      <- if(is.null(runid))     of_main                 else paste(of_main,runid,sep='_')
ofname      <- if(is.null(mod_mimic)) ofname                  else paste(mod_mimic,ofname,sep='_')

# factorial analysis over-rides UQ analysis
#if(!uq) factorial <- T
if(factorial&(uq|mcmc)) {
 uq   <- F
 mcmc <- F
 print('',quote=F)
 print(paste('Both factorial and UQ or MCMC run specified: Factorial ensemble will be run'),quote=F)
}

# select run/ensemble type and output
runtype <-
  if(factorial) 'factorial'      else
  if(procSA)    'SAprocess_ye'   else
  if(salt)      'SApar_saltelli' else
  if(mcmc)      paste0('mcmc_',mcmc_type)

print(paste('Run type:',runtype,'selected') ,quote=F)

if((uq|mcmc)&of_format!='rds') {
  of_format <- 'rds'
  print('',quote=F)
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

# define MCMC run parameters
maat$wpars$mcmc               <- mcmc
maat$wpars$mcmc_type          <- mcmc_type
maat$wpars$mcmc_lklihood      <- mcmc_lklihood
maat$wpars$mcmc_outlier       <- mcmc_outlier
maat$wpars$mcmc_bdry_handling <- mcmc_bdry_handling
maat$wpars$mcmc_prior         <- mcmc_prior
maat$wpars$mcmc_converge      <- mcmc_converge
maat$wpars$mcmc_chains        <- mcmc_chains
maat$wpars$mcmc_maxiter       <- mcmc_maxiter
maat$wpars$mcmc_thin          <- mcmc_thin
maat$wpars$mcmc_thin_obs      <- mcmc_thin_obs
maat$wpars$mcmc_homosced      <- mcmc_homosced
maat$wpars$mcmc_delta         <- mcmc_delta
maat$wpars$mcmc_c_rand        <- mcmc_c_rand
maat$wpars$mcmc_c_ergod       <- mcmc_c_ergod
maat$wpars$mcmc_p_gamma       <- mcmc_p_gamma
maat$wpars$mcmc_n_CR          <- mcmc_n_CR
maat$wpars$mcmc_adapt_CR      <- mcmc_adapt_pCR
maat$wpars$mcmc_CR_burnin     <- mcmc_CR_burnin
maat$wpars$mcmc_check_iter    <- mcmc_check_iter



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
# build maat wrapper and model objects

maat$build(mod_mimic=mod_mimic, mod_out=mod_out )

# set debugging flags
maat$model$cpars$verbose  <- verbose
maat$model$cpars$cverbose <- cverbose



##################################
# check process representation functions specified in input exist

search_fnames <-  function(v, ln ) {
  for( c1 in v ) if(!(c1 %in% ls(pos=1))) stop(paste('The function: ',c1,' , specified in init', ln ,'does not exist')) else 'exists'
}

print('',quote=F)
print('Check static fnames requested exist:',quote=F)
out <- lapply(init_static$fnames,
         function(l) lapply(l, function(l1) if(is.list(l1)) lapply(l1, search_fnames, ln='static') else search_fnames(l1, ln='static' )) )
print('  all static fnames requested exist',quote=F)

print('',quote=F)
print('Check dynamic fnames requested exist:',quote=F)
if(!is.null(unlist(init_dynamic$fnames))) {
  out <- lapply(init_dynamic$fnames,
           function(l) lapply(l, function(l1) if(is.list(l1)) lapply(l1, search_fnames, ln='dynamic') else search_fnames(l1, ln='dynamic' )) )
}
print('  all dynamic fnames requested exist',quote=F)



##################################
# Load meteorological and environment dataset
# - each ensemble member is run over this entire dataframe
# - not used unless specified

kill <- F
if(!is.null(metdata)) {
  # read user defined met data translator
  setwd(pdir)
  if(file.exists(paste0(mod_obj,'_user_met.xml'))) {
    met_trans <- readXML(paste0(mod_obj,'_user_met.xml'))
    if(any(names(met_trans$env)==mod_obj)) met_trans <- met_trans$env[[which(names(met_trans$env)==mod_obj)]]
    else {
      print('',quote=F)
      print('Met translator file:',quote=F)
      print(paste0(mod_obj,'_user_met.xml'),quote=F)
      print('does not contain list for:',quote=F)
      print(mod_obj,quote=F)
      stop()
    }

  } else {
    print('',quote=F)
    print('Met translator file:',quote=F)
    print(paste0(mod_obj,'_user_met.xml'),quote=F)
    print('does not exist in:',quote=F)
    print(pdir,quote=F)
    stop()
  }

  # read metdata file
  print('', quote=F )
  print('Met data directory & filename:' , quote=F )
  print(mdir, quote=F )
  setwd(mdir)
  if(file.exists(metdata)&!kill) {
    print(metdata, quote=F )
    metdf <- read.csv(metdata,strip.white=T)

    # ALJ: if doing a Sphagnum simulation
    #      screen out night-time values from met data file
    #      subset it to remove 0's and negative values
    # sub_idx <- which(metdf$EM_PAR_8100_x > 0)
    # metdf <- metdf[sub_idx, ]

    print(head(metdf), quote=F )

    ###################################
    # future work: add eval data to model object - total hack for now
    # for Sphagnum simulations
    # maat$dataf$obs    <- metdf$GPP.PAR.ecor.real
    # maat$dataf$obsse  <- metdf$GPP.PAR.ecor.real.se
    # for ACi simulations
    maat$dataf$obs <- metdf$A
    ###################################

    # order met data in metfile according to that specified in the <mod_obj>_user_met.XML
    # - need to add a trap to catch met data files that do not contain all the data specified in <mod_obj>_user_met.XML
    cols  <- match(unlist(sapply(met_trans,function(l) l)),names(metdf))
    metdf <- metdf[,cols, drop=F ]

    # if time variable specified put it first and don't rename it
    if('time'%in%names(met_trans)) {
      tcol  <- which(names(met_trans)=='time')
      metdf <- metdf[, c(tcol,c(1:length(metdf))[-tcol]), drop=F ]

      # rename to maat variables as defined in <mod_obj>_user_met.XML and prefix with the model object for compatibility with the configure function
      names(metdf)[1] <- 'time'
      names(metdf)[2:length(metdf)] <- paste(mod_obj,names(met_trans)[-tcol],sep='.')
    } else {
      names(metdf) <- paste(mod_obj,names(met_trans),sep='.')
    }

    ###################################
    # add to MAAT object
    # future work: generalize this
    # for sphagnum simulations
    # maat$dataf$met <- t(as.matrix(metdf[,-1]))
    # for ACii simulations
    maat$dataf$met <- t(as.matrix(metdf))
    ###################################

    # remove met data file
    rm(metdf)

  } else {
    print('',quote=F)
    print('Met data file:',quote=F)
    print(metdata,quote=F)
    print('does not exist in:',quote=F)
    print(mdir,quote=F)
    stop()
  }
  setwd(pdir)
}



##################################
# Load evaluation dataset
# - not used unless specified

kill <- F
if(!is.null(evaldata)&F) {
  # read user defined eval data translator
  setwd(pdir)
  if(file.exists(paste0(mod_obj,'_user_eval.xml'))) {
    eval_trans <- readXML(paste0(mod_obj,'_user_eval.xml'))
    if(any(names(eval_trans$state)==mod_obj)) eval_trans <- eval_trans$env[[which(names(eval_trans$state)==mod_obj)]]
    else {
      print('',quote=F)
      print('Evaluation data translator file:',quote=F)
      print(paste0(mod_obj,'_user_eval.xml'),quote=F)
      print('does not contain list for:',quote=F)
      print(mod_obj,quote=F)
      stop()
    }

  } else {
    print('',quote=F)
    print('Evaluation data translator file:',quote=F)
    print(paste0(mod_obj,'_user_eval.xml'),quote=F)
    print('does not exist in:',quote=F)
    print(pdir,quote=F)
    stop()
  }

  # read eval data file
  if(is.null(edir)) edir <- mdir
  print('', quote=F )
  print('Eval data directory & filename:' , quote=F )
  print(edir, quote=F )
  setwd(edir)

  if(evaldata=='met') evaldata <- metdata
  if(file.exists(evaldata)&!kill) {
    print(evaldata, quote=F )
    evaldf   <- read.csv(evaldata,strip.white=T)
    print(head(evaldf), quote=F )

    # check met and eval data sets are of same length
    if(dim(maat$dataf$met)[2]!=dim(evaldf)[2])
      stop(paste('Eval data not same length as met data: check equal length of files:', evaldata, metdata ))

    # order eval data in evalfile according to that specified in the <mod_obj>_user_eval.XML
    # - need to add a trap to catch eval data files that do not contain all the data specified in <mod_obj>_user_eval.XML
    cols   <- match(unlist(sapply(eval_trans,function(l) l)),names(evaldf))
    evaldf <- evaldf[,cols]

    # add to MAAT object
    maat$dataf$obs <- evaldf

    # Read standard error for eval data
    if(evalse) {
      print('Standard error selected for eval data, variables should be named same as eval data appended by ".se"')

      # extract se data
      cols     <- match(unlist(sapply(eval_trans,function(l) l)),paste(names(evaldf),'se',sep='.'))
      evaldfse <- evaldf[,cols]

      # add to MAAT object
      maat$dataf$obsse <- evaldfse
    }

    # remove eval data file
    rm(evaldf)

  } else {
    print('',quote=F)
    print('Eval data file:',quote=F)
    print(evaldata,quote=F)
    print('does not exist in:',quote=F)
    print(evaldir,quote=F)
    stop()
  }
  setwd(pdir)
}



##################################
###  Run MAAT

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
