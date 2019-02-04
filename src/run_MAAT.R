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
# - load MAAT model objects from source
# - load init scripts
# - Configure and initialise the MAAT model
# - Configure and initialise the MAAT wrapper
# - Run MAAT
# - Write output

###################################################################

#rm(list=ls())

# any one of the below objects to line 116 can be specified as a single string command line argument to this script 
# the sub-arguments in the string (separated by a space) are interpreted individually as R code.

#       Rscript run_MAAT.R "object1<-value1 object2<-value2"
#  e.g. Rscript run_MAAT.R "dir<-'/home/alp/' multic<-T"

##################################
# command line options and defaults

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

# model object to use, currently leaf or gwater_rt
mod_obj <- NULL

# wrapper object options
# multicore the ensemble
multic  <- F  

# number of cores to use if above is true
procs   <- 4   

# run an emsemble that combines variables in factorial
# - if set to TRUE this will over-ride a UQ analysis  
factorial  <- T

# run an SA/UQ style ensemble, or if -uq- is false a fully factorial ensemble 
uq         <- T      

# types of SA/UQ run
# process SA
procSA     <- T
# Saltelli Sobol SA
salt       <- T
# MCMC parameter estimation 
mcmc       <- F

# parameters for SA run
# ensemble number for an SA/UQ style ensemble, not used if -uq- is false 
psa_n      <- 10
# coefficient of variation if used in parameter sampling in process sensitivity analysis
coef_var   <- 0.1 
# multiplier on process ensemble n for Saltelli ensemble n
salt_nmult <- 100      

# parameters for MCMC run
# type of MCMC
mcmc_type     <- 'demc'
# number of MCMC chains to run (min 2x number of parameters estimated)
mcmc_chains   <- 10
# number of iterations / steps in MCMC chain
mcmc_maxiter  <- 100
# MCMC burn in to discard as a proportion of mcmc_maxiter, or number of iterations if converged sooner
mcmc_burnin   <- 0.5
# MCMC thinning for posterior, as a proportion 
mcmc_thin     <- 0.1
# MCMC thinning for observations, as a proportion 
mcmc_thin_obs <- 1
# MCMC option to assume homoscedastic error in measured observations (else, heteroscedastic)
mcmc_homosced <- F

# run options
# meteorological data file name
metdata       <- NULL 

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
initf       <- if(is.null(runid))     paste(init,'R',sep='.')       else paste(init,'_',runid,'.R',sep='')
# prefix for output files
ofname      <- if(is.null(runid))     of_main                       else paste(of_main,runid,sep='_')
ofname      <- if(is.null(mod_mimic)) ofname                        else paste(mod_mimic,ofname,sep='_')
# prefix for saltelli sensitivity output files
sofname     <- if(is.null(runid))     paste(of_main,'salt',sep='_') else paste(runid,of_main,'salt',sep='_')



###################################################################
### start program

##################################
# load MAAT wrapper from source
setwd(srcdir)
source('wrapper_object.R')
# load MAAT object(s) from source
setwd(paste('system_models',mod_obj,sep='/'))
mod_object <- paste(mod_obj,'object',sep='_')
source(paste0(mod_object,'.R'))



##################################
# Configure the MAAT wrapper

# clone and build the maat wrapper and model object
maat         <- as.proto(wrapper_object$as.list()) 
maat$build(model=mod_object, mod_mimic=mod_mimic )
rm(wrapper_object)
rm(mod_object)

# factorial analysis over-rides UQ analysis
if(!uq) factorial <- T
if(factorial&uq) {
 uq <- F
 print('',quote=F)
 print(paste('Both factorial and UQ run specified: Factorial ensemble will be run') ,quote=F)
}

if(uq&of_format!='rds') {
  of_format <- 'rds'
  print('',quote=F)
  print(paste('of_format changed to rds due to high output volume with SA/UQ ensembles'), quote=F)
}

# define run parameters
maat$wpars$multic        <- multic  
maat$wpars$procs         <- procs   
maat$wpars$UQ            <- uq       
maat$wpars$n             <- psa_n       
maat$wpars$coef_var      <- coef_var       
maat$wpars$nmult         <- salt_nmult       
maat$wpars$mcmc_chains   <- mcmc_chains       
maat$wpars$mcmc_maxiter  <- mcmc_maxiter       
maat$wpars$mcmc_burnin   <- mcmc_burnin       
maat$wpars$mcmc_thin     <- mcmc_thin       
maat$wpars$mcmc_thin_obs <- mcmc_thin_obs       
maat$wpars$mcmc_homosced <- mcmc_homosced
maat$wpars$eval_strings  <- eval_strings       
maat$model$cpars$verbose <- F
maat$model$cpars$output  <- mod_out



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


# add init lists to wrapper
maat$init_static  <- init_static
maat$init_dynamic <- init_dynamic

# output static parameters used in simulation
print('',quote=F)
print('Write record of static run variables:',quote=F)
setwd(odir)
listtoXML(paste(ofname,'setup_static.xml',sep='_'),  'static',  sublist=init_static)
listtoXML(paste(ofname,'setup_dynamic.xml',sep='_'), 'dynamic', sublist=init_dynamic)



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
    # met_trans <- evalXMLlist(met_trans)
    #if(any(names(met_trans)==mod_obj)) met_trans <- met_trans[[which(names(met_trans)==mod_obj)]]$env
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
    
    ###################################
    # add eval data to model object - total hack for now   
    maat$dataf$obs    <- metdf$GPP.PAR.ecor.real    
    maat$dataf$obsse  <- metdf$GPP.PAR.ecor.real.se 
    ###################################

    # order met data in metfile according to that specified in the <mod_obj>_user_met.XML 
    # - need to add a trap to catch met data files that do not contain all the data specified in <mod_obj>_user_met.XML 
    cols  <- match(unlist(sapply(met_trans,function(l) l)),names(metdf))
    metdf <- metdf[,cols] 
        
    # if time variable specified put it first and don't rename it
    if('time'%in%names(met_trans)) {
      tcol  <- which(names(met_trans)=='time')
      metdf <- metdf[, c(tcol,c(1:length(metdf))[-tcol]) ]
      
      # rename to maat variables as defined in <mod_obj>_user_met.XML and prefix with the model object for compatibility with the configure function
      names(metdf)[1] <- 'time'               
      names(metdf)[2:length(metdf)] <- paste(mod_obj,names(met_trans)[-tcol],sep='.')
    } else {
      names(metdf) <- paste(mod_obj,names(met_trans),sep='.')
    }
    
    # add to MAAT object
    maat$dataf$met <- metdf 
      
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
###  Run MAAT

# run factorial MAAT, this is a standard setup combining variables in factorial
if(factorial) {
  maat$model$pars$verbose  <- F
  
  for(i in 1:5) print('',quote=F)
  print('Run Factorial',quote=F)
  st <- system.time(
    maat$run()
  )
  
  for(i in 1:3) print('',quote=F)
  print('MAAT runtime:',quote=F)
  print(st,quote=F)
  
  # process & record output
  setwd(odir)
  df_out <- maat$output()
  write_to_file(df_out,ofname,type=of_format)  
  
  rm(df_out)
  maat$clean()
  print('',quote=F)
  print('MAAT system memory used for factorial run:',quote=F)
  gc()
}



### run Ye algorithm for process sensitivity  analysis if requested
if(procSA&uq) {
  maat$model$pars$verbose  <- F
  maat$wpars$UQtype <- 'ye'
  
  for(i in 1:5) print('',quote=F)
  print('Run Process SA',quote=F)
  st <- system.time(
    maat$run()
  )
  
  for(i in 1:3) print('',quote=F)
  print('MAAT runtime:',quote=F)
  print(st,quote=F)
  
  maat$clean()
  print('',quote=F)
  print('MAAT system memory used for process SA:',quote=F)
  gc()
}



### run Saltelli algorithm for Sobol sensitivity  analysis if requested
if(salt&uq) {
  maat$model$pars$verbose  <- F
  maat$wpars$UQtype <- 'saltelli'

  # run MAAT
  for(i in 1:5) print('',quote=F)
  print('Run Saltelli Sobol',quote=F)
  st <- system.time(
    maat$run()
  )

  for(i in 1:3) print('',quote=F)
  print('MAAT Saltelli Sobol runtime:',quote=F)
  print(st,quote=F)

  maat$clean()
  print('',quote=F)
  print('MAAT system memory used for Saltelli SA:',quote=F)
  gc()
}

### run MCMC algorithm parameter estimation if requested
if(mcmc&uq) {
  maat$model$pars$verbose  <- F
  maat$wpars$UQtype <- 'mcmc'

  # run MAAT
  for(i in 1:5) print('',quote=F)
  print('Run MCMC',quote=F)
  st <- system.time(
    maat$run()
  )

  for(i in 1:3) print('',quote=F)
  print('MAAT MCMC runtime:',quote=F)
  print(st,quote=F)

  maat$clean()
  print('',quote=F)
  print('MAAT system memory used for MCMC:',quote=F)
  gc()
}



### END ###
