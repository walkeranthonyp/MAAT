################################
#
# MAAT Leaf Model - example run script
# 
# AWalker (walkerap@ornl.gov) 
# December 2015
#
################################


###################################################################
### The MAAT model (Multi-Assumption Architecture & Testbed)
###################################################################

# This script runs the leaf or canopy photosynthesis version of the MAAT

# This model is written for ultimate flexibility in function used to represent processes, parameters, and environmental driving data
# The model can be configured in a vast number of possible ways and can be run many times, varying functions, parameters, and environmental driving data during runtime of the model

###
# This script initialises the MAAT and runs it in the following steps:
# - set default arguments
# - parse command line arguments
# - set arguments that depend on other arguments
# - load MAAT objects from source
# - load init script
# - Configure and initialise the MAAT model
# - Configure and initialise the MAAT wrapper
# - Run the MAAT
# - Write output

###################################################################

rm(list=ls())

##################################
# command line options and set default arguments

# directory paths - set these on the command line or modify these here
# source directory (full path) 
# - must be added (here or as a commandline option to this script) before model will run
# - can be modified to target source copied to a static directory rather that a reposiotory with version control 
srcdir  <- NULL 
# project directory (full path)
pdir    <- "#PROJECTDIR#"
# meteorological data directory (full path)
mdir    <- NULL 
# evalutaion data directory (full path)
edir    <- NULL 

# wrapper object options
# multicore the ensemble
multic  <- F  

# number of cores to use if above is true
procs   <- 4   

# run an SA/UQ style ensemble, or if -uq- is false a fully factorial ensemble 
uq      <- T      

# ensemble number for an SA/UQ style ensemble, not used if -uq- is false 
psa_n   <- 10      

# types of SA/UQ run
# process SA
procSA  <- T
# Saltelli Sobol SA
salt    <- T

# multiplier on process ensemble n for Saltelli ensemble n
salt_nmult <- 100      

# run options
# model object to use, currently leaf or canopy
mod_obj <- 'leaf'

# meteorological data file name
metdata <- NULL 

# initialisation data file name
init    <- 'init_MAAT' 

# run i.d. - used as suffix/prefix for in/out files
runid   <- NULL 

# basic output file name
of_main <- 'out'

# output file format
of_format <- 'rds'



#############################################################################################

### DO NOT EDIT BELOW THIS LINE 

#############################################################################################

##################################
# parse command line arguments   
# - any one of the above objects can be specified as a command line argument using the syntax:
# - Rscript <nameofthisscript> "<object1><-<value1>" "<object2><-<value2>"
# - e.g. Rscript example_MAAT.R "dir<-'/home/alp/'" "multic<-T"

print('',quote=F)
print('Read command line arguments',quote=F)
if(length(commandArgs(T))>=1) {
  for( ca in 1:length(commandArgs(T)) ) {
    eval(parse(text=commandArgs(T)[ca]))
  }
}

print('',quote=F)
print(paste('Run ID:',runid) ,quote=F)
if(uq&of_format!='rds') {
  of_format <- 'rds'
  print('',quote=F)
  print(paste('of_format changed to rds due to high output volume with SA/UQ ensembles'), quote=F)
}

# set default values if not specified on command line
# - these are set after parsing command line arguments as they depend on other arguments that could be set on the command line
if(is.null(of_main)) of_main <- proj

# create output directory
date   <- Sys.Date()
odir1  <- paste(pdir,'results',sep='/')
odir   <- paste(odir1,date,sep='/')
if(!file.exists(odir1)) dir.create(odir1)
if(!file.exists(odir))  dir.create(odir)

# create input/output filenames
initf   <- if(is.null(runid)) paste(init,'R',sep='.') else paste(init,'_',runid,'.R',sep='')
ofname  <- if(is.null(runid)) of_main else paste(runid,of_main,sep='_')
sofname <- if(is.null(runid)) of_main else paste(runid,of_main,'salt',sep='_')



###################################################################
### start program

##################################
# load MAAT objects from source
setwd(srcdir)
source('wrapper_object.R')
source(paste(mod_obj,'object.R',sep='_'))



##################################
# Configure and initialise the MAAT model

# build and clone the model object
model <- get(paste(mod_obj,'object',sep='_'))$.build()



##################################
# Configure and initialise the MAAT wrapper

# build and clone the maat wrapper object
maat <- wrapper_object$.build(model=model)

# define run parameters
maat$wpars$multic       <- multic  
maat$wpars$procs        <- procs   
maat$wpars$UQ           <- uq       
maat$wpars$n            <- psa_n       
maat$wpars$nmult        <- salt_nmult       
maat$model$pars$verbose <- F

# define number first and second loop parameter samples
# - deprecated, now done automatically during runtime
# n  <- psa_n
# nB <- 3*psa_n^2

# load init list 
setwd(pdir)
source(initf)

# add init list to wrapper
maat$init_ls <- init_ls


##################################
# Load meteorological and environment dataset
# - each ensemble member is run over this entire dataframe
# - not used unless specified
# - the data frame must have a minimum of two columns 
kill <- F
if(!is.null(metdata)) {
  setwd(mdir)
  if(file.exists(metdata)) {
    maat$dataf$met <- read.csv(metdata,strip.white=T)    
  } else {
    print('',quote=F)
    print('File:',quote=F)
    print(metdata,quote=F)
    print('does not exist in:',quote=F)
    print(mdir,quote=F)
    kill <- T
    stop
  }
  setwd(pdir)
}

if(kill) stop



##################################
###  Run MAAT

if(procSA) {
  maat$model$pars$verbose  <- F
  maat$wpars$UQtype <- 'ye'
  
  st <- system.time(
    maat$run()
  )
  print('',quote=F)
  print('MAAT runtime:',quote=F)
  print(st,quote=F)
  
  # process & record output
  setwd(odir)
  df_out <- maat$output()
  write_to_file(df_out,ofname,type=of_format)  
  
  # delete memory hungry dataframes etc
  rm(df_out)
  maat$clean()
  gc()
  
}



##################################
### run Saltelli algorithm for Sobol analysis if requested

if(salt) {
  
  # reconfigure ensemble parameters
  maat$wpars$UQtype <- 'saltelli'
  # maat$wpars$n     <- psa_n * sobol_nmult
  # n <- nB          <- 2 * psa_n * sobol_nmult
  
  # reconfigure parameter samples 
  # setwd(pdir)
  # source(initf)
  
  # add reconfigured init list to wrapper
  # maat$init_ls <- init_ls
  
  # run reconfigured MAAT
  print('',quote=F)
  print('',quote=F)
  print('Run Saltelli Sobol',quote=F)
  st <- system.time(
    maat$run()
  )
  print('',quote=F)
  print('MAAT Saltelli Sobol runtime:',quote=F)
  print(st,quote=F)
  
  # write output
  salt_out <- maat$output_saltelli()
  setwd(odir)
  write_to_file(salt_out,sofname,type=of_format)
  
}


