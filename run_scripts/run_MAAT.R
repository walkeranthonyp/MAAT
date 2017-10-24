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
srcdir  <- "#SOURCEDIR#" 
# project directory (full path)
pdir    <- "#PROJECTDIR#"
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
factorial  <- F

# run an SA/UQ style ensemble, or if -uq- is false a fully factorial ensemble 
uq         <- T      

# types of SA/UQ run
# process SA
procSA     <- T
# Saltelli Sobol SA
salt       <- T

# parameters for SA/UQ run
# ensemble number for an SA/UQ style ensemble, not used if -uq- is false 
psa_n      <- 10
# coefficient of variation if used in parameter sampling in process sensitivity analysis
coef_var   <- 0.1 
# multiplier on process ensemble n for Saltelli ensemble n
salt_nmult <- 100      

# run options
# model object to use, currently leaf or canopy
mod_obj    <- 'leaf'

# meteorological data file name
metdata    <- NULL 

# pass code snippets as strings to generate parameter samples, see init_MAAT.R and wrapper for use
eval_strings <- F 

# initialise in the configuration of the below specified model
# options are: clm40, 
mod_mimic    <- NULL

# static and dynamic initialisation files are XMLs, if false init file is the R script named below
xml          <- F

# initialisation data file name if not an XML
init         <- 'init_MAAT' 

# run i.d. - used as suffix/prefix for in/out files
runid        <- NULL 

# basic output file name
of_main      <- 'out'

# model output switch 
mod_out      <- 'run'

# output file format.  supported: rds, csv (default)
of_format    <- 'csv'



#############################################################################################
### If users want to edit default options, 
### the above object assignments can be modified once this script has been copied to the project directory
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
if(is.null(odir)) {
  date   <- Sys.Date()
  odir1  <- paste(pdir,'results',sep='/')
  odir   <- paste(odir1,date,sep='/')
  if(!file.exists(odir1)) dir.create(odir1)
  if(!file.exists(odir))  dir.create(odir)
}

# create input/output filenames
initf   <- if(is.null(runid))     paste(init,'R',sep='.')       else paste(init,'_',runid,'.R',sep='')
ofname  <- if(is.null(runid))     of_main                       else paste(runid,of_main,sep='_')
ofname  <- if(is.null(mod_mimic)) ofname                        else paste(ofname,mod_mimic,sep='_')
sofname <- if(is.null(runid))     paste(of_main,'salt',sep='_') else paste(runid,of_main,'salt',sep='_')



###################################################################
### start program

##################################
# load MAAT wrapper from source
setwd(srcdir)
source('wrapper_object.R')
# load MAAT object(s) from source
setwd(paste('./',mod_obj,sep=''))
source(paste(mod_obj,'object.R',sep='_'))
# read default model setup
init_default <- readXML(paste(mod_obj,'default.xml',sep='_'))
# read model mimic setup
if(!is.null(mod_mimic)) {
  setwd('./mimic_xmls')
  init_mimic   <- readXML(paste(mod_obj,'_',mod_mimic,'.xml',sep=''))
  init_default <- fuselists(init_default,init_mimic) 
}



##################################
# Configure and initialise the MAAT model

# build and clone the model object
# model <- get(paste(mod_obj,'object',sep='_'))$.build()



##################################
# Configure the MAAT wrapper

# clone and build the maat wrapper and model object
maat <- as.proto(wrapper_object$as.list()) 
maat$build(model=paste(mod_obj,'object',sep='_'))
rm(wrapper_object)

# factorial analysis over-rides UQ analysis
if(!uq) factorial <- T
if(factorial&uq) {
 uq <- F
 print('',quote=F)
 print(paste('Both factorial and UQ run specified: Factorial ensemble will be run') ,quote=F)
}

# define run parameters
maat$wpars$multic        <- multic  
maat$wpars$procs         <- procs   
maat$wpars$UQ            <- uq       
maat$wpars$n             <- psa_n       
maat$wpars$coef_var      <- coef_var       
maat$wpars$nmult         <- salt_nmult       
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
  init_static <- 
    if(file.exists(staticxml)) readXML(staticxml)
    else                       list(leaf = list(fnames=NA,pars=NA,env=NA))
  
  # read user defined XMLs of dynamic variables
  dynamicxml   <- paste(mod_obj,'user','dynamic.xml',sep='_')
  init_dynamic <-
    if(file.exists(dynamicxml)) readXML(dynamicxml)
    else                        list(leaf = list(fnames=NA,pars=NA,env=NA))
  
  # otherwise read init list R script
} else source(initf)

# combine default values with user defined static values
inits <- if(exists('init_static')) fuselists(init_default,init_static) else init_default

# add init lists to wrapper
maat$init_static  <- inits
maat$init_dynamic <- init_dynamic

# write static parameters used in simulation
print('',quote=F)
print('Write record of static run variables',quote=F)
setwd(odir)
listtoXML(paste(runid,'setup_static.xml',sep='_'),  'static',  sublist=inits)
listtoXML(paste(runid,'setup_dynamic.xml',sep='_'), 'dynamic', sublist=init_dynamic)



##################################
# Load meteorological and environment dataset
# - each ensemble member is run over this entire dataframe
# - not used unless specified

kill <- F
if(!is.null(metdata)) {
  # read user defined met data translator
  setwd(pdir)
  if(file.exists('leaf_user_met.xml')) {
    met_trans <- readXML('leaf_user_met.xml')
    # met_trans <- evalXMLlist(met_trans)
    if(any(names(met_trans)==mod_obj)) met_trans <- met_trans[[which(names(met_trans)==mod_obj)]]$env
    else {
      print('',quote=F)
      print('Met translator file:',quote=F)
      print('leaf_user_met.xml',quote=F)
      print('does not contain list for:',quote=F)
      print(mod_obj,quote=F)      
      stop()
    }
      
  } else {
    print('',quote=F)
    print('Met translator file:',quote=F)
    print('leaf_user_met.xml',quote=F)
    print('does not exist in:',quote=F)
    print(pdir,quote=F)
    stop()
  }

  # read metdata file
  setwd(mdir)
  if(file.exists(metdata)&!kill) {
    metdf <- read.csv(metdata,strip.white=T)  
    
    # order met data in metfile according to that specified in the leaf_user_met.XML 
    # - need to add a trap to catch met data files that do not contain all the data specified in leaf_user_met.XML 
    cols  <- match(unlist(sapply(met_trans,function(l) l)),names(metdf))
    metdf <- metdf[,cols] 
        
    # if time variable specified put it first and don't rename it
    if('time'%in%names(met_trans)) {
      tcol  <- which(names(met_trans)=='time')
      metdf <- metdf[, c(tcol,c(1:length(metdf))[-tcol]) ]
      
      # rename to maat variables as defined in leaf_user_met.XML and prefix with the model object for compatibility with the configure function
      names(metdf)[1] <- 'time'               
      names(metdf)[2:length(metdf)] <- paste(mod_obj,names(met_trans)[-tcol],sep='.')
    } else {
      names(metdf) <- paste(mod_obj,names(met_trans),sep='.')
    }
        
    # add to MAAT object
    maat$dataf$met <- metdf 
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
  print('MAAT system memory used for process SA:',quote=F)
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
  print('MAAT system memory used for process SA:',quote=F)
  gc()
}


