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

# run a UQ style ensemble, or if -uq- is false a fully factorial ensemble 
uq      <- F      

# ensemble number for a UQ style ensemble, not used if if -uq- is false 
n       <- 10      

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

# output file format.  supported: rds, csv (default)
of_format <- 'csv'



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
initf  <- if(is.null(runid)) paste(init,'R',sep='.') else paste(init,'_',runid,'.R',sep='')
ofname <- if(is.null(runid)) of_main else paste(runid,of_main,sep='_')



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
# Configure the MAAT wrapper

# build and clone the maat wrapper object
maat <- wrapper_object$.build(model=model)

# define run parameters
maat$wpars$multic       <- multic  
maat$wpars$procs        <- procs   
maat$wpars$UQ           <- uq       
maat$wpars$n            <- n       
maat$model$pars$verbose <- F



##################################
# Initialise the MAAT wrapper

# load init xml's & list 
# if(xml) {

# read default
init_default <- readXML('leaf_default.xml')

# read user defined values of static variables
setwd(pdir)
if(file.exists('leaf_user_static.xml')) {
  init_user    <- readXML('leaf_user_static.xml')
  init_static  <- fuselists(init_default,init_user)
  init_static  <- evalXMLlist(init_static)
} else init_static <- init_default

# write static parameters used in simulation
listtoXML('setup_static.xml','static', sublist=init_static)

# read user defined values of dynamic variables
if(file.exists('leaf_user_dynamic.xml')) {
  init_dynamic <- readXML('leaf_user_dynamic.xml')
  init_dynamic <- evalXMLlist(init_dynamic)
} else init_dynamic <- list(leaf = list(fnames=NA,pars=NA,env=NA))

# otherwise read init list R script (this has not yet been modified to conform with the new init xml/list structure)
# } else source(init)

# add init list to wrapper
# maat$init_ls <- init_ls
maat$init_static  <- init_static
maat$init_dynamic <- init_dynamic



##################################
# Load meteorological and environment dataset
# - each ensemble member is run over this entire dataframe
# - not used unless specified

# # - the data frame must have a minimum of two columns 
# kill <- F
# if(!is.null(metdata)) {
#   setwd(mdir)
#   if(file.exists(metdata)) {
#     maat$dataf$met <- read.csv(metdata,strip.white=T)    
#   } else {
#     print('',quote=F)
#     print('File:',quote=F)
#     print(metdata,quote=F)
#     print('does not exist in:',quote=F)
#     print(mdir,quote=F)
#     kill <- T
#     stop
#   }
#   setwd(pdir)
# }

kill <- F
if(!is.null(metdata)) {
  # read user defined met data translator
  setwd(pdir)
  if(file.exists('leaf_user_met.xml')) {
    met_trans <- readXML('leaf_user_met.xml')
    met_trans <- evalXMLlist(met_trans)
    if(any(names(met_trans)==mod_obj)) met_trans <- met_trans[[which(names(met_trans)==mod_obj)]]$env
    else {
      print('',quote=F)
      print('Met translator file:',quote=F)
      print('leaf_user_met.xml',quote=F)
      print('does not contain list for:',quote=F)
      print(mod_obj,quote=F)      
      kill <- T
    }
      
  } else {
    print('',quote=F)
    print('Met translator file:',quote=F)
    print('leaf_user_met.xml',quote=F)
    print('does not exist in:',quote=F)
    print(pdir,quote=F)
    kill <- T
    stop    
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
    kill <- T
    stop
  }
  setwd(pdir)
}
if(kill) stop



##################################
###  Run MAAT

maat$model$pars$verbose  <- F

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




