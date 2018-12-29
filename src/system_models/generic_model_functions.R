################################
#
# Generic system model object functions 
# 
# AWalker Dec 2018
#
################################



# GENERIC SYSTEM MODEL OBJECT FUNCTIONS
###############################################################################

# build functions 
###########################################################################

build_no_child <- function(., mod_mimic=NULL, ... ) {

  # read default model setup for highest level model
  source('../../functions/general_functions.R')
  init_default <- readXML(paste(.$name,'default.xml',sep='_'))
 
  # read model mimic setup
  if(!is.null(mod_mimic)) {
    setwd('mimic_xmls')
    print(paste(.$name, 'mimic:', mod_mimic ))
    init_mimic   <- readXML(paste(.$name,'_',mod_mimic,'.xml',sep=''))
    init_default <- fuselists(init_default,init_mimic)
    setwd('..')
  }

  # assign default and mod mimic values to data structure
  .$configure(vlist='fnames', df=unlist(init_default$fnames)) 
  .$configure(vlist='pars',   df=unlist(init_default$pars)) 
  .$configure(vlist='env',    df=unlist(init_default$env)) 

  # generate methods list within object
  .$fns <- proto(envir=., expr = {
      fns = list()
  })

  # assign methods to methods list
  .$fns$test <- function(.) print(.$fnames$sys)
  .$fns$sys  <- get(.$fnames$sys)
}



# main run functions
###########################################################################

run <- function(.) {
 
  # call system model 
  .$fns$test()
  .$fns$sys()

  # print to screen
  if(.$cpars$verbose) print(.$state)

  # output
  .$output()
} 



# run met function
# wrapper function called from an lapply function to run model over every row of a meteorology dataframe
###########################################################################

run_met <- function(.,l) {
  # assumes that each row of the dataframe are sequential
  # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
  # typically used to run the model using data collected at a specific site and to compare against observations
  
  # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
  # any "env" variables specified in the "dataf$env" dataframe but also specified in .$dataf$met will be overwritten by the .$dataf$met values 
 
  # met data assignment
  .$configure(vlist='env', df=.$dataf$met[l,] )
 
  # run model
  .$run()              
}



# output functions
###########################################################################

output <- function(.){
    
  if(.$cpars$output=='full') {
    
    lout <- unlist(c(.$state, .$state_pars ))
    
  } 
    
  lout
}    


# configure functions
###########################################################################

configure_no_child <- function(., vlist, df, o=T ) {
  # This function is called from any of the run functions, or during model initialisation
  # - sets the values within vlist (i.e. .$fnames / .$pars / .$env / .$state ) to the values passed in df 

  # split variable names at . 
  listnames <- vapply( strsplit(names(df),'.', fixed=T), function(cv) {cv3<-character(3); cv3[1:length(cv)]<-cv; t(cv3)}, character(3) )

  # df subscripts for model object 
  mss <- which(listnames[1,]==.$name)

  # variable list subscripts in model object data structure 
  vlss   <- match(listnames[2,mss], names(.[[vlist]]) )

  # remove NAs in vlss from vlss and mss
  if(any(is.na(vlss))) {
    mss  <- mss[-which(is.na(vlss))]
    vlss <- vlss[-which(is.na(vlss))]
  }

  # df subscripts for sublist variables (slmss) and non-sublist variables (nslmss) 
  slss   <- which(listnames[3,mss]!='')
  if(length(slss)>0) {
    slmss  <- mss[slss] 
    nslmss <- mss[-slss]
    vlss   <- vlss[-slss] 
  } else {
    slmss  <- NULL 
    nslmss <- mss 
  }

  # print configure setup if requested
  if(.$cpars$cverbose) {
    print('', quote=F )
    print(paste(.$name,'configure:'), quote=F )
    print(df, quote=F )
    print(listnames, quote=F )
    print(mss, quote=F )
    print(slmss, quote=F )
    print(nslmss, quote=F )
    print(vlss, quote=F )
    print(which(is.na(vlss)), quote=F )
    print(.[[vlist]], quote=F )
  }

  # assign UQ variables
  #print(paste(.$name,'conf:', vlist, names(df), df ))
  if(length(slss)>0)    vapply( slmss, .$configure_sublist, numeric(1), vlist=vlist, df=df ) 
  if(length(nslmss)>0) .[[vlist]][vlss] <- df[nslmss]
  #print(paste(df[vlmoss],.[[vlist]][vlss])) 
}
 

# configure a list variable 
configure_sublist <- function(., ss, vlist, df ) {
  lnames <- strsplit(names(df)[ss], '.', fixed=T )
  ss1    <- which(names(.[[vlist]])==lnames[[1]][2])
  ss2    <- which(names(.[[vlist]][[ss1]])==lnames[[1]][3])
  .[[vlist]][[ss1]][ss2] <- df[ss] 
  return(1) 
} 



### END ###    
