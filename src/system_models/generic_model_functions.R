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

# build function that initialises the object and calls build functions of all child objects
build <- function(., mod_mimic=NULL, mod_out='run', child=F, switches=c(diag=F,verbose=F,cverbose=F), ... ) {

  # read default model setup for highest level model
  source('../../functions/general_functions.R')
  init_default <- readXML(paste(.$name,'default.xml',sep='_'))

  # set switches
  .$cpars$diag     <- switches[1] & !child
  .$cpars$verbose  <- switches[2]
  .$cpars$cverbose <- switches[3]
  .$cpars$mod_out  <- mod_out

  # read model mimic setup
  if(!is.null(mod_mimic)) {
    if(any(grepl(mod_mimic,list.files('./mimic_xmls')))) {
      setwd('mimic_xmls')
      print(paste(.$name,'mimic:',mod_mimic), quote=F )
      init_mimic   <- readXML(paste(.$name,'_',mod_mimic,'.xml',sep=''))
      init_default <- fuselists(init_default,init_mimic)
      setwd('..')
    } else print(paste('mimic:',mod_mimic,', for:',.$name,'not found.'), quote=F )
  }

  # print init_default if requested
  if(.$cpars$verbose) {
    print('',quote=F)
    print('Build, init_default',quote=F)
    print(init_default,quote=F)
  }

  # assign model output function
  if(.$cpars$diag) mod_out <- 'full' ### this will assign full to all child objects too could add a child switch
  .$output <- get(paste('f', 'output', .$name, .$cpars$mod_out, sep='_' ))
  # could add a catch here to test that requested eval variables exist.
  # will need to be done once latest soil model development has been merged
  # to work with soil_decomp will need to be called after structure build

  # assign default and mod mimic values to data structure
  .$configure(vlist='pars',   df=unlist(init_default$pars))
  .$configure(vlist='env',    df=unlist(init_default$env))
  .$configure(vlist='fnames', df=unlist(init_default$fnames), init=T )

  # build child objects
  if(!is.null(.$child_list)) vapply(.$child_list, .$build_child, numeric(0), mod_mimic=mod_mimic )
}


# function that calls child object build function
build_child <- function(., child_obj, mod_mimic=NULL, ... ) {

  # load child object
  child_obj_name  <- paste0(child_obj, '_object' )
  setwd(paste0('../',child_obj))
  source(paste0(child_obj_name, '.R' ))

  # build child object into parent object
  .[[child_obj]] <- as.proto( get(child_obj_name)$as.list(), parent=. ) # should work but may need environment setting
  rm(list=child_obj_name, pos=1 )
  .[[child_obj]]$build(mod_mimic=mod_mimic, mod_out=.$name, child=T )
  setwd(paste0('../',.$name))

  # return nothing
  numeric(0)
}


# function that reads object and all child objects configuration
read_config <- function(.) {

  # create config list
  co1 <- as.list(.)[c('fnames', 'pars', 'env' )]
  l1  <- list()
  
  # invert first two layers of list 
  co2 <- lapply(co1, 
                function(l) {
                  l2 <- list( l )
                  names(l2) <- .$name
                  l1 <- c(l1, l2 )
                })

  # call child objects
  if(!is.null(.$child_list)) {
    
    # co3 is a list of co2 type lists
    co3 <- lapply(.$child_list, 
                  function(c) {
                    .[[c]]$read_config()
                  })

    # fuse individually with co2
    for(i in 1:length(co3))  {
      for(n in names(co2)) co2[[n]] <- c(co2[[n]], co3[[i]][[n]] )
    rm(co3)
  }}

  co2
}




# main run functions
###########################################################################

run <- function(.) {

  # call system model
  .$fns$sys()

  # print to screen
  if(.$cpars$verbose) print(.$state)

  # output
  .$output()
}



# init and run met functions
# wrapper function called from an lapply function to run model over every row of a meteorology dataframe
###########################################################################

init_state <- function(.) {
  .$state      <- rapply(.$state, function(v) numeric(length(v)), how='replace' )
  .$state_pars <- rapply(.$state_pars, function(v) numeric(length(v)), how='replace' )
}


# this currently works both when called from unit testing and from the wrapper 
# - not 100 % sure why as when called from the wrapper .$dataf should read .super$dataf
# - maybe a result of being called from the wrapper and maybe . represents the object within which the function is called?
run_met <- function(.,l) {

  if(!is.null(.$init)) .$init()
  t(vapply(1:.$dataf$lm, .$run_met1, .$dataf$mout )) 
}


run_met1 <- function(.,l) {
  # assumes that each row of the dataframe are sequential
  # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
  # typically used to run the model using data collected at a specific site and to compare against observations

  # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
  # any "env" variables specified in the "dataf$env" dataframe but also specified in .$dataf$met will be overwritten by the .$dataf$met values

  # met data assignment
  .$configure_met(df=.$dataf$met[,l])

  # run model
  .$run()
}



## output functions
############################################################################

f_output_eval <- function(.){

  lout <- unlist(c(.$state, .$state_pars ))

  # NAs in this match operation means requested eval variables have been mis-labled - needs a trap 
  lout[match(names(.$dataf$out), names(lout) )]
}


# configure functions
###########################################################################

# check character names on row vectors of dataf functions
configure_check <- function(., vlist='env', df ) {
  if(!is.character(names(df))) {
    print('', quote=F )
    print(paste('Check dataf names:',.$name,',',vlist,'.'), quote=F )
    print('dft:', quote=F )
    print(df, quote=F )
    print('names(df):', quote=F )
    print(names(df), quote=F )
    stop('FATAL ERROR: names(df) is not a character vector, will cause strsplit to fail.')
  }
}


# This function is called from any of the run functions, or during model initialisation
# - sets the values within .$fnames / .$pars / .$env / .$state to the values passed in df
configure <- function(., vlist, df, init=F, o=T ) {

  # error catch
  # APW: move this to the wrapper object - it only needs tested once
  if(!is.character(names(df))) {
    print('', quote=F )
    print(paste('Configure:',.$name,',',vlist,'.'), quote=F )
    print('df passed to configure:', quote=F )
    print(df, quote=F )
    print('names(df):', quote=F )
    print(names(df), quote=F )
    stop('FATAL ERROR: names(df) is not a character vector, will cause strsplit to fail.')
  }

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
  vlss_full <- vlss
  slss      <- which(listnames[3,mss]!='')
  if(length(slss)>0) {
    slmss  <- mss[slss]
    nslmss <- mss[-slss]
    vlss   <- vlss[-slss]
  } else {
    slmss  <- NULL
    nslmss <- mss
  }

  # print configure setup if requested
  if(.$cpars$cverbose&o) {
    print('', quote=F )
    print(paste('Configure:',.$name,',',vlist,'.'), quote=F )
    print('df passed to configure:', quote=F )
    print(df, quote=F )
    print('listnames:', quote=F )
    print(listnames, quote=F )
    print(paste('  subscripts for this model object:',.$name), quote=F )
    print(mss, quote=F )
    print('  subscripts for variables that are lists:', quote=F )
    print(slmss, quote=F )
    print('  subscripts for variables that are not lists:', quote=F )
    print(nslmss, quote=F )
    print(paste('  subscripts in:',.$name, vlist ,' for variables that are not lists:'), quote=F )
    print(vlss, quote=F )
    print(names(.[[vlist]]), quote=F )
  }

  # assign UQ variables
  #print(paste(.$name,'configure:', vlist, names(df), df ))
  if(length(slss)>0)    vapply( slmss, .$configure_sublist, numeric(1), vlist=vlist, df=df )
  if(length(nslmss)>0) .[[vlist]][vlss] <- df[nslmss]

  # assign methods to methods list
  if(vlist=='fnames') {

    # assign all methods to methods list
    if(init) {
      fnslist <- as.list(rapply(.$fnames, function(c) get(c, pos=1 ) ))
      .$fns   <- as.proto(fnslist, parent=. )

      # specific methods assignment (for methods not included in fnames)
      if(!is.null(.$configure_unique)) .$configure_unique(init=T, flist=unlist(.$fnames) )

    # assign methods only updatred in fnames to methods list
    } else if(length(mss)>0) {
      flist <- unlist(.$fnames[vlss_full])
      for(n in 1:length(flist)) {
        .$fns[[names(flist[n])]]              <- get(flist[n], pos=1 )
        environment(.$fns[[names(flist[n])]]) <- .$fns
      }

      # specific methods assignment (for methods not included in fnames)
      if(!is.null(.$configure_unique)) .$configure_unique(flist=flist)
    }
  }

  # call child configure
  #print(paste('conf:',vlist, names(df), df, length(mss) ))
  if(!is.null(.$child_list) & any(listnames[1,]!=.$name)) {
    dfc <- if(length(mss)>0) df[-which(listnames[1,]==.$name)] else df
    vapply( .$child_list, .$configure_child , numeric(0), vlist=vlist, df=dfc )
  }
}


# configure a list variable
configure_sublist <- function(., ss, vlist, df ) {
  lnames <- strsplit(names(df)[ss], '.', fixed=T )
  ss1    <- which(names(.[[vlist]])==lnames[[1]][2])
  ss2    <- which(names(.[[vlist]][[ss1]])==lnames[[1]][3])
  .[[vlist]][[ss1]][ss2] <- df[ss]
  return(1)
}


# call a child configure function
configure_child <- function(., child, vlist, df ) {

  .[[child]]$configure(vlist=vlist, df=df )
  numeric(0)
}


# configure function for meteorological / boundary conditions
configure_met <- function(., df ) {

  #print('configure met')
  #print(df)
  #print(class(names(df)))

  # split variable names at .
  listnames <- vapply( strsplit(names(df),'.', fixed=T), function(cv) {cv3<-character(3); cv3[1:length(cv)]<-cv; t(cv3)}, character(3) )

  # df subscripts for model object
  mss  <- 1:length(df)

  # variable list subscripts in model object data structure
  vlss <- match(listnames[2,], names(.[['env']]) )

  # remove NAs in vlss from vlss and mss
  if(any(is.na(vlss))) {
    mss  <- mss[-which(is.na(vlss))]
    vlss <- vlss[-which(is.na(vlss))]
  }

  # assign UQ variables
  .[['env']][vlss] <- df[mss]
}


# configure function for .test functions to configure fns from reassigned fnames
configure_test <- function(.) {

  # configure methods
  fnslist <- as.list(rapply(.$fnames, function(c) get(c, pos=1 ) ))
  .$fns   <- as.proto(fnslist, parent=. )
  if(!is.null(.$configure_unique)) .$configure_unique(init=T, flist=unlist(.$fnames) )
  if(!is.null(.$init))             .$init()
}



### END ###
