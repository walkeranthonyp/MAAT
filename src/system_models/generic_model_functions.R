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
  # APW: could add a catch here to test that requested eval variables exist.
  # APW: will need to be done once latest soil model development has been merged
  # APW: to work with soil_decomp will need to be called after structure build

  # build model pool structure
  # APW: prob can delete n_pools case in latest version
  if(!is.null(.$pars$n_outfluxes)) {
    .$build_pool_structure(init_default$pars[[.$name]]$n_pools, init_default$pars[[.$name]]$n_outfluxes )
  } else if(!is.null(.$pars$n_pools)) 
    .$build_pool_structure(init_default$pars[[.$name]]$n_pools)

  # assign default and mod mimic values to data structure
  .$configure(vlist='pars',   df=unlist(init_default$pars))
  .$configure(vlist='env',    df=unlist(init_default$env))
  .$configure(vlist='fnames', df=unlist(init_default$fnames), init=T )

  # build child objects
  if(!is.null(.$child_list)) vapply(.$child_list, .$build_child, numeric(0), mod_mimic=mod_mimic )
}


# function to build lists that vary in length depending on pool size
# - takes the number of n_pools/n_outfluxes from init files (max n if in dynamic init)
# - and replaces fnames' decomp, outflux, transfer lists and lists of pool_pars with NA to that n
# APW: note this does not recreate the lists to n, just wipes them, could add that functionality
# - also creates state matrix and outflux list  
build_pool_structure <- function(., init_n_pools, init_n_outfluxes=NULL ) {

  # generate lists indexed by pool 
  # number of model pools
  n_pools <- 
    if(!is.null(.$.super[['init_dynamic']][['pars']][[.$name]][['n_pools']])) {
      max(.$init_dynamic$pars[[.$name]]$n_pools)
    } else if(!is.null(.$.super[['init_static']][['pars']][[.$name]][['n_pools']])) { 
      .$init_static$pars[[.$name]]$n_pools
    #} else .$pars$n_pools 
    } else init_n_pools 
 
  print('', quote=F )
  print(paste0('Building ', .$name, ' model pool structure with: ', n_pools, ' pools.'), quote=F )

  if(!is.null(.$fnames$decomp)) {
    print('',quote=F)
    print('Indexing pool fnames lists ...', quote=F ) 

    # generate fnames decomp list
    print('decomp')
    .$fnames$decomp <- list() 
    lnames <- paste0('d',1:n_pools) # may need a variable maxnpools
    .$fnames$decomp[lnames] <- NA  
    #print(.$fnames$decomp, quote=F )

    # generate transfer list
    if(!is.null(.$fnames$transfer)) {
      print('transfer')
      .$fnames$transfer <- list() 
      tn     <- expand.grid(1:n_pools,1:n_pools)
      # remove transfers from/to the same pool
      tn     <- tn[-((1:n_pools -1)*n_pools + 1:n_pools),]
      lnames <- apply(tn, 1, function(v) paste0('t',v[1],'_to_',v[2]) )
      .$fnames$transfer[lnames] <- NA  
      #print(.$fnames$transfer, quote=F )
    }

    print('',quote=F)
    print('Indexing pool_state_pars lists ...', quote=F ) 
    for(l in .$pool_state_pars) { 
      if(is.list(.$fnames[[l]])) {
        print(paste(l,'fname'))
        .$fnames[[l]] <- list()
        lnames        <- paste0(l,1:n_pools) 
        .$fnames[[l]][lnames] <- paste('f', l, 'constant', sep='_' )  
      }
      if(is.list(.$state_pars[[l]])) {
        print(paste(l,'state_par'))
        .$state_pars[[l]] <- list()
        lnames            <- paste0(l,1:n_pools) 
        .$state_pars[[l]][lnames] <- NA 
      }
    }

    # - for parameters that are indexed by the number of pools
    print('',quote=F)
    print('Indexing pool_pars lists ...',quote=F)  
    for(l in .$pool_pars) {
      if(is.list(.$pars[[l]])) {
        print(l)
        .$pars[[l]] <- list()
        lnames      <- paste0(l,1:n_pools) 
        .$pars[[l]][lnames] <- NA 
      }
    }
    # assign default values for input_coef & cstate0 lists
    .$pars$input_coef[] <- 0
    .$pars$cstate0[]    <- 1
  } 


  # generate lists indexed by outflux
  # number of model outfluxes
  if(!is.null(init_n_outfluxes)) {
    n_outfluxes <- 
      if(!is.null(.$.super[['init_dynamic']][['pars']][[.$name]][['n_outfluxes']])) {
        max(.$init_dynamic$pars[[.$name]]$n_outfluxes)
      } else if(!is.null(.$.super[['init_static']][['pars']][[.$name]][['n_outfluxes']])) { 
        .$init_static$pars[[.$name]]$n_outfluxes
      #} else .$pars$n_outfluxes 
      } else init_n_outfluxes 
   
    print('', quote=F )
    print(paste0('Building ', .$name, ' model outflux structure with: ', n_outfluxes, ' fluxes.'), quote=F )
  
    if(!is.null(.$fnames$outflux)) {
      # - for fnames that are indexed by the number of outfluxes
      print('',quote=F)
      print('Indexing outflux_fnames lists ...', quote=F ) 
      print('outflux')
      .$fnames$outflux <- list() 
      lnames <- paste0('of',1:n_outfluxes) 
      .$fnames$outflux[lnames] <- NA  
      #print('  outflux list:', quote=F )
      #print(.$fnames$outfluxes, quote=F )
  
      # - for fnames that are indexed by the number of outfluxes
      for(l in c('tcor', 'wcor', 'scor')) { 
        if(is.list(.$fnames[[l]])) {
          print(l)
          .$fnames[[l]] <- list()
          lnames        <- paste0(l,1:n_outfluxes) 
          .$fnames[[l]][lnames] <- paste('f', l, 'none', sep='_' ) 
        }
      }
      for(l in .$outflux_state_pars) { 
        if(is.list(.$fnames[[l]])) {
          print(l)
          .$fnames[[l]] <- list()
          lnames        <- paste0(l,1:n_outfluxes) 
          .$fnames[[l]][lnames] <- paste('f', l, 'constant', sep='_' )  
        }
      }
  
      # - for state_pars that are indexed by the number of outfluxes
      print('',quote=F)
      print('Indexing outflux_state_pars lists ...', quote=F ) 
      for(l in .$outflux_state_pars) { 
        if(is.list(.$state_pars[[l]])) {
          print(l)
          .$state_pars[[l]] <- list()
          lnames            <- paste0(l,1:n_outfluxes) 
          .$state_pars[[l]][lnames] <- NA 
        }
      }

      # - for parameters that are indexed by the number of outfluxes
      print('',quote=F)
      print('Indexing outflux_pars lists ...',quote=F)  
      for(l in .$outflux_pars) { 
        if(is.list(.$pars[[l]])) {
          print(l)
          .$pars[[l]] <- list()
          lnames      <- paste0(l,1:n_outfluxes) 
          .$pars[[l]][lnames] <- NA 
        }
      }
    } 
  }


  # generate state matrix and outflux list
  # APW: this seems redundant with init function in model object 
  .$state$cpools  <- matrix(1.01, ncol=1, nrow=n_pools )
  .$state$outflux <- as.list(numeric(.$pars$n_outflux))
  names(.$state$outflux) <- paste0('of',1:.$pars$n_outflux)
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



# main run functions
###########################################################################

run <- function(.) {

  print('')
  print(.$env) 
  print(.$state) 

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
# - maybe a result of being called from the wrapper and maybe . represents the object within which the function is called rather than to which it belongs
run_met <- function(.,l) {

  #print('')
  #print(.$fnames$decomp)
  #lapply( grep('decomp\\.',names(.$fns),value=T), function(char) {print(char); print(.$fns[[char]])} )
  #print(.$fnames$transfer)
  #lapply( grep('transfer\\.',names(.$fns),value=T), function(char) {print(char); print(.$fns[[char]])} )
  #print('')

  # initialize: pool structure etc
  if(!is.null(.$init)) .$init()
  
  # call steady state system model if it exists and doesn't return a null value
  # APW Matt: this is the additional function call to initialise the model at steady state
  #           if you want to turn this off just set the fnames$steadystate value to f_steadystate_null (a dummy function that just returns NULL)
  # APW : need to add this as an option to run_MAAT, if !init_steady then fnames$steady <- null
  if(!is.null(.$fns$steadystate)) if(!is.null(.$fns$steadystate() )) .$fns$steadystate()  

  # run over an input/meteorological dataset 
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
  listnames <- vapply( strsplit(names(df),'.', fixed=T), 
                 function(cv) {
                   #print(cv)
                   cv3<-character(3); cv3[1:length(cv)]<-cv; t(cv3)
                   },
                 character(3) )

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
      fnslist <- as.list(rapply(.$fnames, function(c) if(is.na(c)) NA else get(c, pos=1 ) ))
      #fnslist <- as.list(rapply(.$fnames, function(c) if(is.na(c)) NA else {print(c); get(c, pos=1 )} ))
      .$fns   <- as.proto(fnslist, parent=. )

      # specific methods assignment (for methods not included in fnames)
      if(!is.null(.$configure_unique)) .$configure_unique(init=T, flist=unlist(.$fnames) )

    # assign methods only updated in fnames to methods list
    } else if(length(mss)>0) {
      flist <- unlist(.$fnames[vlss_full])
      for(n in 1:length(flist)) {
        .$fns[[names(flist[n])]]              <- if(is.na(flist[n])) NA else get(flist[n], pos=1 )
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
  #print(vlist)
  #print(ss)
  #print(names(df)[ss])
  lnames <- strsplit(names(df)[ss], '.', fixed=T )
  ss1    <- which(names(.[[vlist]])==lnames[[1]][2])
  #print(names(.[[vlist]]))
  #print(lnames)
  #print(ss1)
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
  listnames <- vapply( strsplit(names(df),'.',fixed=T), function(cv) {cv3<-character(3); cv3[1:length(cv)]<-cv; t(cv3)}, character(3) )

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
  fnslist <- as.list(rapply(.$fnames, function(c) if(is.na(c)) NA else get(c, pos=1 ) ))
  .$fns   <- as.proto(fnslist, parent=. )
  if(!is.null(.$configure_unique)) .$configure_unique(init=T, flist=unlist(.$fnames) )
  if(!is.null(.$init))             .$init()
}



### END ###
