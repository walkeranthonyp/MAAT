################################
#
# Canopy object 
# 
# AWalker December 2015
#
################################

source('canopy_functions.R')
source('canopy_system_functions.R')
source('../leaf/leaf_object.R')



# CANOPY OBJECT
###############################################################################

canopy_object <-
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'canopy'
    
    # expected child objects
    # the 'leaf_object' object named 'leaf'
    child_list <- list('leaf') 
    leaf <- NULL
    
    # build function
    build <- function(.,model) {
      .$leaf <- as.proto( leaf_object$as.list() )
    }
    
    
    
    ###########################################################################
    # main run functions
    
    # Alternate run function    
    run <- function(.){
      if(.$pars$verbose) print('canopy_run')
      
      # initialise canopy
      .$state$lai <- get(.$fnames$lai)(.) # this also could point to a higher level plant object  
      get(.$fnames$can_pars)(.)
      
      # run canopy model
      get(.$fnames$cansys)(.)
      
      #output
      .$output()      
    }
    

    
    ###########################################################################
    # Output function
    
    # output processing function
    # -- returns a vector of outputs
    output <- function(.){
      if(.$pars$output=='run') {
        list(A=.$state$integrated$A,gs=.$state$integrated$gs,respiration=.$state$integrated$respiration)
        
      } else if(.$pars$output=='leaf') {
        list(A=.$state$integrated$A,cc=.$state$integrated$cc,ci=.$state$integrated$ci,
             gi=.$state$integrated$gi,gs=.$state$integrated$gs,respiration=.$state$integrated$respiration,lim=NA)
        
      } else if(.$pars$output=='all_lim') {
        list(A=.$state$integrated$A,cc=.$state$integrated$cc,ci=.$state$integrated$ci,
             gi=.$state$integrated$gi,gs=.$state$integrated$gs,respiration=.$state$integrated$respiration,lim=NA,
             wc_lim=.$state$integrated$wc_lim,
             wj_lim=.$state$integrated$wj_lim,
             wp_lim=.$state$integrated$wp_lim,
             layers_wc_lim=.$state$integrated$layers_wc_lim,
             layers_wj_lim=.$state$integrated$layers_wj_lim,
             layers_wp_lim=.$state$integrated$layers_wp_lim
        )
        
      } else if(.$pars$output=='full') {
        c(.$state$integrated,.$state_pars)
      }
    }    

    
    
    ###########################################################################
    # Variables
    
    # function names
    fnames <- list(
      cansys          = 'f_cansys_multilayer',
      can_pars        = 'f_canlight_pars',
      can_scale_light = 'f_canlight_beerslaw',
      can_scale_N     = 'f_leafN_CLMuniform',
      can_scale_Ca    = 'f_Ca_uniform',
      can_scale_vpd   = 'f_vpd_uniform',
      lai             = 'f_constant'
    )
    
    # Environment
    env <- list(
      par      = numeric(0),      
      par_dir  = numeric(0),      
      par_diff = numeric(0),      
      ca_conc  = numeric(0),
      vpd      = numeric(0),
      zenith   = 0
    )
    
    # state
    state <- list(
      # External
      lai     = numeric(0),      # 1.5 for Sphagnum Williams & Flannagan, 1998
      mass_a  = 10,
      C_to_N  = 40,
      totalN  = numeric(0),
      
      # Calculated state
      # canopy layer vectors
      # variable canopy environment etc
      vert    = data.frame(
        par_dir         = numeric(0),
        par_diff        = numeric(0),
        apar_sun        = numeric(0),
        apar_shade      = numeric(0),
        f_sun           = numeric(0),
        f_shade         = numeric(0),
        leaf.ca_conc    = numeric(0),
        leaf.vpd        = numeric(0),
        leaf.par        = numeric(0),
        leaf.leafN_area = numeric(0)
        ),
      # canopy state - logically these should be in the above vert dataframe 
      A           = numeric(0),
      respiration = numeric(0),
      ci          = numeric(0),
      cc          = numeric(0),
      gs          = numeric(0),
      gi          = numeric(0),
      lim         = character(0),
      
      # integrated canopy values
      integrated = list(
        A             = numeric(0),        # canopy assimilation rate                         (umol m-2s-1)
        wc_lim        = numeric(0),        # assimilation rate of canopy layers wc limited    (umol m-2s-1)
        wj_lim        = numeric(0),        # assimilation rate of canopy layers wj limited    (umol m-2s-1)        
        wp_lim        = numeric(0),        # assimilation rate of canopy layers wp limited    (umol m-2s-1)        
        layers_wc_lim = numeric(0),        # number of canopy layers wc limited        
        layers_wj_lim = numeric(0),        # number of canopy layers wj limited        
        layers_wp_lim = numeric(0),        # number of canopy layers wp limited
        cb            = numeric(0),        # canopy mean boundary layer CO2                   (Pa)
        ci            = numeric(0),        # canopy mean leaf internal CO2                    (Pa) 
        cc            = numeric(0),        # canopy mean chloroplast CO2                      (Pa)
        gb            = numeric(0),        # canopy boundary conductance                      (mol m-2s-1)
        gs            = numeric(0),        # canopy stomatal conductance                      (mol m-2s-1) 
        gi            = numeric(0),        # canopy leaf internal conductance                 (mol m-2s-1)
        rb            = numeric(0),        # canopy boundary resistance                       (m2s mol-1)
        rs            = numeric(0),        # canopy stomatal resistance                       (m2s mol-1) 
        ri            = numeric(0),        # canopy leaf internal resistance                  (m2s mol-1)
        respiration   = numeric(0)         # canopy respiration rate                          (umol m-2s-1)        
      )
    )
    
    # state parameters
    state_pars <- list(
      m            = numeric(0),    
      G_dir        = numeric(0),
      k_dir        = numeric(0),
      k_diff       = numeric(0),
      k_dirprime   = numeric(0),
      k_diffprime  = numeric(0),
      lscattering  = numeric(0),
      alb_dir      = numeric(0),
      alb_diff     = numeric(0),
      alb_dir_can  = numeric(0),
      alb_diff_can = numeric(0)
    )
    
    # parameters
    pars <- list(
      verbose    = F,
      output     = 'run',
      lai        = 10,
      lai_max    = 4,
      lai_curve  = 0.5,
      leaf_cores = 1,
      G          = 0.5,        # light extinction coefficient assuming leaves are black bodies and randomly distributed horizontally, 0.5 assumes random or spherical leaf orientation, 1.5 for Sphagnum Williams & Flannagan, 1998
      can_clump  = 1,          # canopy clumping coefficient, 1 - random horizontal distribution, leaves become more clumped as coefficient goes towards zero.
      k_layer    = 0,          # used by some to determine light scaling, not the correct solution to the simplifying assumption of Beer's law (Wang 2003) 
      alb_soil   = 0.15,       # soil albedo
      leaf_reflectance = 0.075 # leaf reflectance
    )
      
    
    
    ###########################################################################
    # Run & configure functions
    
    configure <- function(.,vlist,df,o=T) {
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
     
      # process UQ variables
      uqvars <- names(df)
      prefix <- vapply( strsplit(uqvars,'.', fixed=T), function(cv) cv[1], 'character' )
      modobj <- .$name
      dfss   <- which(prefix==modobj)
      vlss   <- match(uqvars[dfss], paste0(modobj,'.',names(.[[vlist]])) )

      # catch NAs in vlss
      if(any(is.na(vlss))) stop(paste('names mismatch between model object variables and input list variable:', uqvars[which(is.na(vlss))] ))

      # assign UQ variables
      .[[vlist]][vlss] <- df[dfss]

      # call child (leaf) configure
      child_configure <- function(., child ) if(any(prefix==child)) .[[child]]$configure(.,vlist,df,F) 
      vapply( .$child_list, .$child_configure , NULL )     
 
      if(.$cpars$cverbose&o) {
        print('',quote=F)
        print('Canopy configure:',quote=F)
        print(prefix,quote=F)
        print(df,quote=F)
        print(.[vlist],quote=F)
      }
    }
    
    
    run_met <- function(.,l){
      # This wrapper function is called from an lapply function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are sequential
      # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
      # typically used to run the model using data collected at a specific site and to compare against observations
      
      # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
      # any "env" variables specified in the "drv$env" dataframe and specified here will be overwritten by the values specified here 
      
      # met data assignment
      .$configure(func='write_env',df=.$dataf$met[l,],F)
      
      # run model
      .$run()              
    }
    

    # initialise the number of layers in the canopy
    init_vert <- function(.,l) {
      .$state$vert <- data.frame(
        par_dir         = numeric(l),
        par_diff        = numeric(l),
        apar_sun        = numeric(l),
        apar_shade      = numeric(l),
        f_sun           = numeric(l),
        f_shade         = numeric(l),
        leaf.ca_conc    = numeric(l),
        leaf.vpd        = numeric(l),
        leaf.par        = numeric(l),
        leaf.leafN_area = numeric(l)
      )
    }
    
    
    # function to run the leaves within the canopy
    run_leaf <- function(.,ii){
      # This wrapper function is called from an lapply or mclapply function to run over each leaf in the canopy
      # assumes that each row of the dataframe are independent and non-sequential
      
      # expects .$state$vert to exist in this object
      .$leaf$configure(func='write_env',  df=data.frame(.$state$vert[ii,]),F)
      .$leaf$configure(func='write_state',df=data.frame(.$state$vert[ii,]),F)
      
      # run leaf
      .$leaf$run()        
    }
    
    
    
    #######################################################################           
    # Test functions
    
    .test <- function(.,verbose=T){
      
      # Child Objects
      .$leaf <- as.proto(leaf_object$as.list(),all.names=T)

      # parameter settings
      .$pars$verbose       <- verbose
      .$leaf$pars$verbose  <- F
      .$pars$outfull       <- T
      
      .$env$par_dir    <- 2000
      .$env$ca_conc    <- 200
      .$pars$lai       <- 10
      .$state$mass_a   <- 175
      .$state$C_to_N   <- 40
      
      .$run()
    }
    
    .test_aca <- function(.,verbose=F,verbose_loop=F,canopy.par_dir=c(100,1000),canopy.ca_conc=seq(50,1200,50)){
      
      # Child Objects
      .$leaf <- as.proto(leaf_object$as.list(),all.names=T)

      .$pars$verbose      <- verbose
      .$leaf$pars$verbose  <- F
      .$pars$outfull       <- F
      
      .$env$par_dir    <- 2000
      .$env$ca_conc    <- 200
      .$pars$lai       <- 10
      .$state$mass_a   <- 175
      .$state$C_to_N   <- 40
      
      if(verbose) str.proto(canopy_object)
      
      .$dataf       <- list()
      .$dataf$met   <- expand.grid(mget(c('canopy.ca_conc','canopy.par_dir')))
      
      .$dataf$out  <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      print(cbind(.$dataf$met,.$dataf$out))
      p1 <- xyplot(A~.$dataf$met$canopy.ca_conc|as.factor(.$dataf$met$canopy.par_dir),.$dataf$out,abline=0,
                   ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'))
      print(p1)
    }
    
    #######################################################################           
    # End canopy object    
})



