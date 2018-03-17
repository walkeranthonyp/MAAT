###############################
#
# Canopy object 
# 
# AWalker December 2015
#
################################

setwd('../leaf')
source('leaf_object.R')
setwd('../canopy')
source('canopy_functions.R')
source('canopy_system_functions.R')



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
    # main run function
    
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
        list(A=.$state$integrated$A, gs=.$state$integrated$gs, respiration=.$state$integrated$respiration)
        
      } else if(.$pars$output=='leaf') {
        list(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
             gi=.$state$integrated$gi, gs=.$state$integrated$gs, respiration=.$state$integrated$respiration, lim=NA)
        
      } else if(.$pars$output=='all_lim') {
        list(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
             gi=.$state$integrated$gi, gs=.$state$integrated$gs, respiration=.$state$integrated$respiration, lim=NA, 
             Acg_lim=.$state$integrated$Acg_lim, 
             Ajg_lim=.$state$integrated$Ajg_lim, 
             Apg_lim=.$state$integrated$Apg_lim, 
             layers_Acg_lim=.$state$integrated$layers_Acg_lim, 
             layers_Ajg_lim=.$state$integrated$layers_Ajg_lim, 
             layers_Apg_lim=.$state$integrated$layers_Apg_lim
        )
        
      } else if(.$pars$output=='full') {
        c(.$state$integrated, .$state_pars)
      }
    }    

    
    
    ###########################################################################
    # Variables
    
    # function names
    fnames <- list(
      cansys          = 'f_cansys_multilayer',
      can_pars        = 'f_canlight_pars',
      can_scale_light = 'f_canlight_beerslaw_goudriaan',
      can_scale_N     = 'f_leafN_CLMuniform',
      can_scale_Ca    = 'f_Ca_uniform',
      can_scale_vpd   = 'f_vpd_uniform',
      lai             = 'f_constant'
    )
    
    # parameters
    pars <- list(
      verbose    = F,
      output     = 'run',
      layers     = 10,
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
    
    # Environment
    env <- list(
      par      = numeric(1),      
      par_dir  = numeric(1),      
      par_diff = numeric(1),      
      ca_conc  = numeric(1),
      vpd      = numeric(1),
      zenith   = 0
    )
    
    # state parameters
    state_pars <- list(
      m            = numeric(1),    
      G_dir        = numeric(1),
      k_dir        = numeric(1),
      k_diff       = numeric(1),
      k_dirprime   = numeric(1),
      k_diffprime  = numeric(1),
      lscattering  = numeric(1),
      alb_dir      = numeric(1),
      alb_diff     = numeric(1),
      alb_dir_can  = numeric(1),
      alb_diff_can = numeric(1)
    )
    
    # state
    state <- list(
      # External
      lai     = numeric(1),      # 1.5 for Sphagnum Williams & Flannagan, 1998
      mass_a  = 10,
      C_to_N  = 40,
      totalN  = numeric(1),
      
      # Calculated state
      # canopy layer vectors
      vert    = data.frame(
        # variable canopy environment etc
        par_dir         = numeric(1),
        par_diff        = numeric(1),
        apar_sun        = numeric(1),
        apar_shade      = numeric(1),
        f_sun           = numeric(1),
        f_shade         = numeric(1),
        leaf.ca_conc    = numeric(1),
        leaf.vpd        = numeric(1),
        leaf.par        = numeric(1),
        leaf.leafN_area = numeric(1),
        # variable canopy physiology
        A               = numeric(1),
        respiration     = numeric(1),
        ci              = numeric(1),
        cc              = numeric(1),
        rb              = numeric(1),
        rs              = numeric(1),
        ri              = numeric(1),
        lim             = numeric(1)
      ),
      
      # integrated canopy values
      integrated = list(
        A              = numeric(1),        # canopy assimilation rate                         (umol m-2s-1)
        Acg_lim        = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
        Ajg_lim        = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
        Apg_lim        = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
        layers_Acg_lim = numeric(1),        # number of canopy layers Ac limited        
        layers_Ajg_lim = numeric(1),        # number of canopy layers Aj limited        
        layers_Apg_lim = numeric(1),        # number of canopy layers Ap limited
        cb             = numeric(1),        # canopy mean boundary layer CO2                   (Pa)
        ci             = numeric(1),        # canopy mean leaf internal CO2                    (Pa) 
        cc             = numeric(1),        # canopy mean chloroplast CO2                      (Pa)
        #gb             = numeric(1),        # canopy boundary conductance                      (mol m-2s-1)
        #gs             = numeric(1),        # canopy stomatal conductance                      (mol m-2s-1) 
        #gi             = numeric(1),        # canopy leaf internal conductance                 (mol m-2s-1)
        rb             = numeric(1),        # canopy boundary resistance                       (m2s mol-1)
        rs             = numeric(1),        # canopy stomatal resistance                       (m2s mol-1) 
        ri             = numeric(1),        # canopy leaf internal resistance                  (m2s mol-1)
        respiration    = numeric(1)         # canopy respiration rate                          (umol m-2s-1)        
      )
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
      .$configure(vlist='env', df=.$dataf$met[l,], F )
      
      # run model
      .$run()              
    }
    

    # initialise the number of layers in the canopy
    init_vert <- function(.,l) {
      .$state$vert <- lapply(.$state$vert, function(v,leng) numeric(leng), leng=l )
    }
    
    
    # function to run the leaves within the canopy
    run_leaf <- function(.,ii,df){
      # This wrapper function is called from an lapply or mclapply function to run over each leaf in the canopy
      # assumes that each row of the dataframe are independent and non-sequential
      
      .$leaf$configure(vlist='env',   df=df[ii,], F )
      .$leaf$configure(vlist='state', df=df[ii,], F )
      
      # run leaf
      .$leaf$run()        
    }
    
    
    
    #######################################################################           
    # Test functions
    
    .test <- function(.,verbose=T){
      
      # Child Objects
      .$leaf <- as.proto(leaf_object$as.list(),all.names=T)
      .$leaf$cpars$output <- 'all_lim'

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



### END ###
