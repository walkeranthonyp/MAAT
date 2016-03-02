################################
#
# Canopy object 
# 
# AWalker December 2015
#
################################

source('leaf_object.R')
source('canopy_functions.R')


# CANOPY OBJECT
###############################################################################

canopy_object <-
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects & build function
    
    name <- 'canopy'
    
    # expected child objects
    # the 'leaf_object' object named 'leaf'
    leaf <- NULL
    
    # build function
    .build <- function(.) {
      model      <- as.proto(.$as.list())
      model$leaf <- as.proto(leaf_object$as.list(),parent=model)
      model
    }
    
    
    
    ###########################################################################
    # main run functions
    
    run <- function(.){
      if(.$pars$verbose) print('canopy_run')
      
      # initialise canopy
      .$state$lai <- get(.$fnames$lai)(.) # this also could point to a higher level plant object  
      layers      <- ceiling(.$state$lai) # this could be a function specifying either no. of layers or the below lai assignment   
      .$init_vert(l=layers)
      get(.$fnames$can_pars)(.)
      
      # canopy layer 
      # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
      # populate layer dataframe      
      .$state$vert$leaf.ca_conc    <- get(.$fnames$can_scale_Ca)(.,l=1:layers,layers=layers)
      .$state$vert$leaf.par        <- get(.$fnames$can_scale_light)(.,l=1:layers,layers=layers)
      .$state$vert$leaf.leafN_area <- get(.$fnames$can_scale_N)(.,l=1:layers,layers=layers)
      
      # run leaf
      lout   <- as.data.frame(do.call(rbind,lapply(1:layers,.$run_leaf)))
#       lout   <- mclapply(1:layers,.$leaf$wrapper,df=leafdf,mc.cores=.$pars$leaf_cores) # this can be used to run each canopy layer in parallel, in practise this was found to be slower
      
      # assign data
      .$state$A           <- unlist(lout$A)
      .$state$cc          <- unlist(lout$cc)
      .$state$ci          <- unlist(lout$ci)
#       .$state$gi          <- unlist(lout$gi)
#       .$state$gs          <- unlist(lout$gs)
      .$state$ri          <- unlist(lout$ri)
      .$state$rs          <- unlist(lout$rs)
      .$state$lim         <- unlist(lout$lim)
      .$state$respiration <- unlist(lout$respiration)      
      
      #scale bottom layer fluxes/pars to leaf area in that layer
      # - need to generalise this to a multilayer canopy
      if(ceiling(.$state$lai)-floor(.$state$lai)==1){
        scale <- .$state$lai %% 1
        .$state$A[layers]           <- .$state$A[layers]           * scale        
#         .$state$gi[layers]          <- .$state$gi[layers]          * scale        
#         .$state$gs[layers]          <- .$state$gs[layers]          * scale        
        .$state$ri[layers]          <- .$state$ri[layers]          * scale        
        .$state$rs[layers]          <- .$state$rs[layers]          * scale        
        .$state$respiration[layers] <- .$state$respiration[layers] * scale      
        .$state$cc[layers]          <- .$state$cc[layers]          * scale        
        .$state$ci[layers]          <- .$state$ci[layers]          * scale        
      }      
      
      #integrate canopy layers
      # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
      # canopy sum values
      .$state$integrated$A             <- sum(.$state$A)
      .$state$integrated$respiration   <- sum(.$state$respiration)
      .$state$integrated$wc_lim        <- sum(.$state$A * (.$state$lim=='wc')) 
      .$state$integrated$wj_lim        <- sum(.$state$A * (.$state$lim=='wj'))
      .$state$integrated$wp_lim        <- sum(.$state$A * (.$state$lim=='wp'))
      .$state$integrated$layers_wc_lim <- sum(.$state$lim=='wc')
      .$state$integrated$layers_wj_lim <- sum(.$state$lim=='wj')
      .$state$integrated$layers_wp_lim <- sum(.$state$lim=='wp')
#       .$state$integrated$gi            <- sum(.$state$gi)
#       .$state$integrated$gs            <- sum(.$state$gs)
      .$state$integrated$ri            <- 1 / sum(1/.$state$ri)
      .$state$integrated$rs            <- 1 / sum(1/.$state$rs)
      # canopy mean values
      .$state$integrated$cc            <- sum(.$state$cc) / .$state$lai
      .$state$integrated$ci            <- sum(.$state$ci) / .$state$lai
      
      #output
      .$output()      
    }

    
    # Alternate run function    
    run_alt   <- function(.){
      
      # initialise canopy
      .$state$lai <- get(.$fnames$lai)(.) # this also could point to a higher level plant object  
      get(.$fnames$can_pars)(.)
      
      # run canopy model
      get(.$fnames$canopy)(.)
      
#       #integrate canopy layers
#       # - need to generalise this to allow direct light only, both direrct and diffuse (and sunlit and shade leaves)
#       # canopy sum values
#       .$state$integrated$A             <- sum(.$state$A)
#       .$state$integrated$respiration   <- sum(.$state$respiration)
#       .$state$integrated$Ac_lim        <- sum(.$state$A * (.$state$lim=='ac')) 
#       .$state$integrated$Aj_lim        <- sum(.$state$A * (.$state$lim=='aj'))
#       .$state$integrated$Ap_lim        <- sum(.$state$A * (.$state$lim=='ap'))
#       .$state$integrated$layers_ac_lim <- sum(.$state$lim=='ac')
#       .$state$integrated$layers_aj_lim <- sum(.$state$lim=='aj')
#       .$state$integrated$layers_ap_lim <- sum(.$state$lim=='ap')
#       .$state$integrated$gi            <- sum(.$state$gi)
#       .$state$integrated$gs            <- sum(.$state$gs)
#       # canopy mean values
#       .$state$integrated$cc            <- sum(.$state$cc) / .$state$lai
#       .$state$integrated$ci            <- sum(.$state$ci) / .$state$lai
      
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
      can_pars        = 'f_canlight_pars',
      canopy          = 'f_multilayer',
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
    
    # initialisation function
    init <- function(.,init_ls) {
      # expects to have the wrapper object named 'maat' as parent
      
      comb_init_list <- function(.,lls,cls) {
        nas <- sum(is.na(lls))
        if(nas>1) stop else if(nas==0) names(lls) <- paste('leaf',names(lls),sep='.') # need to write error message in here
        comb_ls <- if(nas==0) lls else NA

        nas <- sum(is.na(cls))
        if(nas>1) stop else if(nas==0) names(cls) <- paste('canopy',names(cls),sep='.') # need to write error message in here
        if(nas==0&is.na(comb_ls)) cls else c(comb_ls,cls) 
        
#         if(!is.na(lls)) names(lls) <- paste('leaf',names(lls),sep='.')
#         if(!is.na(cls)) names(cls) <- paste('canopy',names(cls),sep='.')
#         comb_ls <- if(!is.na(lls)) lls else NA
#         if(!is.na(cls)&is.na(comb_ls)) cls else c(comb_ls,cls) 
      }
      
      maat$static$fnames <- comb_init_list(lls=.$init_ls$lfs,cls=.$init_ls$cfs)
      maat$static$pars   <- comb_init_list(lls=.$init_ls$lps,cls=.$init_ls$cps)
      maat$static$env    <- comb_init_list(lls=.$init_ls$les,cls=.$init_ls$ces)

      maat$vars$fnames   <- comb_init_list(lls=.$init_ls$lfv,cls=.$init_ls$cfv)
      maat$vars$pars     <- comb_init_list(lls=.$init_ls$lpv,cls=.$init_ls$cpv)
      maat$vars$env      <- comb_init_list(lls=.$init_ls$lev,cls=.$init_ls$cev)
      
#       # Static during runtime
#       # add prefix to list names
#       if(!is.na(initls$lfs)) names(initls$lfs) <- paste('leaf',names(initls$lfs),sep='.')
#       if(!is.na(initls$cfs)) names(initls$cfs) <- paste('canopy',names(initls$cfs),sep='.')
#       if(!is.na(initls$lps)) names(initls$lps) <- paste('leaf',names(initls$lps),sep='.')
#       if(!is.na(initls$cps)) names(initls$cps) <- paste('canopy',names(initls$cps),sep='.')
#       if(!is.na(initls$les)) names(initls$les) <- paste('leaf',names(initls$les),sep='.')
#       if(!is.na(initls$ces)) names(initls$ces) <- paste('canopy',names(initls$ces),sep='.')
#       
#       # assign lists
#       fs <- if(!is.na(initls$lfs)) initls$lfs else NA
#       ps <- if(!is.na(initls$lps)) initls$lps else NA
#       es <- if(!is.na(initls$les)) initls$les else NA
#       fs <- if(!is.na(initls$cfs)&is.na(fs)) initls$cfs else c(fs,initls$cfs) 
#       ps <- if(!is.na(initls$cps)&is.na(ps)) initls$cps else c(ps,initls$cps) 
#       es <- if(!is.na(initls$ces)&is.na(es)) initls$ces else c(es,initls$ces)
#       .$static$fnames <- fs
#       .$static$pars   <- ps
#       .$static$env    <- es
#       
#       # Dynamic during runtime
#       # add prefix to list names
#       if(!is.na(initls$lfv)) names(initls$lfv) <- paste('leaf',names(initls$lfv),sep='.')
#       if(!is.na(initls$cfv)) names(initls$cfv) <- paste('canopy',names(initls$cfv),sep='.')
#       if(!is.na(initls$lpv)) names(initls$lpv) <- paste('leaf',names(initls$lpv),sep='.')
#       if(!is.na(initls$cpv)) names(initls$cpv) <- paste('canopy',names(initls$cpv),sep='.')
#       if(!is.na(initls$lev)) names(initls$lev) <- paste('leaf',names(initls$lev),sep='.')
#       if(!is.na(initls$cev)) names(initls$cev) <- paste('canopy',names(initls$cev),sep='.')
#       
#       # assign lists
#       fv <- if(!is.na(initls$lfv)) initls$lfv else NA
#       pv <- if(!is.na(initls$lpv)) initls$lpv else NA
#       ev <- if(!is.na(initls$lev)) initls$lev else NA
#       fv <- if(!is.na(initls$cfv)&is.na(fv)) initls$cfv else c(fv,initls$cfv) 
#       pv <- if(!is.na(initls$cpv)&is.na(pv)) initls$cpv else c(pv,initls$cpv) 
#       ev <- if(!is.na(initls$cev)&is.na(ev)) initls$cev else c(ev,initls$cev)
#       .$vars$fnames <- fv
#       .$vars$pars   <- pv
#       .$vars$env    <- ev
      
    }
    
    configure <- function(.,func,df,o=T) {
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames, .$pars, .$env, .$state to the values passed in df 
      # - df is a single row dataframe
      
      # name and assign the UQ variables
      uqvars     <- names(df)
      prefix     <- substr(uqvars,1,str_locate(uqvars,'\\.')[,2]-1)
      lapply(uqvars[which(prefix=='canopy')], func, .=., df=df)
      lapply(uqvars[which(prefix=='leaf')],   get(func,envir=.$leaf), .=.$leaf, df=df)
      
      if(.$pars$verbose&o) {
        print('',quote=F)
        print('Canopy configure:',quote=F)
        print(prefix,quote=F)
        print(df,quote=F)
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
    
    # variable assignment functions - called from the above configuration functions
    write_fnames <- function(.,var,df) .$fnames[which(names(.$fnames)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))] <- df[which(names(df)==var)]
    write_pars   <- function(.,var,df) .$pars[which(names(.$pars)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))]     <- df[which(names(df)==var)]
    write_env    <- function(.,var,df) .$env[which(names(.$env)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))]       <- df[which(names(df)==var)]
    write_state  <- function(.,var,df) .$state[which(names(.$state)==substr(var,str_locate(var,'\\.')[,2]+1,nchar(var)))]   <- df[which(names(df)==var)]

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



