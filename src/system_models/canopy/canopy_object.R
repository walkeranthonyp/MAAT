###############################
#
# Canopy object 
# 
# AWalker December 2015
#
################################

library(proto)

source('canopy_functions.R')
source('canopy_system_functions.R')



# CANOPY OBJECT
###############################################################################

canopy_object <-
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects, & build function
    
    name <- 'canopy'
    
    # expected child objects
    # the 'leaf_object' object named 'leaf'
    child_list <- list('leaf') 
    leaf <- NULL
    
    # build function
    build <- function(., mod_mimic=NULL, ... ) {
    
      # read default model setup for highest level model
      source('../../functions/general_functions.R')
      init_default <- readXML(paste(.$name,'default.xml',sep='_'))
     
      # read model mimic setup
      if(!is.null(mod_mimic)&F) {
        setwd('mimic_xmls')
        print(paste(.$name,'mimic:', mod_mimic ))
        init_mimic   <- readXML(paste(.$name,'_',mod_mimic,'.xml',sep=''))
        init_default <- fuselists(init_default,init_mimic)
        setwd('..')
      }

      # assign default and mod mimic values to data structure
      .$configure(vlist='fnames', df=unlist(init_default$fnames))
      .$configure(vlist='pars',   df=unlist(init_default$pars))
      .$configure(vlist='env',    df=unlist(init_default$env))

      # build child objects
      setwd(paste0('../',.$child_list[[1]]))
      source(paste0(.$child_list[[1]],'_object.R'))
      .$leaf     <- as.proto( leaf_object$as.list() )
      rm(leaf_object, pos=1 )
      .$leaf$build(mod_mimic=mod_mimic)
      .$leaf$cpars$output <- .$name 
      setwd(paste0('../',.$name))
    }
    
    
    
    ###########################################################################
    # main run function
    
    run <- function(.) {
      if(.$cpars$verbose) print('canopy run()')
     
      # initialise canopy
      .$state$lai <- get(.$fnames$lai)(.) # this also could point to a higher level plant object  
      get(.$fnames$pars_init)(.)

      # assign canopy environment to leaf environment 
      # any canopy scaling of these variables will overwrite the values written here 
      envss     <- which(names(.$env) %in% names(.$leaf$env) )
      df        <- as.data.frame(.$env[envss])
      names(df) <- paste0('leaf.',names(df))
      .$leaf$configure(vlist='env', df=df ) 
  
      # set leaf absorptance to 1 as all cansys functions (should) account for leaf scattering
      .$leaf$pars$a <- 1.0

      # calculate water status
      get(.$fnames$water_status)(.)

      # calculate diffuse and direct radiation
      get(.$fnames$par_partition)(.)      
      
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
    
      if(.$cpars$output=='run') {
        c(A=.$state$integrated$A, gs=.$state$integrated$gs, rd=.$state$integrated$rd)
        
      } else if(.$cpars$output=='leaf') {
        c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
          gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd, lim=NA)
        
      } else if(.$cpars$output=='all_lim') {
        c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
          gi=.$state$integrated$gi, gs=.$state$integrated$gs, rd=.$state$integrated$rd, lim=NA, 
          Acg_lim=.$state$integrated$Acg_lim, 
          Ajg_lim=.$state$integrated$Ajg_lim, 
          Apg_lim=.$state$integrated$Apg_lim 
        )
        
      } else if(.$cpars$output=='full') {
        c(.$state$integrated, .$state_pars)

      } else if(.$cpars$output=='leaf_dem') {
        vapply(.$state$vert$layer, function(v) v, .$state$vert$layer[[1]] )

      }
    }    

    
    
    ###########################################################################
    # Variables
    
    # function names
    fnames <- list(
      cansys        = 'f_cansys_multilayer',
      pars_init     = 'f_pars_init',
      rt            = 'f_rt_beerslaw_goudriaan',
      scale_n       = 'f_scale_n_CLMuniform',
      scale_vcmax   = 'f_scale_vcmax_beerslaw',
      vcmax0        = 'f_vcmax0_constant',
      k_vcmax       = 'f_k_vcmax_constant',
      scale_ca      = 'f_scale_ca_uniform',
      scale_vpd     = 'f_scale_vpd_uniform',
      lai           = 'f_lai_constant',
      par_partition = 'f_par_partition_spitters',
      water_status  = 'f_water_status_none',
      fwdw          = 'f_fwdw_wtd_lin'
    )
    
    # parameters
    pars <- list(
      layers           = 10,
      lai              = 10,
      lai_max          = 4,
      lai_curve        = 0.5,
      leaf_cores       = 1,
      G                = 0.5,     # light extinction coefficient assuming leaves are black bodies and randomly distributed horizontally, 0.5 assumes random or spherical leaf orientation, 1.5 for Sphagnum Williams & Flannagan, 1998
      can_clump        = 1,       # canopy clumping coefficient, 1 - random horizontal distribution, leaves become more clumped as coefficient goes towards zero.
      k_layer          = 0.5,     # for multilayer canopy, where in the layer to calculate physiology, 0 - bottom, 0.5 - midway, 1 - top; not the correct solution to the simplifying assumption of Beer's law (Wang 2003) 
      alb_soil         = 0.15,    # soil albedo
      leaf_reflectance = 0.075,   # leaf reflectance
      vcmax0           = 35,      # vcmax at extreme top of canopy
      k_vcmax          = 0.2,     # scaling exponent for vcmax through canopy
      k_vcmax_expa     = -2.43,   # intercept parameter in exponnent to calculate scaling exponent for vcmax through canopy
      k_vcmax_expb     = 9.63e-3, # slope parameter in exponnent to calculate scaling exponent for vcmax through canopy
      fwdw_wl_slope    = -0.022,  # delta sphagnum fwdw ratio per mm of decrease in water level      (mm-1), currently from Adkinson & Humpfries 2010, Rydin 1985 has similar intercept but slope seems closer to -0.6 
      fwdw_wl_sat      = 16,      # sphagnum fwdw ratio at 0 water level, currently from Adkinson & Humpfries 2010     
      fwdw_wl_exp_a    = -0.037,  # decrease in sphagnum fwdw ratio as an exponential f of water level (cm), currently from Strack & Price 2009
      fwdw_wl_exp_b    = 3.254    # decrease in sphagnum fwdw ratio as an exponential f of water level (cm) 
    )
    
    # Environment
    env <- list(
      temp      = 25,      
      par       = 1000,      
      par_dir   = numeric(1),      
      par_diff  = numeric(1),      
      ca_conc   = 400,
      vpd       = 1,
      clearness = 1,
      zenith    = 0,
      water_td  = numeric(1),
      sphag_h   = numeric(1)
    )
    
    # state parameters
    state_pars <- list(
      m            = numeric(1),    
      G_dir        = numeric(1),
      k_dir        = numeric(1),
      k_diff       = numeric(1),
      k_dirprime   = numeric(1),
      k_diffprime  = numeric(1),
      k_vcmax      = numeric(1),
      lscattering  = numeric(1),
      alb_dir      = numeric(1),
      alb_diff     = numeric(1),
      alb_dir_can  = numeric(1),
      alb_diff_can = numeric(1)
    )
    
    # state
    state <- list(
      # External
      lai     = numeric(1), 
      mass_a  = 10,
      C_to_N  = 40,
      totalN  = 7,
      vcmax0  = numeric(1),      
 
      # Calculated state
      # canopy layer vectors
      vert    = list(
        # variable canopy environment etc
        leaf = list( 
          leaf.ca_conc     = numeric(1),
          leaf.vpd         = numeric(1),
          leaf.par         = numeric(1),
          leaf.atref.vcmax = numeric(1), 
          leaf.leafN_area  = numeric(1)
        ),
        # variable canopy light & physiology by sun and shade leaves
        sun = list( 
          apar     = numeric(1),
          fraction = numeric(1),
          A        = numeric(1),
          rd       = numeric(1),
          ci       = numeric(1),
          cc       = numeric(1),
          gb       = numeric(1),
          gs       = numeric(1),
          gi       = numeric(1),
          g        = numeric(1),
          lim      = numeric(1)
        ),
        shade = list( 
          apar     = numeric(1),
          fraction = numeric(1),
          A        = numeric(1),
          rd       = numeric(1),
          ci       = numeric(1),
          cc       = numeric(1),
          gb       = numeric(1),
          gs       = numeric(1),
          gi       = numeric(1),
          g        = numeric(1),
          lim      = numeric(1)
        ),
        layer = list( 
          apar     = numeric(1),
          A        = numeric(1),
          rd       = numeric(1),
          ci       = numeric(1),
          cc       = numeric(1),
          gb       = numeric(1),
          gs       = numeric(1),
          gi       = numeric(1),
          g        = numeric(1),
          Acg_lim  = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
          Ajg_lim  = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
          Apg_lim  = numeric(1)         # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
        )
      ),
      
      # integrated canopy values
      integrated = list(
        apar       = numeric(1),        # canopy absorbed PAR
        A          = numeric(1),        # canopy assimilation rate                         (umol m-2s-1)
        Acg_lim    = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
        Ajg_lim    = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
        Apg_lim    = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
        cb         = numeric(1),        # canopy mean boundary layer CO2                   (Pa)
        ci         = numeric(1),        # canopy mean leaf internal CO2                    (Pa) 
        cc         = numeric(1),        # canopy mean chloroplast CO2                      (Pa)
        gb         = numeric(1),        # canopy leaf boundary conductance                 (mol H2O m-2 s-1)
        gs         = numeric(1),        # canopy stomatal conductance                      (mol H2O m-2 s-1)
        gi         = numeric(1),        # canopy leaf internal conductance                 (mol CO2 m-2 s-1)
        g          = numeric(1),        # canopy total conductance                         (mol H2O m-2 s-1)
        rd         = numeric(1)         # canopy respiration rate                          (umol m-2s-1)        
      )
    )

    # run control parameters
    cpars <- list(
      verbose       = F,          # write diagnostic output during runtime 
      cverbose      = F,          # write configuration output during runtime 
      output        = 'run'       # type of output from run function
    )
    
      
    
    ###########################################################################
    # Run & configure functions
    
    configure <- function(., vlist, df, o=T ) {
      # This function is called from any of the run functions, or during model initialisation
      # - sets the values within .$fnames / .$pars / .$env / .$state to the values passed in df 

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
      if(.$cpars$cverbose&o) {
        print('', quote=F )
        print('Leaf configure:', quote=F )
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
      #print(paste('Leaf conf:', vlist, names(df), df ))
      if(length(slss)>0)    vapply( slmss, .$configure_sublist, numeric(1), vlist=vlist, df=df ) 
      if(length(nslmss)>0) .[[vlist]][vlss] <- df[nslmss]
      #print(paste(df[nslmss],.[[vlist]][vlss])) 

      # call child (leaf) assign 
      #print(paste('conf:',vlist, names(df), df, length(mss) ))
      if(any(listnames[1,]!=.$name)) {
        dfc <- if(length(mss)>0) df[-which(listnames[1,]==.$name)] else df 
        vapply( .$child_list, .$child_configure , 1, vlist=vlist, df=dfc )
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
    child_configure <- function(., child, vlist, df ) { 
      #print(paste('child conf:',vlist, names(df), df ))
      .[[child]]$configure(vlist=vlist, df=df ) ; return(1) 
    }
    
    
    run_met <- function(.,l){
      # This wrapper function is called from an lapply function to run this model over every row of a dataframe
      # assumes that each row of the dataframe are sequential
      # allows the system state at t-1 to affect the system state at time t if necessary (therefore mclapply cannot be used)
      # typically used to run the model using data collected at a specific site and to compare against observations
      
      # expects .$dataf$met to exist in the object, usually in a parent "wrapper" object
      # any "env" variables specified in the "drv$env" dataframe and specified here will be overwritten by the values specified here 
      
      # met data assignment
      .$configure(vlist='env', df=.$dataf$met[l,] )
      
      # run model
      .$run()              
    }
    

    # initialise the number of layers in the canopy
    init_vert <- function(.,l) {
      .$state$vert$leaf  <- lapply(.$state$vert$leaf,  function(v,leng) numeric(leng), leng=l )
      .$state$vert$sun   <- lapply(.$state$vert$sun,   function(v,leng) numeric(leng), leng=l )
      .$state$vert$shade <- lapply(.$state$vert$shade, function(v,leng) numeric(leng), leng=l )
      .$state$vert$layer <- lapply(.$state$vert$layer, function(v,leng) numeric(leng), leng=l )
    }
    
    
    # function to run the leaves within the canopy
    run_leaf <- function(.,ii,df){
      # This wrapper function is called from an (v/l)apply function to run over each leaf in the canopy
      # assumes that each row of the dataframe are independent and non-sequential
      
      .$leaf$configure(vlist='env',   df=df[ii,] )
      .$leaf$configure(vlist='state', df=df[ii,] )
      .$leaf$configure(vlist='pars',  df=df[ii,] )
      
      # run leaf
      .$leaf$run()        
    }
    
    
    
    #######################################################################           
    # Test functions
    
    .test <- function(., verbose=T, par=2000, ca_conc=400, lai=6 ) {
      
      # Child Objects
      #.$leaf <- as.proto(leaf_object$as.list(),all.names=T)
      #.$leaf$cpars$output <- 'all_lim'
      .$build()

      # parameter settings
      .$cpars$verbose       <- verbose
      .$leaf$cpars$verbose  <- F
      
      .$env$par        <- par
      .$env$ca_conc    <- ca_conc
      .$pars$lai       <- lai
      .$state$mass_a   <- 175
      .$state$C_to_N   <- 40
      
      .$run()
    }
    
    .test_aca <- function(., verbose=F, verbose_loop=F, canopy.par=c(100,1000), canopy.ca_conc=seq(50,1200,50),
                          rs = 'f_r_zero' ) {
      
      # Child Objects
      #.$leaf <- as.proto(leaf_object$as.list(),all.names=T)
      #.$leaf$cpars$output <- 'all_lim'
      .$build()
      .$leaf$fnames$rs    <- rs

      .$cpars$verbose       <- verbose
      .$leaf$cpars$verbose  <- F
      
      .$env$par        <- 2000
      .$env$ca_conc    <- 200
      .$pars$lai       <- 10
      .$state$mass_a   <- 175
      .$state$C_to_N   <- 40
      
      if(verbose) str.proto(canopy_object)
      
      .$dataf       <- list()
      .$dataf$met   <- expand.grid(mget(c('canopy.ca_conc','canopy.par')))
      
      .$dataf$out  <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      print(cbind(.$dataf$met,.$dataf$out))
      p1 <- xyplot(A~.$dataf$met$canopy.ca_conc|as.factor(.$dataf$met$canopy.par),.$dataf$out,abline=0,
                   ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'))
      print(p1)
    }
    
    #######################################################################           
    # End canopy object    
})



### END ###
