###############################
#
# Canopy structure object 
# 
# AWalker October 2018
#
################################

library(proto)

source('canopy_structure_functions.R')
source('canopy_structure_system_functions.R')



# CANOPY OBJECT
###############################################################################

canopy_structure_object <-
  proto(expr={
    
    ###########################################################################
    # Object name, expected child objects, & build function
    
    name <- 'canopy_structure'
    
    # expected child objects
    # the 'canopy_object' object named 'canopy'
    child_list <- list('canopy') 
    canopy     <- NULL
    
    # build function
    build <- function(., mod_mimic=NULL, ... ) {
    
      # read default model setup for highest level model
      source('../../functions/general_functions.R')
      init_default <- readXML(paste(.$name,'default.xml',sep='_'))
     
      # read model mimic setup
      if(!is.null(mod_mimic)&F) {
        setwd('mimic_xmls')
        print(paste('Canopy structure mimic:', mod_mimic ))
        init_mimic   <- readXML(paste(.$name,'_',mod_mimic,'.xml',sep=''))
        init_default <- fuselists(init_default,init_mimic)
        setwd('..')
      }

      # build child objects
      setwd('../canopy')
      source('canopy_object.R')
      .$canopy     <- as.proto( canopy_object$as.list() )
      rm(canopy_object, pos=1 )
      init_child <- .$canopy$build(mod_mimic=mod_mimic)
      .$canopy$cpars$output <- 'all_lim'
      setwd(paste0('../',.$name))

      # build full init list
      c(init_default, init_child )
    }
    
    
    
    ###########################################################################
    # main run function
    
    run <- function(.) {
      if(.$cpars$verbose) print('canopy_structure run()')
    
      # assign canopy_structure environment to canopy environment 
      # any canopy_structure scaling of these variables will overwrite the values written here 
      envss     <- which(names(.$env) %in% names(.$canopy$env) )
      df        <- as.data.frame(.$env[envss])
      names(df) <- paste0('canopy.',names(df))
      .$canopy$configure(vlist='env', df=df ) 

      # run canopy_structure model
      get(.$fnames$canstrctsys)(.)
      
      #output
      .$output()      
    }
    

    
    ###########################################################################
    # Output function
    
    # output processing function
    # -- returns a vector of outputs
    output <- function(.){
      if(.$cpars$output=='run') {
        c(A=.$state$integrated$A, gs=.$state$integrated$gs, respiration=.$state$integrated$respiration)
        
      } else if(.$cpars$output=='canopy') {
        c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
          gi=.$state$integrated$gi, gs=.$state$integrated$gs, respiration=.$state$integrated$respiration, lim=NA)
        
      } else if(.$cpars$output=='all_lim') {
        c(A=.$state$integrated$A, cc=.$state$integrated$cc, ci=.$state$integrated$ci, 
          gi=.$state$integrated$gi, gs=.$state$integrated$gs, respiration=.$state$integrated$respiration, lim=NA 
          #Acg_lim=.$state$integrated$Acg_lim, 
          #Ajg_lim=.$state$integrated$Ajg_lim, 
          #Apg_lim=.$state$integrated$Apg_lim, 
          #layers_Acg_lim=.$state$integrated$layers_Acg_lim, 
          #layers_Ajg_lim=.$state$integrated$layers_Ajg_lim, 
          #layers_Apg_lim=.$state$integrated$layers_Apg_lim
        )
        
      } else if(.$cpars$output=='full') {
        c(unlist(.$state$integrated), unlist(.$state_pars) )
      }
    }    

    
    
    ###########################################################################
    # Variables
    
    # function names
    fnames <- list(
      canstrctsys    = 'f_canstrctsys_leafdem',
      pars_init      = 'f_pars_init',
      lai            = 'f_lai_leafdem_env',
      leafdem_upper  = 'f_leafdem_upper_wu2017',
      leafdem_lower  = 'f_leafdem_lower_wu2017',
      leafdem_vcmax0 = 'f_leafdem_vcmax0_constant',
      water_status   = 'f_water_status_none',
      fwdw           = 'f_fwdw_wtd_lin'
    )
    
    # parameters
    pars <- list(
      canopy_cores     = 1,
      layers           = 30,
      lai              = 10,
      lai_upper        = 2.5,
      lai_ftop         = 0.7,
      vcmax0 = list(
        young  = 7,
        mature = 35,
        old    = 20
      ),
      lai_max          = 4,
      lai_curve        = 0.5,
      fwdw_wl_slope    = -0.022, # delta sphagnum fwdw ratio per mm of decrease in water level      (mm-1), currently from Adkinson & Humpfries 2010, Rydin 1985 has similar intercept but slope seems closer to -0.6 
      fwdw_wl_sat      = 16,     # sphagnum fwdw ratio at 0 water level, currently from Adkinson & Humpfries 2010     
      fwdw_wl_exp_a    = -0.037, # decrease in sphagnum fwdw ratio as an exponential f of water level (cm), currently from Strack & Price 2009
      fwdw_wl_exp_b    = 3.254   # decrease in sphagnum fwdw ratio as an exponential f of water level (cm) 
    )
    
    # Environment
    env <- list(
      temp       = 25,      
      par        = 1000,      
      ca_conc    = 400,
      vpd        = 1,
      zenith     = 0,
      lai_young  = 1, 
      lai_mature = 2,
      lai_old    = 3,
      water_td   = numeric(1),
      sphag_h    = numeric(1)
    )
    
    # state parameters
    state_pars <- list(
      alb_dir_can  = numeric(1),
      alb_diff_can = numeric(1)
    )
    
    # state
    state <- list(
      # External
      lai            = numeric(1),    
      lai_young      = numeric(1),    
      lai_mature     = numeric(1),    
      lai_old        = numeric(1),    
      upper_can_prop = numeric(3), 
      lower_can_prop = numeric(3), 
      
      #mass_a  = 10,
      #C_to_N  = 40,
      #totalN  = 7,
     
      # leaf demography associated traits
      leafdem_traits <- list(
        vcmax0 = numeric(3)
      ),
 
      # Calculated state
      # canopy_structure layer vectors
      vert    = list(
        # canopy_structure leaf demography 
        young = list( 
          apar    = numeric(1),        # canopy absorbed PAR
          A       = numeric(1),
          rd      = numeric(1),
          Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
          Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
          Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
          cb      = numeric(1),
          ci      = numeric(1),
          cc      = numeric(1),
          g       = numeric(1),  
          gb      = numeric(1),
          gs      = numeric(1),
          gi      = numeric(1)
        ),
        mature = list( 
          apar    = numeric(1),        # canopy absorbed PAR
          A       = numeric(1),
          rd      = numeric(1),
          Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
          Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
          Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
          cb      = numeric(1),
          ci      = numeric(1),
          cc      = numeric(1),
          g       = numeric(1),  
          gb      = numeric(1),
          gs      = numeric(1),
          gi      = numeric(1)
        ),
        old = list( 
          apar    = numeric(1),        # canopy absorbed PAR
          A       = numeric(1),
          rd      = numeric(1),
          Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
          Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
          Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
          cb      = numeric(1),
          ci      = numeric(1),
          cc      = numeric(1),
          g       = numeric(1),  
          gb      = numeric(1),
          gs      = numeric(1),
          gi      = numeric(1)
        ),
        layer = list( 
          apar    = numeric(1),        # canopy absorbed PAR
          A       = numeric(1),
          rd      = numeric(1),
          Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
          Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
          Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
          cb      = numeric(1),
          ci      = numeric(1),
          cc      = numeric(1),
          g       = numeric(1),  
          gb      = numeric(1),
          gs      = numeric(1),
          gi      = numeric(1)
        )
      ),
      
      # integrated canopy values
      integrated = list(
        apar    = numeric(1),        # canopy absorbed PAR
        A       = numeric(1),        # canopy assimilation rate                         (umol m-2s-1)
        rd      = numeric(1),        # canopy respiration rate                          (umol m-2s-1)        
        Acg_lim = numeric(1),        # assimilation rate of canopy layers Ac limited    (umol m-2s-1)
        Ajg_lim = numeric(1),        # assimilation rate of canopy layers Aj limited    (umol m-2s-1)        
        Apg_lim = numeric(1),        # assimilation rate of canopy layers Ap limited    (umol m-2s-1)        
        cb      = numeric(1),        # canopy mean boundary layer CO2                   (Pa)
        ci      = numeric(1),        # canopy mean canopy internal CO2                  (Pa) 
        cc      = numeric(1),        # canopy mean chloroplast CO2                      (Pa)
        g       = numeric(1),        # canopy conductance                               (mol H2O m-2 s-1 )
        gb      = numeric(1),        # canopy leaf boundary conductance                 (mol H2O m-2 s-1 )
        gs      = numeric(1),        # canopy stomatal conductance                      (mol H2O m-2 s-1 )
        gi      = numeric(1)         # canopy leaf internal conductance                 (mol CO2 m-2 s-1 ) 
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
    

    # initialise the number of layers in the canopy_structure
    init_vert <- function(.,l) {
      .$state$vert$young  <- lapply(.$state$vert$young,  function(v,leng) numeric(leng), leng=l )
      .$state$vert$mature <- lapply(.$state$vert$mature, function(v,leng) numeric(leng), leng=l )
      .$state$vert$old    <- lapply(.$state$vert$old,    function(v,leng) numeric(leng), leng=l )
      .$state$vert$layer  <- lapply(.$state$vert$layer,  function(v,leng) numeric(leng), leng=l )
    }
    
    
    # function to run the canopy within the canopy_structure
    run_canopy <- function(.,ii,df){
      # This wrapper function is called from an (v/l)apply function to run over each canopy in the canopy_structure
      # assumes that each row of the dataframe are independent and non-sequential
      
      #.$canopy$configure(vlist='env',   df=df[ii,] )
      #.$canopy$configure(vlist='state', df=df[ii,] )
      .$canopy$configure(vlist='pars', df=df[ii,] )
      
      # run canopy
      .$canopy$run()        
    }
    
    
    
    #######################################################################           
    # Test functions
    
    .test <- function(.,verbose=T){
      
      # Child Objects
      .$build()

      # parameter settings
      .$cpars$verbose       <- verbose
      .$canopy$cpars$verbose  <- F
      
      .$env$par        <- 2000
      .$env$ca_conc    <- 400
      .$env$lai        <- 10
      #.$state$mass_a   <- 175
      #.$state$C_to_N   <- 40
      
      .$run()
    }
    
    .test_aca <- function(., verbose=F, verbose_loop=F, 
                          canopy_structure.par=c(100,1000), 
                          canopy_structure.ca_conc=seq(50,1200,50),
                          rs = 'f_r_zero' ) {
      
      # Child Objects
      .$build()
      .$canopy$leaf$fnames$rs <- rs

      .$cpars$verbose         <- verbose
      .$canopy$cpars$verbose  <- F
      
      #.$env$par        <- 2000
     # .$env$ca_conc    <- 200
     # .$pars$lai       <- 10
     # .$state$mass_a   <- 175
     # .$state$C_to_N   <- 40
      
      if(verbose) str.proto(canopy_structure_object)
     
      # generate met data 
      .$dataf       <- list()
      .$dataf$met   <- expand.grid(mget(c('canopy_structure.ca_conc','canopy_structure.par')))
      
      # run maat
      .$dataf$out  <- data.frame(do.call(rbind,lapply(1:length(.$dataf$met[,1]),.$run_met)))
      
      # run output
      print(cbind(.$dataf$met,.$dataf$out))
      p1 <- 
        if(length(canopy_structure.ca_conc) > 2) 
          xyplot(A~.$dataf$met$canopy_structure.ca_conc|as.factor(.$dataf$met$canopy_structure.par),
                 .$dataf$out,abline=0,
                 ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),
                 xlab=expression(C[a]*' ['*mu*mol*' '*mol^-1*']'))
        else if(length(canopy_structure.par) > 2) 
          xyplot(A~.$dataf$met$canopy_structure.par|as.factor(.$dataf$met$canopy_structure.ca_conc),
                 .$dataf$out,abline=0,
                 ylab=expression('A ['*mu*mol*' '*m^-2*s-1*']'),
                 xlab=expression('PAR ['*mu*mol*' '*m^-2*s^-1*']'))
      
      print(p1)
    }
    
    #######################################################################           
    # End canopy_structure object    
})



### END ###
