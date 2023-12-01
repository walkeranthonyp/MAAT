################################
#
# MAAT allometry model object 
#
# AWalker  November 2023 
#
################################

library(proto)
source('allometry_functions.R')
source('allometry_system_functions.R')



# allometry OBJECT
###############################################################################

# use generic allometry
setwd('..')
source('generic_model_object.R')
allometry_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('allometry')



# assign object functions
###########################################################################
allometry_object$name <- 'allometry'
# if this new object will contain nested child objects uncomment and edit the below lines of code and delete this comment 
# - otherwise delete all 5 lines
#allometry_object$child_list      <- list('child_name1') 
#allometry_object$build_child     <- build_child  
#allometry_object$configure_child <- configure_child  



# assign unique run function
###########################################################################
# if run function needs to be modified - add new function here



# functions unique to object that do not live in fnames/fns, i.e. do not vary ever
###########################################################################
# add structural functions (i.e not alternative process functions)  unique to model object here



# assign object variables 
###########################################################################

# function names
####################################
allometry_object$fnames <- list(
  sys                = 'f_sys_1',
  crwnallom          = 'f_crwnallom_canopy_mod',
  height             = 'f_height_hdbh_obrien1995',
  leaf_biomass       = 'f_leaf_biomass_dbh_salda',
  crown_area         = 'f_crown_area_fates',
  crown_radius       = 'f_crown_radius',
  crown_depth        = 'f_crown_depth_fates',
  abg_biomass        = 'f_abg_biomass_dbh_saldariaga1998',
  blg_biomass        = 'f_blg_biomass_dbh_fracabg',
  wood_biomass       = 'f_wood_biomass',
  lai                = 'f_lai_kovenock_exponential_sla',
  lai_linear         = 'f_lai_linear_constant_sla',
  lai_exponential    = 'f_lai_exponential_sla',
  fineroot_biomass   = 'f_fineroot_biomass_leaffrac'
  #coarseroot_biomass = 'f_coarseroot_biomass'
)


# environment
####################################
allometry_object$env <- list(
  dbh              = numeric(1),        # diameter at breast height (DBH, 1.3 m) (cm)      
  canopy_openness  = numeric(1),        # metric of canopy openness (0-1; 0 is fully closed, 1 is fully open) 
  canopy_lai_above = 0                  # mean LAI of canopy above this individual (m2 m-2) 
)


# state
####################################
allometry_object$state <- list(
  dbh                     = numeric(1),      # diameter at breast height (DBH, 1.3 m) (cm) 
  height                  = numeric(1),      # height (m) 
  leaf_biomass            = numeric(1),      # leaf biomass (kg C individual-1) 
  crown_area              = numeric(1),      # crown area (m2) 
  crown_radius            = numeric(1),      # crown radius (m) 
  crown_depth             = numeric(1),      # crown depth (m) 
  crown_lai               = numeric(1),      # crown (m2 leaf m-2 crown/ground area) 
  abg_wood_biomass        = numeric(1),      # above-ground wood biomass (kg C individual-1) 
  blg_wood_biomass        = numeric(1),      # below-ground biomass (kg C individual-1) 
  wood_biomass            = numeric(1),      # wood biomass (kg C individual-1) 
  fineroot_biomass        = numeric(1)       # fine-root biomass (kg C individual-1) 
)


# state parameters (i.e. calculated parameters)
####################################
allometry_object$state_pars <- list(
  crown_leaf_area              = numeric(1),         # crown leaf area (m2)  
  leaf_biomass_perarea         = numeric(1),         # leaf biomass per unit ground area (g C m-2, note this is in grams ) 
  leaf_biomass_perarea_slamax  = numeric(1),         # crown leaf biomass per unit ground area where slamax reached (kg C m-2)  
  crwnarea_norm                = numeric(1),         # normalisation parameter for carown area to DBH allometry, influenced by canopy packing
  dbh_effective                = numeric(1)          # effective DBH for allom functions that saturate  
)


# parameters
####################################
allometry_object$pars   <- list(
  
  # height
  dbh_max      = 90,                  # DBH where max height reached (cm)         
  hdbh_p1      = 0.64,
  hdbh_p2      = 0.37,
  hdbh_p3      = NA,
        
  # biomass
  wood_density = 0.7,                 # wood density (g cm-3)  
  frac_abg     = 0.6,                 # fraction woody biomass aboveground 
  fineroot_to_leaf_ratio = 0.6,       # ratio of fine-root to leaf biomass  
  abg_p1       = 0.06896, 
  abg_p2       = 0.572, 
  abg_p3       = 1.94, 
  abg_p4       = 0.931, 
  carbon_to_biomass = 2,              # conversion factor for carbon units to biomass dry weight (DW) 

  # leaf, crown area & depth
  slatop            = 0.0120,         # specific leaf area (SLA) at top of canopy (m2 g-1 C)
  slamax            = 0.0954,         # maximum specific leaf area (SLA) (m2 g-1 C)
  k                 = 0.2,            # extinction coefficient for SLA decline through crown (proportion) 
  crwnarea_p2       = 1.56,           # crown area to DBH exponent
  crwndepth_frac    = 0.5,            # crown depth fraction of height 
  crwnarea_dbh_max  = 0.6568464,      # crown depth fraction of height 
  crwnarea_dbh_min  = 0.3381119,      # crown depth fraction of height 
  ldbh_p1           = 0.07,
  ldbh_p2           = 1.3,
  ldbh_p3           = 0.55

)


# run control parameters
####################################
allometry_object$cpars <- list(
  verbose  = F,          # write diagnostic output during runtime 
  cverbose = F,          # write diagnostic output from configure function 
  output   = 'run'       # type of output from run function
)



# output functions
#######################################################################        

f_output_allometry_eval <- f_output_eval 

f_output_allometry_run <- function(.) {
  unlist(.$state)
}

f_output_allometry_state <- function(.) {
  unlist(.$state)
}

f_output_allometry_full <- function(.) {
  c(unlist(.$state),unlist(.$state_pars))
}



# test functions
#######################################################################        

allometry_object$.test <- function(., verbose=F, diag=F, cverbose=F, 
                                   dbh=5, allometry.lai=NULL ) {

  if(verbose) { str(.); print(.$env) }
  .$build(mod_out='full', switches=c(diag,verbose,cverbose) )

  if(!is.null(allometry.lai)) .$fnames$lai <- allometry.lai 
  
  .$configure_test()

  .$env$dbh <- dbh

  .$run()
}


allometry_object$.test_pars <- function(., verbose=F, diag=F, cverbose=F, 
                                   dbh=5,  
                                   dbh_max      = 90,                  # DBH where max height reached (cm)         
                                   hdbh_p1      = 0.64,
                                   hdbh_p2      = 0.37,
                                   hdbh_p3      = NA,
                                   # biomass
                                   wood_density = 0.7,           # wood density (g cm-3)  
                                   frac_abg     = 0.6,           # fraction woody biomass aboveground 
                                   fineroot_to_leaf_ratio = 0.6, # ratio of fine-root to leaf biomass  
                                   abg_p1       = 0.06896, 
                                   abg_p2       = 0.572, 
                                   abg_p3       = 1.94, 
                                   abg_p4       = 0.931, 
                                   carbon_to_biomass = 2,       # conversion factor for carbon units to biomass dry weight (DW) 
                                   # leaf, crown area & depth
                                   slatop       = 0.01995827,     # specific leaf area (SLA) at top of canopy (m2 g-1)
                                   crwnarea_p1  = 1.56,           # crown area to DBH exponent
                                   crwndepth_frac  = 0.5,         # crown depth fraction of height 
                                   crwnarea_dbh_max  = 0.6568464, # crown depth fraction of height 
                                   crwnarea_dbh_min  = 0.3381119, # crown depth fraction of height 
                                   ldbh_p1      = 0.07,
                                   ldbh_p2      = 1.3,
                                   ldbh_p3      = 0.55
                                   ) {

  if(verbose) { str(.); print(.$env) }
  .$build(mod_out='full', switches=c(diag,verbose,cverbose) )

  .$env$dbh <- dbh
  .$pars$dbh_max      <- dbh_max                  # DBH where max height reached (cm)         
  .$pars$hdbh_p1      <- hdbh_p1
  .$pars$hdbh_p2      <- hdbh_p2
  .$pars$hdbh_p3      <- hdbh_p3
  .$pars$wood_density <- wood_density           # wood density (g cm-3)  
  .$pars$frac_abg     <- frac_abg           # fraction woody biomass aboveground 
  .$pars$fineroot_to_leaf_ratio <- fineroot_to_leaf_ratio  # ratio of fine-root to leaf biomass  
  .$pars$abg_p1       <- abg_p1 
  .$pars$abg_p2       <- abg_p2 
  .$pars$abg_p3       <- abg_p3 
  .$pars$abg_p4       <- abg_p4 
  .$pars$slatop       <- slatop     # specific leaf area (SLA) at top of canopy (m2 g-1)
  .$pars$crwnarea_p1  <- crwnarea_p1           # crown area to DBH exponent
  .$pars$crwndepth_frac    <- crwndepth_frac         # crown depth fraction of height 
  .$pars$crwnarea_dbh_max  <- crwnarea_dbh_max  # crown depth fraction of height 
  .$pars$crwnarea_dbh_min  <- crwnarea_dbh_min  # crown depth fraction of height 
  .$pars$ldbh_p1      <- ldbh_p1
  .$pars$ldbh_p2      <- ldbh_p2
  .$pars$ldbh_p3      <- ldbh_p3

  .$run()
}


allometry_object$.test_dbh <- function(., verbose=F, diag=F, cverbose=F, 
                                   allometry.dbh=c(0.5,1:9,seq(10,200,5)),
                                   fnames.abg_biomass   = NULL,
                                   fnames.leaf_biomass  = NULL,
                                   fnames.height        = NULL,
                                   dbh_max      = 90,            # DBH where max height reached (cm)         
                                   hdbh_p1      = 0.64,
                                   hdbh_p2      = 0.37,
                                   hdbh_p3      = NA,
                                   # biomass
                                   wood_density = 0.7,           # wood density (g cm-3)  
                                   frac_abg     = 0.6,           # fraction woody biomass aboveground 
                                   fineroot_to_leaf_ratio = 0.6, # ratio of fine-root to leaf biomass  
                                   abg_p1       = 0.06896, 
                                   abg_p2       = 0.572, 
                                   abg_p3       = 1.94, 
                                   abg_p4       = 0.931, 
                                   carbon_to_biomass = 2,         # conversion factor for carbon units to biomass dry weight (DW) 
                                   # leaf, crown area & depth
                                   slatop       = 0.012,          # specific leaf area (SLA) at top of canopy (m2 g-1)
                                   crwnarea_p1  = 1.56,           # crown area to DBH exponent
                                   crwndepth_frac  = 0.5,         # crown depth fraction of height 
                                   crwnarea_dbh_max  = 0.6568464, # crown depth fraction of height 
                                   crwnarea_dbh_min  = 0.3381119, # crown depth fraction of height 
                                   ldbh_p1      = 0.07,
                                   ldbh_p2      = 1.3,
                                   ldbh_p3      = 0.55
                                   ) {

  if(verbose) { str(.); print(.$env) }
  .$build(mod_out='full', switches=c(diag,verbose,cverbose) )

  # configure fnames

  if(!is.null(fnames.abg_biomass))  .$fnames$abg_biomass = fnames.abg_biomass 
  if(!is.null(fnames.leaf_biomass)) .$fnames$leaf_biomass = fnames.leaf_biomass 
  if(!is.null(fnames.height))       .$fnames$height = fnames.height 
  .$configure_test() 

  .$pars$dbh_max      <- dbh_max                  # DBH where max height reached (cm)         
  .$pars$hdbh_p1      <- hdbh_p1
  .$pars$hdbh_p2      <- hdbh_p2
  .$pars$hdbh_p3      <- hdbh_p3
  .$pars$wood_density <- wood_density           # wood density (g cm-3)  
  .$pars$frac_abg     <- frac_abg           # fraction woody biomass aboveground 
  .$pars$fineroot_to_leaf_ratio <- fineroot_to_leaf_ratio  # ratio of fine-root to leaf biomass  
  .$pars$abg_p1       <- abg_p1 
  .$pars$abg_p2       <- abg_p2 
  .$pars$abg_p3       <- abg_p3 
  .$pars$abg_p4       <- abg_p4 
  .$pars$slatop       <- slatop     # specific leaf area (SLA) at top of canopy (m2 g-1)
  .$pars$crwnarea_p1  <- crwnarea_p1           # crown area to DBH exponent
  .$pars$crwndepth_frac    <- crwndepth_frac         # crown depth fraction of height 
  .$pars$crwnarea_dbh_max  <- crwnarea_dbh_max  # crown depth fraction of height 
  .$pars$crwnarea_dbh_min  <- crwnarea_dbh_min  # crown depth fraction of height 
  .$pars$ldbh_p1      <- ldbh_p1
  .$pars$ldbh_p2      <- ldbh_p2
  .$pars$ldbh_p3      <- ldbh_p3

  # configure met data and run
  .$dataf          <- list()
  .$dataf$met      <- t(as.matrix(expand.grid(mget(c('allometry.dbh')))))      
  .$dataf$lm       <- dim(.$dataf$met)[2]
  .$dataf$mout     <- .$output()
  .$dataf$out      <- .$run_met() 
  .$dataf$out_full <- as.data.frame(cbind(t(.$dataf$met), .$dataf$out ))
  print(head(.$dataf$out_full))
  .$dataf$out_full$BA <- pi * (allometry.dbh/2)^2
  
#  p1 <- xyplot(crown_area + crown_leaf_area ~ allometry.dbh , .$dataf$out_full, 
#	       type='l', abline=list(h=c(0,1)), auto.key=list(points=F, lines=T, x=0, y=1),  
#               ylab=expression('CA, LA ['*m^2*']'), xlab=expression(DBH*' ['*m*']') )
  p1 <- xyplot(crown_area ~ allometry.dbh , .$dataf$out_full, 
	       type='l', abline=list(h=c(0)), auto.key=list(points=F, lines=T, x=0, y=1),  
               ylab=expression('CA ['*m^2*']'), xlab=expression(DBH*' ['*cm*']') )
  p2 <- xyplot(height + crown_depth + I(10*crown_lai) + crown_radius ~ allometry.dbh , .$dataf$out_full, 
	       type='l', abline=list(h=c(0)), auto.key=list(points=F, lines=T, x=0, y=1),  
               ylab=expression('H, CD [m]; 10*CLAI ['*m^2*m^-2*']'), xlab=expression(DBH*' ['*cm*']') )
  p3 <- xyplot(abg_wood_biomass + blg_wood_biomass + wood_biomass + BA ~ allometry.dbh , .$dataf$out_full, 
	       type='l', abline=list(h=c(0)), auto.key=list(points=F, lines=T, x=0, y=1),  
               ylab=expression('Biomass ['*kg *' C]; BA ['*cm^2*']'), xlab=expression(DBH*' ['*cm*']') )
  p4 <- xyplot(leaf_biomass + fineroot_biomass~ allometry.dbh , .$dataf$out_full, 
	       type='l', abline=list(h=c(0)), auto.key=list(points=F, lines=T, x=0, y=1),  
               ylab=expression('Biomass ['*kg *' C]'), xlab=expression(DBH*' ['*cm*']') )
  print(p2, split=c(1,1,2,2), more=T )
  print(p1, split=c(2,1,2,2), more=T )
  print(p3, split=c(1,2,2,2), more=T )
  print(p4, split=c(2,2,2,2), more=F )
}



### END ###
