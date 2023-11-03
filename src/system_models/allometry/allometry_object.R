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
  canopy_packing     = 'f_canopy_packing',
  height             = 'f_height_hdbh_obrien1995',
  leaf_biomass       = 'f_leaf_biomass_dbh_salda',
  crown_area         = 'f_crown_area_fates',
  crown_depth        = 'f_crown_depth_fates',
  abg_biomass        = 'f_abg_biomass_dbh_saldariaga1998',
  blg_biomass        = 'f_blg_biomass_dbh_fracabg',
  wood_biomass       = 'f_wood_biomass',
  crown_lai          = 'f_crown_lai_linear',
  fineroot_biomass   = 'f_fineroot_biomass_leaffrac'
  #coarseroot_biomass = 'f_coarseroot_biomass'
)


# environment
####################################
allometry_object$env <- list(
  dbh            = numeric(1),        # diameter at breast height (DBH, 1.3 m) (cm)      
  canopy_packing = numeric(1)         # metric of canopy packing (0-1), 0 is fully packed, 1 is fully open 
)


# state
####################################
allometry_object$state <- list(
  dbh              = numeric(1),      # diameter at breast height (DBH, 1.3 m) (cm) 
  height           = numeric(1),      # height (m) 
  leaf_biomass     = numeric(1),      # leaf biomass (kg C individual-1) 
  crown_area       = numeric(1),      # crown area (m2) 
  crown_depth      = numeric(1),      # crown depth (m) 
  crown_lai        = numeric(1),      # crown (m2 leaf m-2 crown/ground area) 
  abg_wood_biomass = numeric(1),      # above-ground wood biomass (kg C individual-1) 
  blg_wood_biomass = numeric(1),      # below-ground biomass (kg C individual-1) 
  wood_biomass     = numeric(1),      # wood biomass (kg C individual-1) 
  fineroot_biomass = numeric(1)       # fine-root biomass (kg C individual-1) 
)


# state parameters (i.e. calculated parameters)
####################################
allometry_object$state_pars <- list(
  leaf_area     = numeric(1),         # crown leaf area (m2)  
  crwnarea_norm = numeric(1)          # normalisation parameter for carown area to DBH allometry, influenced by canopy packing
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
                                   dbh=5 ) {

  if(verbose) { str(.); print(.$env) }
  .$build(mod_out='full', switches=c(diag,verbose,cverbose) )

  .$env$dbh <- dbh

  .$run()
}


allometry_object$.test1 <- function(., verbose=F,
                                  allometry.text='f_text_combine',
                                  allometry.calcval='f_calcval_product',
                                  allometry.print_out='f_print_out_textonly' ) {
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  
  .$fnames$text      <- allometry.text
  .$fnames$calcval   <- allometry.calcval
  .$fnames$print_out <- allometry.print_out

  # configure_test must be called if any variables in the fnames lista re reassigned
  .$configure_test()  
  .$run()
}



### END ###
