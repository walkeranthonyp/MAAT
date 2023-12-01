################################
#
# MAAT allometry system representation functions (SRFs) 
# 
# AWalker November 2023 
#
################################



################################
# allometry system function one
# - this is a variant on the hello world program

f_sys_1 <- function(.) {

  # initialise
  .super$state$height              <- 0.0 
  .super$state$leaf_biomass        <- 0.0 
  .super$state$crown_area          <- 0.0 
  .super$state$crown_depth         <- 0.0 
  .super$state$crown_lai           <- 0.0
  .super$state$abg_wood_biomass    <- 0.0 
  .super$state$blg_wood_biomass    <- 0.0
  .super$state$wood_biomass        <- 0.0
  .super$state$fineroot_biomass    <- 0.0 
  .super$state$dbh                 <- .super$env$dbh
  .super$state_pars$dbh_effective  <- min(.super$state$dbh,.super$pars$dbh_max)
  .super$state_pars$crwnarea_norm  <- .$crwnallom() 

  # define state 
  .super$state$height               <- .$height()
  .super$state$leaf_biomass         <- .$leaf_biomass()
  .super$state$crown_area           <- .$crown_area()
  .super$state$crown_radius         <- .$crown_radius()
  .super$state$crown_depth          <- .$crown_depth()
  .super$state_pars$leaf_biomass_perarea <- 1e3 * .super$state$leaf_biomass / .super$state$crown_area
  .super$state$crown_lai            <- .$lai()
  .super$state_pars$crown_leaf_area <- .super$state$crown_area * .super$state$crown_lai 
  .super$state$abg_wood_biomass     <- .$abg_biomass()
  .super$state$blg_wood_biomass     <- .$blg_biomass()
  .super$state$wood_biomass         <- .$wood_biomass()
  .super$state$fineroot_biomass     <- .$fineroot_biomass()
#  .super$state$coarseroot_biomass <- .$coarseroot_biomass() --- synonymous with belowground wood in FATES


  # call print function
  #.$print_out()  
}



### END ###
