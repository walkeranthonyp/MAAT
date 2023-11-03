################################
#
# MAAT allometry process representation functions (PRFs)
# 
# AWalker  November 2023
#
################################

# info:
# https://fates-docs.readthedocs.io/en/latest/fates_tech_note.html#cohort-state-variables



### FUNCTIONS
################################

# - fates default lmode: 1, BCI hmode: 5



# height ~ dbh
#########################

# Obrien et al 1995
# FATES case 1
f_height_hdbh_obrien1995 <- function(.) {
  10^(log10(min(.super$state$dbh,.super$pars$dbh_max)) * (.super$pars$hdbh_p1+.super$pars$hdbh_p2)) 
}

# Poorter et al 2006
# FATES case 2 
f_height_hdbh_poorter2006 <- function(.) {
  .super$pars$hdbh_p1 * 
    (1 - exp(.super$pars$hdbh_p2*min(.super$state$dbh,.super$pars$dbh_max)^.super$pars$hdbh_p3) )  
}

# Power Law
# FATES case 3
f_height_hdbh_power <- function(.) {
  .super$pars$hdbh_p1 * min(.super$state$dbh,.super$pars$dbh_max)^.super$pars$hdbh_p2 
}

# Chave et al 2014 
# FATES case 4 
f_height_hdbh_chave2014 <- function(.) {
  #p1e = p1   ! -eclim (assumed that p1 already has eclim removed)
  dbh_eff <- min(.super$state$dbh,.super$pars$dbh_max)
  exp(.super$pars$hdbh_p1 + .super$pars$hdbh_p2*log(dbh_eff) + .super$pars$hdbh_p3*log(dbh_eff)^2 ) 
}

# Martinez Cano et al 2019
# FATES case 5 
f_height_hdbh_martinezcano2019 <- function(.) {
  (.super$pars$hdbh_p1 * min(.super$state$dbh,.super$pars$dbh_max)^.super$pars$hdbh_p2) / 
    (.super$pars$hdbh_p3 + .super$state$dbh^.super$state$hdbh_p2)  
}



# leaf biomass ~ dbh 
#########################

# - leaf biomass in kg C, dbh in cm

# Saldarriaga et al 1998 ???? -- doesn't seem to be
# FATES case 1
# FATES d2blmax_salda, biogeochem/FatesAllometryMod.F90
f_leaf_biomass_dbh_salda <- function(.) {
  .super$pars$ldbh_p1 * min(.super$state$dbh,.super$pars$dbh_max)^.super$pars$ldbh_p2 * 
    .super$pars$wood_density^.super$pars$ldbh_p3
}

# power 
# FATES case 2
# FATES d2blmax_pwr, biogeochem/FatesAllometryMod.F90
f_leaf_biomass_dbh_power <- function(.) {
  .super$pars$ldbh_p1 * .super$state$dbh^.super$pars$ldbh_p2 / .super$pars$carbon_to_biomass
}

# FATES dh2blmax_2pwr, biogeochem/FatesAllometryMod.F90
# FATES case 3
f_leaf_biomass_dbhmax_power <- function(.) {
  .super$pars$ldbh_p1 * 
    pmin(.super$state$dbh,.super$pars$dbh_max)^.super$pars$ldbh_p2 / .super$pars$carbon_to_biomass
}



# crown area ~ dbh
#########################

# - FATES subroutine: carea_2pwr, biogeochem/FatesAllometryMod.F90
# - for some other spread related process, current min and max are equal at 0.768654
# - spreadterm = spread * d2ca_max + (1._r8 - spread) * d2ca_min
# - 1.56 is the leaf biomass dbh exponent, maintains constant lai / canopy depth, but what about trimming / sla variation
# - dbh in cm, carea in m2

f_crown_area_fates <- function(.) {
  .super$state_pars$crwnarea_norm * (.super$state$dbh)^.super$pars$crwnarea_p1  
}

f_crown_area_purves2008 <- function(.) {
  pi * (.super$state$dbh*.super$pars$spread)^.super$pars$crwnarea_p1  
}

f_crown_radius <- function(.) {
  (.super$state$crown_area^0.5) / pi
}

f_crown_depth_fates <- function(.) {
  .super$state$height * .super$pars$crwndepth_frac  
}

# from FATES - calculates norm constant for crown area allometry considering a canopy packing term
# -- in FATES calculated parameter is called spreadterm and canopy_packing here is called spread 
f_canopy_packing <- function(.) {
  .super$env$canopy_packing * .super$pars$crwnarea_dbh_max + (1.0 - .super$env$canopy_packing) * .super$pars$crwnarea_dbh_max
}



# agb ~ dbh, and other biomass fractions 
#########################

# Saldarriaga et al 1998
# FATES case 1
f_abg_biomass_dbh_saldariaga1998 <- function(.) {
  .super$pars$frac_abg * .super$pars$abg_p1 * 
    .super$state$height^.super$pars$abg_p2 * .super$state$dbh^.super$pars$abg_p3 * 
    .super$pars$wood_density^.super$pars$abg_p4  
}

# power law
# FATES case 2
f_abg_biomass_dbh_power <- function(.) {
  .super$state$dbh^.super$pars$abg_p2 * .super$pars$abg_p1 / .super$pars$carbon_to_biomass   
}

# Chave et al 2014
# FATES case 3 
f_abg_biomass_dbh_chave2014 <- function(.) {
  (.super$pars$wood_density * .super$state$dbh^2 * .super$state$height )^.super$pars$abg_p2 * 
    .super$pars$abg_p1 / .super$pars$carbon_biomass   
}



# other biomass fractions with less functional variation in FATES 
#########################

f_fineroot_biomass_leaffrac <- function(.) {
  .super$state$leaf_biomass * .super$pars$fineroot_to_leaf_ratio
}

f_blg_biomass_dbh_fracabg <- function(.) {
  .super$state$abg_wood_biomass * (1-.super$pars$abg_frac) / .super$pars$abg_frac
}

f_wood_biomass <- function(.) {
  .super$state$abg_wood_biomass + .super$state$blg_wood_biomass 
}



# LAI
#########################

# - leaf biomass in kg C, area in m2, sla in g C m-2
f_crown_lai_linear <- function(.) {
   # this could be a function that calculates LAI including sla gradient
  .super$state_pars$leaf_area  <- .super$state$leaf_biomass * 1e3 * .super$pars$slatop   
  .super$state_pars$leaf_area / .super$state$crown_area
}

# biomass in g C, sla in m2 g-1 C 
f_lai_constantsla <- function(.) leaf_biomass * .super$pars$slatop  

f_lai_slagradient <- function(.) {
  leaf_biomass * .super$pars$slatop
} 





# print values
f_print_out_textonly <- function(.) {
  print(.super$state$text, quote=F )
  #print(.super$state$text)
}

f_print_out_textandval <- function(.) {
  print(paste(.super$state$text,.super$state$calcval), quote=F )
}



### END ###
