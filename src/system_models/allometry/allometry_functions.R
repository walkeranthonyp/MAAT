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
  10^(log10(.super$state_pars$dbh_effective) * (.super$pars$hdbh_p1+.super$pars$hdbh_p2)) 
}

# Poorter et al 2006
# FATES case 2 
f_height_hdbh_poorter2006 <- function(.) {
  .super$pars$hdbh_p1 * 
    (1 - exp(.super$pars$hdbh_p2 * .super$state_pars$dbh_effective^.super$pars$hdbh_p3) )  
}

# Power Law
# FATES case 3
f_height_hdbh_power <- function(.) {
  .super$pars$hdbh_p1 * .super$state_pars$dbh_effective^.super$pars$hdbh_p2 
}

# Chave et al 2014 
# FATES case 4 
f_height_hdbh_chave2014 <- function(.) {
  #p1e = p1   ! -eclim (assumed that p1 already has eclim removed)
  dbh_eff <- .super$state_pars$dbh_effective
  exp(.super$pars$hdbh_p1 + .super$pars$hdbh_p2*log(dbh_eff) + .super$pars$hdbh_p3*log(dbh_eff)^2 ) 
}

# Martinez Cano et al 2019
# FATES case 5 
f_height_hdbh_martinezcano2019 <- function(.) {
  (.super$pars$hdbh_p1 * .super$state_pars$dbh_effective^.super$pars$hdbh_p2) / 
    (.super$pars$hdbh_p3 + .super$state_pars$dbh_effective^.super$pars$hdbh_p2)  
}



# leaf biomass ~ dbh 
#########################

# - leaf biomass in kg C, dbh in cm

# Saldarriaga et al 1998 ???? -- doesn't seem to be
# FATES case 1
# FATES d2blmax_salda, biogeochem/FatesAllometryMod.F90
f_leaf_biomass_dbh_salda <- function(., dbh=.super$state$dbh ) {
  .super$pars$ldbh_p1 * dbh^.super$pars$ldbh_p2 * 
    .super$pars$wood_density^.super$pars$ldbh_p3
}

# power 
# FATES case 2
# FATES d2blmax_pwr, biogeochem/FatesAllometryMod.F90
f_leaf_biomass_dbh_power <- function(., dbh=.super$state$dbh ) {
  .super$pars$ldbh_p1 * dbh^.super$pars$ldbh_p2 / .super$pars$carbon_to_biomass
}

# FATES dh2blmax_2pwr, biogeochem/FatesAllometryMod.F90
# FATES case 3
f_leaf_biomass_dbhmax_power <- function(., dbh=.super$state$dbh ) {
  .super$pars$ldbh_p1 * 
    dbh^.super$pars$ldbh_p2 / .super$pars$carbon_to_biomass
}



# crown area ~ dbh
#########################

# - FATES subroutine: carea_2pwr, biogeochem/FatesAllometryMod.F90
# - for some other spread related process, current min and max are equal at 0.768654
# - spreadterm = spread * d2ca_max + (1._r8 - spread) * d2ca_min
# - 1.56 is the leaf biomass dbh exponent, maintains constant lai / canopy depth, but what about trimming / sla variation
# - dbh in cm, carea in m2

f_crown_area_fates <- function(., dbh=.super$state$dbh ) {
  #.super$state_pars$crwnarea_norm * dbh^.super$pars$crwnarea_p2  
  .super$state_pars$crwnarea_norm * dbh^.super$pars$ldbh_p2  
}

# check this APW
f_crown_area_purves2008 <- function(., dbh=.super$state$dbh ) {
  pi * (dbh*.super$state_pars$crownarea_norm)^.super$pars$crwnarea_p1  
}

f_crown_radius <- function(.) {
  (.super$state$crown_area / pi)^0.5
}

f_crown_depth_fates <- function(.) {
  .super$state$height * .super$pars$crwndepth_frac  
}

# from FATES - calculates norm constant for crown area allometry considering a canopy openness term
# -- in FATES calculated parameter is called spreadterm and canopy_openness here is called spread 
f_crwnallom_canopy_mod <- function(.) {
  .super$env$canopy_openness * .super$pars$crwnarea_dbh_max + (1.0 - .super$env$canopy_openness) * .super$pars$crwnarea_dbh_min
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
    .super$pars$abg_p1 / .super$pars$carbon_to_biomass   
}



# other biomass fractions with less functional variation in FATES 
#########################

f_fineroot_biomass_leaffrac <- function(.) {
  .super$state$leaf_biomass * .super$pars$fineroot_to_leaf_ratio
}

f_blg_biomass_dbh_fracabg <- function(.) {
  .super$state$abg_wood_biomass * (1-.super$pars$frac_abg) / .super$pars$frac_abg
}

f_wood_biomass <- function(.) {
  .super$state$abg_wood_biomass + .super$state$blg_wood_biomass 
}



# LAI
#########################

# leaf_biomass in g C, sla in m2 g-1 C 
f_lai_linear_constant_sla <- function(., leaf_biomass=.super$state_pars$leaf_biomass_perarea, sla=.super$pars$slatop ) {
  leaf_biomass * sla 
}

# leaf_biomass in g C, sla in m2 g-1 C 
f_lai_exponential_sla <- function(., leaf_biomass=.super$state_pars$leaf_biomass_perarea ) {

  # if leaf_biomass becomes too large, lai becomes an imaginary number
  # -- because lai equation would require log of < 0
  clim <- (exp(-1.0 * .super$pars$k * .super$env$canopy_lai_above)) / 
    (.super$pars$k * .super$pars$slatop)
  if(leaf_biomass >= clim) 
    stop('FATAL ERROR: allometry, f_lai_kovenock_sla_exponential: leaf_biomass too high ')

  (log(exp(-1.0 * .super$pars$k * .super$env$canopy_lai_above) - 
    .super$pars$k * .super$pars$slatop * leaf_biomass) + 
    (.super$pars$k * .super$env$canopy_lai_above)) / (-1.0 * .super$pars$k)
} 

# leaf_biomass in g C, sla in m2 g-1 C 
f_lai_kovenock_exponential_sla <- function(.) {

  # Leaf biomass per unit area at which slamax is reached 
  .super$state_pars$leaf_biomass_perarea_slamax <- 
    (.super$pars$slatop - .super$pars$slamax * exp(-1.0 * .super$pars$k * .super$env$canopy_lai_above)) / 
    (-1.0 * .super$pars$k * .super$pars$slatop * .super$pars$slamax)
  if(.super$state_pars$leaf_biomass_perarea_slamax < 0.0) .super$state_pars$leaf_biomass_perarea_slamax <- 0.0

  # calculate exponential and linear components separately and sum
  if(.super$state_pars$leaf_biomass_perarea > .super$state_pars$leaf_biomass_perarea_slamax) {
    lai_exp <- .$lai_exponential(leaf_biomass=.super$state_pars$leaf_biomass_perarea_slamax)
    lai_lin <- .$lai_linear(leaf_biomass=(.super$state_pars$leaf_biomass_perarea - .super$state_pars$leaf_biomass_perarea_slamax), 
                            sla=.super$pars$slamax )
  } else {
    lai_exp <- .$lai_exponential(leaf_biomass=.super$state_pars$leaf_biomass_perarea)
    lai_lin <- 0
  }

  lai_exp + lai_lin
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
