###########################
#
# allometry for MAAT object unit testing
# 
# AWalker March 2017
#
###########################



### Load model scripts 
###############################

# allometry
source('allometry_object.R')
allometry_object$.test()
allometry_object$.test(dbh=1)
allometry_object$.test(allometry.lai='f_lai_linear_constant_sla')

allometry_object$.test(verbose=T)

allometry_object$fnames
allometry_object$fns
allometry_object$fns$lai
allometry_object$state
allometry_object$state_pars
allometry_object$pars
allometry_object$env
allometry_object$dataf
allometry_object$run()


source('allometry_object.R')
allometry_object$.test_pars()
allometry_object$.test_pars(dbh=1)
allometry_object$.test_pars(wood_density=1)

library(lattice)
source('allometry_object.R')
allometry_object$.test_dbh()
allometry_object$.test_dbh(wood_density=1)
allometry_object$dataf
allometry_object$output()
allometry_object$state
allometry_object$fns$crwnallom()


# https://fates-docs.readthedocs.io/en/latest/fates_tech_note.html#cohort-state-variables
# 
#  FATES defaults:
#  all default to 1, can be varied by PFT
#  amode -- ABG biomass --
#  cmode -- coarse root biomass --
#  fmode -- fine root biomass -- 
#  lmode -- leaf biomass -- 
#  hmode -- height -- 
#  smode -- sapwood biomass -- 
#  stmode -- storage allometry -- 
#  
#  PFT
#  fates_pftname =
#    "broadleaf_evergreen_tropical_tree          ",
#    "needleleaf_evergreen_extratrop_tree        ",
#    "needleleaf_colddecid_extratrop_tree       ",
#    "broadleaf_evergreen_extratrop_tree         ",
#    "broadleaf_hydrodecid_tropical_tree         ",
#    "broadleaf_colddecid_extratrop_tree        ",
#    "broadleaf_evergreen_extratrop_shrub        ",
#    "broadleaf_hydrodecid_extratrop_shrub       ",
#    "broadleaf_colddecid_extratrop_shrub       ",
#    "arctic_c3_grass                            ",
#    "cool_c3_grass                              ",
#    "c4_grass                                   " ;
# fates_alloc_storage_cushion = 1.2, 1.2, 1.2, 1.2, 2.4, 1.2, 1.2, 2.4, 1.2, 1.2, 1.2, 1.2 
# fates_alloc_store_priority_frac = 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 ;
# fates_allom_agb1          = 0.06896, 0.06896, 0.06896, 0.06896, 0.06896, 0.06896, 0.06896, 0.06896, 0.06896, 0.01, 0.01, 0.01 ;
# fates_allom_agb2          = 0.572, 0.572, 0.572, 0.572, 0.572, 0.572, 0.572, 0.572, 0.572, 0.572, 0.572, 0.572 ;
# fates_allom_agb3          = 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94 ;
# fates_allom_agb4          = 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931 ;
# fates_allom_agb_frac      = 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 ;
# fates_allom_blca_expnt_diff = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
# fates_allom_crown_depth_frac = 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.95, 0.95, 0.95, 1, 1, 1 ;
# fates_allom_d2bl1         = 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07 ;
# fates_allom_d2bl2         = 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3 ;
# fates_allom_d2bl3         = 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55 ;
# fates_allom_d2ca_coefficient_max = 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464, 0.6568464 ;
# fates_allom_d2ca_coefficient_min = 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119, 0.3381119 ;
# fates_allom_d2h1          = 0.64, 0.64, 0.64, 0.64, 0.64, 0.64, 0.64, 0.64, 0.64, 0.64, 0.64, 0.64 ;
# fates_allom_d2h2          = 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37 ; 
# fates_allom_d2h3          = -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9 ;
# fates_allom_dbh_maxheight = 90, 80, 80, 80, 90, 80, 3, 3, 2, 0.35, 0.35, 0.35 ;
# fates_wood_density        = 0.7, 0.4, 0.7, 0.53, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 ;
# fates_leaf_slamax         = 0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.012, 0.03, 0.03, 0.03, 0.03, 0.03 ;
# fates_leaf_slatop         = 0.012, 0.005, 0.024, 0.009, 0.03, 0.03, 0.012, 0.03, 0.03, 0.03, 0.03, 0.03 ;
# 

source('allometry_object.R')
# default FATES broadleaf evergreen tropical tree
allometry_object$.test_dbh()
# FATES needleleaf evergreen extratropical tree
allometry_object$.test_dbh(dbh_max=80, slatop=0.005, wood_density=0.4 )
# FATES needleleaf cold-decid extratropical tree
allometry_object$.test_dbh(dbh_max=80, slatop=0.024 )
# FATES broadleaf evergreen extratropical tree
allometry_object$.test_dbh(dbh_max=80, slatop=0.009, wood_density=0.54 )
# FATES broadleaf cold-decid extratropical tree
allometry_object$.test_dbh(dbh_max=80, slatop=0.03 )



# Rob Tunison -- Puerto Rico
# 	fates_pftname 
# 	"earlysucc_evergreen_tropical_tree", 
# 	"midsucc_evergreen_tropical_tree", 
# 	"latesucc_evergreen_tropical_tree", 
# 	"c4_grass" 

# 	fates_allom_amode   = 3,3,3,1 
#   fates_allom_agb1   = 0.09136, 0.1192, 0.07597, 0.01 
#   fates_allom_agb2   = 0.985, 0.897, 0.95, 0.572 
#   fates_allom_agb3   = 1.94, 1.94, 1.94, 1.94
#   fates_allom_agb4   = 0.931, 0.931, 0.931, 0.931 
#   fates_allom_lmode   = 3,3,3,1 
#   fates_allom_d2bl1   = 0.1266844, 0.1266844, 0.1266844, 0.07
#   fates_allom_d2bl2   = 1.281329, 1.281329, 1.281329, 1.3 
#   fates_allom_d2bl3   = 0.55, 0.55, 0.55, 0.55 
#   fates_allom_agb_frac   = 0.8, 0.7, 0.6, 0.6 
#   fates_allom_crown_depth_frac   = 0.3, 0.5, 0.5, 1 
#   fates_allom_hmode   = 5,5,5,1 
#   fates_allom_d2h1   = 22.38, 57.6, 30.6, 0.64
#   fates_allom_d2h2   = 2.275, 0.74, 0.661, 0.37
#   fates_allom_d2h3   = 4.804, 21.6, 46.795, -999.9 
#   fates_allom_dbh_maxheight   = 70, 80, 90, 0.35 
#   fates_allom_la_per_sa_int   = 1.0, 0.8, 0.8, 0.8 
#   fates_allom_l2fr   = 0.8, 0.8, 0.7, 0.8 
#   fates_wood_density   = 0.3, 0.5, 0.6, 0.7
#   fates_leaf_slatop   = 0.012, 0.02, 0.024, 0.03


source('allometry_object.R')
# default FATES broadleaf evergreen tropical tree
allometry_object$.test_dbh()
# FATES-RT earlysucc_evergreen_tropical_tree -- Cecropia
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  fnames.height='f_height_hdbh_martinezcano2019',
  abg_p1=0.09136, abg_p2=0.985, abg_p3=1.94, abg_p4=0.931, 
  ldbh_p1=0.1266844, ldbh_p2=1.281329, ldbh_p3=0.55,
  hdbh_p1=22.38, hdbh_p2=2.275, hdbh_p3=4.804, 
  dbh_max=70, fineroot_to_leaf_ratio=0.8, wood_density=0.3, slatop=0.012,
  crwndepth_frac=0.3, frac_abg=0.8
  )

# FATES-RT midsucc_evergreen_tropical_tree 
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  fnames.height='f_height_hdbh_martinezcano2019',
  abg_p1=0.01192, abg_p2=0.897, abg_p3=1.94, abg_p4=0.931, 
  ldbh_p1=0.1266844, ldbh_p2=1.281329, ldbh_p3=0.55,
  hdbh_p1=57.6, hdbh_p2=0.74, hdbh_p3=21.6, 
  dbh_max=80, fineroot_to_leaf_ratio=0.8, wood_density=0.5, slatop=0.02,
  crwndepth_frac=0.5, frac_abg=0.7
)

# FATES-RT latesucc_evergreen_tropical_tree 
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  fnames.height='f_height_hdbh_martinezcano2019',
  abg_p1=0.07597, abg_p2=0.95, abg_p3=1.94, abg_p4=0.931, 
  ldbh_p1=0.1266844, ldbh_p2=1.281329, ldbh_p3=0.55,
  hdbh_p1=30.6, hdbh_p2=0.661, hdbh_p3=46.795, 
  dbh_max=90, fineroot_to_leaf_ratio=0.7, wood_density=0.6, slatop=0.024,
  crwndepth_frac=0.5, frac_abg=0.6
)



### END ###