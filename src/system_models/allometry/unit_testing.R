###########################
#
# allometry for MAAT object unit testing
# 
# AWalker March 2017
#
###########################


library(lattice)

### Load model scripts 
###############################

# allometry
source('allometry_object.R')
allometry_object$.test()
allometry_object$.test(dbh=1)
allometry_object$.test(dbh=20)
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

source('allometry_object.R')
allometry_object$.test_dbh()
allometry_object$.test_dbh(wood_density=1)
allometry_object$dataf
allometry_object$output()
allometry_object$state
allometry_object$fns$crwnallom
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

# FATES 11 March 2026
#  fates_pftname =
# "broadleaf_evergreen_tropical_tree", "needleleaf_evergreen_extratrop_tree", "needleleaf_colddecid_extratrop_tree", 
# "broadleaf_evergreen_extratrop_tree", "broadleaf_hydrodecid_tropical_tree", "broadleaf_colddecid_extratrop_tree", 
# "broadleaf_evergreen_extratrop_shrub", "broadleaf_hydrodecid_extratrop_shrub", "broadleaf_colddecid_extratrop_shrub", 
# "broadleaf_evergreen_arctic_shrub", "broadleaf_colddecid_arctic_shrub", "arctic_c3_grass", "cool_c3_grass", "c4_grass"
# fates_alloc_storage_cushion      = [1.2, 1.2, 1.2, 1.2, 2.4, 1.2, 1.2, 2.4, 1.2, 1.5, 1.4, 1.2, 1.2, 1.2]
# fates_alloc_store_priority_frac  = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.6, 0.6, 0.8, 0.8]
# fates_allom_agb1                 = [0.0673, 0.1364012, 0.0393057, 0.2653695, 0.0673, 0.0728698, 0.06896, 0.06896, 0.06896, 0.06896, 0.06896, 0.001, 0.001, 0.003]
# fates_allom_agb2                 = [0.976, 0.9449041, 1.087335, 0.8321321, 0.976, 1.0373211, 0.572, 0.572, 0.572, 0.5289883, 0.6853945, 1.6592, 1.6592, 1.3456]  
# fates_allom_agb3                 = [1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 2.1010352, 1.7628613, 1.248, 1.248, 1.869]
# fates_allom_agb4                 = [0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, -999.9, -999.9, -999.9]
# fates_allom_agb_frac             = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 1.0, 1.0, 1.0]
# fates_allom_blca_expnt_diff      = [-0.12, -0.34, -0.32, -0.22, -0.12, -0.35, 0.0, 0.0, 0.0, 0.0, 0.0, -0.487, -0.487, -0.259]
# fates_allom_crown_depth_frac     = 
# fates_allom_d2bl1                = [0.04, 0.07, 0.07, 0.01, 0.04, 0.07, 0.07, 0.07, 0.07, 0.0481934, 0.0481934, 0.0004, 0.0004, 0.0012]
# fates_allom_d2bl2                = [1.6019679, 1.5234373, 1.3051237, 1.9621397, 1.6019679, 1.3998939, 1.3, 1.3, 1.3, 1.0600586, 1.7176758, 1.7092, 1.7092, 1.5879]
# fates_allom_d2bl3                = [0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.3417, 0.3417, 0.9948]
# fates_allom_d2ca_coefficient_max = [0.2715891, 0.3693718, 1.0787259, 0.0579297, 0.2715891, 1.1553612, 0.6568464, 0.6568464, 0.6568464, 0.4363427, 0.3166497, 0.0408, 0.0408, 0.0862]
# fates_allom_d2ca_coefficient_min = [0.2715891, 0.3693718, 1.0787259, 0.0579297, 0.2715891, 1.1553612, 0.6568464, 0.6568464, 0.6568464, 0.4363427, 0.3166497, 0.0408, 0.0408, 0.0862]
# fates_allom_d2h1                 = [78.4087704, 306.842667, 106.8745821, 104.3586841, 78.4087704, 31.4557047, 0.64, 0.64, 0.64, 0.8165625, 0.778125, 0.1812, 0.1812, 0.3353]
# fates_allom_d2h2                 = [0.8124383, 0.752377, 0.9471302, 1.1146973, 0.8124383, 0.9734088, 0.37, 0.37, 0.37, 0.2316113, 0.4027002, 0.6384, 0.6384, 0.4235]
# fates_allom_d2h3                 = [47.6666164, 196.6865691, 93.9790461, 160.6835089, 47.6666164, 16.5928174, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9]
# fates_allom_dbh_maxheight        = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 3.0, 3.0, 2.0, 2.4, 1.9, 20.0, 20.0, 30.0]
# fates_wood_density               = [0.548327, 0.44235, 0.454845, 0.754336, 0.548327, 0.566452, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7]
# fates_recruit_height_min         = [1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 0.2, 0.2, 0.2, 0.8, 0.8, 0.11, 0.2, 0.2]
# fates_allom_amode                = [3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 5, 5, 5]
# fates_allom_lmode                = [2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 5, 5, 5]
# fates_allom_hmode                = [5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 3, 3, 3]
# fates_allom_cmode                = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# fates_allom_dmode                = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# fates_allom_fmode                = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# fates_allom_la_per_sa_int        = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
# fates_allom_la_per_sa_slp        = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# fates_leaf_slamax                = [0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.012, 0.03, 0.03, 0.012, 0.032, 0.05, 0.05, 0.05]
# fates_leaf_slatop                = [0.012, 0.005, 0.024, 0.009, 0.03, 0.03, 0.012, 0.03, 0.03, 0.01, 0.032, 0.027, 0.05, 0.05]

# amode: 1, Saldariaga 1998; 2, 2-param power law; 3, Chave 2014; 4, 3-param power law; 5, 3-param power law grass
# lmode: 1, Saldariaga 1998; 2, 2-param power law; 3, 2-param power law max; 4, 3-param power law max; 5, 3-param power law max grass
# -- ca: 1, 1,3,5, 2-param power law max dbh; 2, 2-param power law; 4, 3-param power law max DBH & height
# hmode: 1, O'Brien 1995; 2, Poorter 2006; 3, 2-param power law; 4, Chave 2014; 5, Martinez-Cano 2019
# cmode: 1, constant ratio 
# dmode: 1, 
# fmode: 1,  


source('allometry_object.R')
# default FATES broadleaf evergreen tropical tree
allometry_object$.test_dbh(write_data='FATES_260311_BlEv_trop_tree.csv')
# Lambir modified FATES broadleaf evergreen tropical tree
# - lmode = 3 
# - fmode = 2 (what is that?)
allometry_object$.test_dbh(
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  # froot mode?
  ldbh_p1=0.1267, ldbh_p2=1.2813, crwnarea_dbh_min=0.7687,
  hdbh_p1=90, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
  wood_density=0.84, slatop=0.02, slamax=0.04, frac_abg=0.8,  
  write_data='FATES_260311_BlEv_trop_tree_Lambir.csv')
# FATES needleleaf evergreen extra-tropical tree
allometry_object$.test_dbh(abg_p1=0.1364012, abg_p2=0.9449041, ldbh_edif=-0.34, 
                           crwnarea_dbh_min=0.3693718, ldbh_p1=0.07, ldbh_p2=1.5234373, 
                           hdbh_p1=306.84, hdbh_p2=0.752377, hdbh_p3=196.6866, 
                           wood_density=0.44235, slatop=0.005, write_data='FATES_260311_NlEv_extrop_tree.csv' )
# FATES needleleaf cold-deciduous extra-tropical tree
# @200 cm DBH: LAI ~0.14, crown radius ~43 m, height ~66 m, AGBw ~80,000; start LAI ~1.1
allometry_object$.test_dbh(abg_p1=0.0393057, abg_p2=1.087335, ldbh_edif=-0.32, 
                           crwnarea_dbh_min=1.0787259, ldbh_p1=0.07, ldbh_p2=1.3051237, 
                           hdbh_p1=106.8745821, hdbh_p2=0.9471302, hdbh_p3=93.9790461, 
                           wood_density=0.454845, slatop=0.024  )
# FATES broadleaf evergreen extratropical tree
# @200 cm DBH: LAI ~0.25, crown radius ~28 m, height ~73 m, AGBw ~25,000; start LAI ~1
allometry_object$.test_dbh(abg_p1=0.2653695, abg_p2=0.8321321, ldbh_edif=-0.22, 
                           crwnarea_dbh_min=0.0579297, ldbh_p1=0.01, ldbh_p2=1.9621397, 
                           hdbh_p1=104.3586841, hdbh_p2=1.1146973, hdbh_p3=160.6835089, 
                           wood_density=0.754336, slatop=0.009, write_data='FATES_260311_BlEv_extrop_tree.csv'  )
# default FATES broadleaf hydro-deciduous tropical tree
# @200 cm DBH: LAI ~1.3, crown radius ~28 m, height ~48 m, AGBw ~25,000; start LAI ~3.2
allometry_object$.test_dbh(slatop=0.03)
# FATES broad-leaf cold-deciduous extra-tropical tree
# @200 cm DBH: LAI ~0.15, crown radius ~62 m, height ~29 m, AGBw ~39,000; start LAI ~1.3
allometry_object$.test_dbh(abg_p1=0.0728698, abg_p2=1.0373211, ldbh_edif=-0.35, 
                           crwnarea_dbh_min=1.1553612, ldbh_p1=0.07, ldbh_p2=1.3998939, 
                           hdbh_p1=31.4557047, hdbh_p2=0.9734088, hdbh_p3=16.5928174, 
                           wood_density=0.566452, slatop=0.03, write_data='FATES_260311_BlDc_cold_extrop_tree.csv' )

# FATES broadleaf evergreen tropical tree -- RK BCI e2
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  abg_p1=0.0673,abg_p2=0.967,abg_p3=-9,abg_p4=-9,
  frac_abg=0.8, wood_density=0.4,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.1267, ldbh_p2=1.2813, ldbh_p3=-9,
  slatop=0.02, slamax=0.04,
  crwnarea_dbh_min=0.7687, crwnarea_dbh_max=0.7687,
  crwnarea_p2=1.2813,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=57.6, hdbh_p2=0.74, hdbh_p3=21.6, dbh_max=200   
)

# FATES broadleaf evergreen tropical tree -- RK BCI e3
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  abg_p1=0.0673,abg_p2=0.967,abg_p3=-9,abg_p4=-9,
  frac_abg=0.8, wood_density=0.4,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.1267, ldbh_p2=1.2813, ldbh_p3=-9,
  slatop=0.02, slamax=0.04,
  crwnarea_dbh_min=0.7687, crwnarea_dbh_max=0.7687,
  crwnarea_p2=1.2813,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=85, hdbh_p2=0.99, hdbh_p3=70, dbh_max=300   
)

# FATES broadleaf evergreen tropical tree -- RK BCI e4
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  abg_p1=0.0673,abg_p2=0.967,abg_p3=-9,abg_p4=-9,
  frac_abg=0.8, wood_density=0.4,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.1267, ldbh_p2=1.2813, ldbh_p3=-9,
  slatop=0.02, slamax=0.04,
  crwnarea_dbh_min=0.7687, crwnarea_dbh_max=0.7687,
  crwnarea_p2=1.2813,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=85, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300   
)

# FATES broadleaf evergreen tropical tree -- RK BCI e4 tester
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  abg_p1=0.0673,abg_p2=0.967,abg_p3=-9,abg_p4=-9,
  frac_abg=0.8, wood_density=0.4,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.25, ldbh_p2=1.2813, ldbh_p3=-9,
  slatop=0.02, slamax=0.04,
  crwnarea_dbh_min=0.7687, crwnarea_dbh_max=0.7687,
  crwnarea_p2=1.2813,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=85, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300   
)


# for ensemble
source('allometry_object.R')
# default FATES broadleaf evergreen tropical tree
# @200 cm DBH: LAI ~0.5, crown radius ~28 m, height ~48 m, AGBw ~25,000
allometry_object$.test_dbh()  
# @200 cm DBH: LAI ~0.15, crown radius ~ 46 m, height ~48 m, AGBw ~25,000
allometry_object$.test_dbh(crwnarea_dbh_min=0.75)  
# @200 cm DBH: LAI ~1, crown radius ~21 m, height ~48 m, AGBw ~25,000; start LAI ~4
allometry_object$.test_dbh(crwnarea_dbh_min=0.1, ldbh_edif=-0.2 )  
# @200 cm DBH: LAI ~30, crown radius ~11 m, height ~48 m, AGBw ~25,000; start LAI ~2.5
allometry_object$.test_dbh(crwnarea_dbh_min=0.1, ldbh_edif=0.2 )  
# @200 cm DBH: LAI ~5, crown radius ~11 m, height ~48 m, AGBw ~25,000; start LAI ~3
allometry_object$.test_dbh(crwnarea_dbh_min=0.1, ldbh_edif=0.05 )  
# @200 cm DBH: LAI ~1.2, crown radius ~18 m, height ~48 m, AGBw ~25,000; start LAI ~0.9
allometry_object$.test_dbh(ldbh_edif=0.05 )  
# @200 cm DBH: LAI ~1.2, crown radius ~18 m, height ~52 m, AGBw ~28,000; start LAI ~0.9
allometry_object$.test_dbh(ldbh_edif=0.05, hdbh_p3=35 )  
# @200 cm DBH: LAI ~1.2, crown radius ~18 m, height ~35 m, AGBw ~19,000; start LAI ~0.9
allometry_object$.test_dbh(ldbh_edif=0.05, hdbh_p1=57 )  
# @200 cm DBH: LAI ~1, crown radius ~21 m, height ~35 m, AGBw ~19,000; start LAI ~1
allometry_object$.test_dbh(ldbh_edif=0.0, hdbh_p1=57 )  
# @200 cm DBH: LAI ~1.2, crown radius ~18 m, height ~35 m, AGBw ~25,000; start LAI ~0.9
allometry_object$.test_dbh(ldbh_edif=0.05, hdbh_p1=57, abg_p1=0.09 )  
# @200 cm DBH: LAI ~3.3, crown radius ~12.4 m, height ~48 m, AGBw ~25,000; start LAI ~3.3
allometry_object$.test_dbh(crwnarea_dbh_min=0.1, ldbh_edif=0 )  
# @200 cm DBH: LAI ~0.85, crown radius ~21.7 m, height ~35 m, AGBw ~25,000; start LAI ~0.70
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=0.05, hdbh_p1=57, abg_p1=0.09 )  
# @200 cm DBH: LAI ~1.86, crown radius ~18 m, height ~35 m, AGBw ~25,000; start LAI ~1.3
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=0.05, hdbh_p1=57, abg_p1=0.09 )  
# @200 cm DBH: LAI ~0.68 crown radius ~24 m, height 34.6 m, AGBw ~25,000; start LAI ~1.5
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=-0.12, hdbh_p1=57, abg_p1=0.09 )  
# @200 cm DBH: LAI 2.17, crown radius ~21.7 m, height ~35 m, AGBw ~25,000; start LAI ~1.51
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=0.05, hdbh_p1=57, abg_p1=0.09, ldbh_p1=0.09 )  
# @200 cm DBH: LAI 1.72, crown radius ~24.0 m, height ~54.7 m, AGBw ~25,000; start LAI ~4.42
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=-0.12, hdbh_p1=90, abg_p1=0.059, ldbh_p1=0.09 )  
# @200 cm DBH: LAI 16, crown radius ~12.7 m, height ~54.7 m, AGBw ~25,000; start LAI ~3.4
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=0.12, hdbh_p1=90, abg_p1=0.059, ldbh_p1=0.09 )  
# @200 cm DBH: LAI 3, crown radius ~12.7 m, height ~54.7 m, AGBw ~25,000; start LAI ~1.2
allometry_object$.test_dbh(crwnarea_dbh_min=0.02, ldbh_edif=0.12, hdbh_p1=90, abg_p1=0.059, ldbh_p1=0.04 )  
# @200 cm DBH: LAI 0.77, crown radius ~34.0 m, height ~34.6 m, AGBw ~24,700; start LAI ~1.7
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=-0.12, hdbh_p1=57, abg_p1=0.09, ldbh_p1=0.09 )  
# @200 cm DBH: LAI 1.06, crown radius ~34.0 m, height ~34.6 m, AGBw ~24,700; start LAI ~2.5
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=-0.12, hdbh_p1=57, abg_p1=0.09, ldbh_p1=0.12 )  


source('allometry_object.R')
# test 1
# min
# @200 cm DBH: LAI 1.06, crown radius ~34.0 m, height ~34.6 m, AGBw ~25,000; @1.3m height LAI ~2.16
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=-0.12, hdbh_p1=57, 
                           abg_p1=0.091149, ldbh_p1=0.12, frac_abg=0.8, write_data='test1_min.csv' )  
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=-0.12, hdbh_p1=57, 
                           abg_p1=0.091149, ldbh_p1=0.12, frac_abg=0.8, 
                           fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
                           ldbh_p2=1.2813, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
                           wood_density=0.84, slatop=0.02, slamax=0.04,
                           write_data='test1_min_lambir.csv' )  
# max
# @200 cm DBH: LAI 3, crown radius ~12.7 m, height ~54.7 m, AGBw ~25,000; @1.3m height LAI ~1.28
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=0.12, hdbh_p1=90, 
                           abg_p1=0.058364, ldbh_p1=0.04, write_data='test1_max.csv'  )  
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=0.12, hdbh_p1=90, 
                           abg_p1=0.058364, ldbh_p1=0.04, 
                           fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
                           ldbh_p2=1.2813, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
                           wood_density=0.84, slatop=0.02, slamax=0.04,
                           write_data='test1_max_lambir.csv'  )  
# mean
# @200 cm DBH: LAI 1.93, crown radius ~12.7 m, height ~44.7 m, AGBw ~26,300; @1.3m height LAI ~1.93
allometry_object$.test_dbh(crwnarea_dbh_min=0.3, ldbh_edif=0, hdbh_p1=73.5, 
                           abg_p1=0.07476, ldbh_p1=0.08, write_data='test1_mean.csv'  )  
allometry_object$.test_dbh(crwnarea_dbh_min=0.3, ldbh_edif=0, hdbh_p1=73.5, 
                           abg_p1=0.07476, ldbh_p1=0.08, 
                           fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
                           ldbh_p2=1.2813, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
                           wood_density=0.84, slatop=0.02, slamax=0.04,
                           write_data='test1_mean_lambir.csv'  )  

# test 2
# min
# @200 cm DBH: LAI 1.06, crown radius ~34.0 m, height ~34.6 m, AGBw ~25,000; @1.3m height LAI ~2.16
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=-0.12, hdbh_p1=57, 
                           abg_p1=0.091149, ldbh_p1=0.08, frac_abg=0.8, write_data='test2_min.csv'  )  
allometry_object$.test_dbh(crwnarea_dbh_min=0.4, ldbh_edif=-0.12, hdbh_p1=57, 
                           abg_p1=0.091149, ldbh_p1=0.08, frac_abg=0.8, 
                           fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
                           ldbh_p2=1.2813, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
                           wood_density=0.84, slatop=0.02, slamax=0.04,
                           write_data='test2_min_lambir.csv'  )  
# max
# @200 cm DBH: LAI 3, crown radius ~12.7 m, height ~54.7 m, AGBw ~25,000; @1.3m height LAI ~1.28
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=0.12, hdbh_p1=90, 
                           abg_p1=0.058364, ldbh_p1=0.02, write_data='test2_max.csv'  )  
allometry_object$.test_dbh(crwnarea_dbh_min=0.2, ldbh_edif=0.12, hdbh_p1=90, 
                           abg_p1=0.058364, ldbh_p1=0.02, 
                           fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
                           ldbh_p2=1.2813, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
                           wood_density=0.84, slatop=0.02, slamax=0.04,
                           write_data='test2_max_lambir.csv'  )  
# mean
# @200 cm DBH: LAI 1.93, crown radius ~12.7 m, height ~44.7 m, AGBw ~26,300; @1.3m height LAI ~1.93
allometry_object$.test_dbh(crwnarea_dbh_min=0.3, ldbh_edif=0, hdbh_p1=73.5, 
                           abg_p1=0.07476, ldbh_p1=0.05, write_data='test2_mean.csv'  )  
allometry_object$.test_dbh(crwnarea_dbh_min=0.3, ldbh_edif=0, hdbh_p1=73.5, 
                           abg_p1=0.07476, ldbh_p1=0.05, 
                           fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
                           ldbh_p2=1.2813, hdbh_p2=0.99, hdbh_p3=20, dbh_max=300, 
                           wood_density=0.84, slatop=0.02, slamax=0.04,
                           write_data='test2_mean_lambir.csv'  )  


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
  abg_p1=0.09136, abg_p2=0.985, abg_p3=1.94, abg_p4=0.931, 
  wood_density=0.3, frac_abg=0.8,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.1266844, ldbh_p2=1.281329, ldbh_p3=0.55,
  slatop=0.012,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=22.38, hdbh_p2=2.275, hdbh_p3=4.804, dbh_max=70, 
  fineroot_to_leaf_ratio=0.8, 
  crwndepth_frac=0.3
  )

# FATES-RT midsucc_evergreen_tropical_tree 
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  abg_p1=0.01192, abg_p2=0.897, abg_p3=1.94, abg_p4=0.931, 
  wood_density=0.5, frac_abg=0.7,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.1266844, ldbh_p2=1.281329, ldbh_p3=0.55,
  slatop=0.02,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=57.6, hdbh_p2=0.74, hdbh_p3=21.6, dbh_max=80, 
  fineroot_to_leaf_ratio=0.8, 
  crwndepth_frac=0.5
)

# FATES-RT latesucc_evergreen_tropical_tree 
allometry_object$.test_dbh(
  fnames.abg_biomass='f_abg_biomass_dbh_chave2014',
  abg_p1=0.07597, abg_p2=0.95, abg_p3=1.94, abg_p4=0.931, 
  wood_density=0.6, frac_abg=0.6,
  fnames.leaf_biomass='f_leaf_biomass_dbhmax_power',
  ldbh_p1=0.1266844, ldbh_p2=1.281329, ldbh_p3=0.55,
  slatop=0.024,
  fnames.height='f_height_hdbh_martinezcano2019',
  hdbh_p1=30.6, hdbh_p2=0.661, hdbh_p3=46.795, dbh_max=90, 
  fineroot_to_leaf_ratio=0.7, 
  crwndepth_frac=0.5 
)



### END ###