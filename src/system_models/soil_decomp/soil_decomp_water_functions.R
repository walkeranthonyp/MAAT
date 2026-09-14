################################
#
# MAAT soil_decomp process representation functions (PRFs)
# water and soil/litter environment scaling functions 
# 
# Matt Craig, AWalker October 2019 
#
################################

f_wcor_none <- function(...) 1
f_scor_none <- function(...) 1


################################
# soil/litter trait scalar functions
 
# power function of clay % (APW: I think) 
f_scor_sulman <- function(., i ) (.super$env$clay/.super$pars$clayref) ^ .super$pars$qslope_exp


# linear function of sand %
# Equation B3, Abramoff et al. 2022
f_scor_century_texture <- function(., i ) 
  .super$pars$scor_texture_int + .super$pars$scor_texture_slope*(1.0-.super$env$sand)


# exponential function of litter lignin proportion  
# Equation B4, Abramoff et al. 2022
f_scor_century_quality <- function(., i ) 
  exp(.super$pars$scor_quality_exp*.super$env$lignin)


# exponential function of litter quality 
f_scor_wieder_quality <- function(., i )
  .super$state_pars$anpp_mod * exp(.super$pars$k_exp[[i]]*.super$state_pars$quality_mod) 
 
 
# exponential function of soil clay 
f_scor_wieder_texture <- function(., i )
  .super$pars$k_desorp_tuning * exp(.super$pars$k_exp[[i]]*.super$env$clay) 



################################
# water scalar functions

# diffusion limitation power law 
# Ghezzehei et al. 2018, MillennialV2
f_wcor_ghezzehei_diffusion <- function(., i ) {
  (.super$env$vwc/.super$env$porosity) ^ .super$pars$wcor_diffusion_exp
}


# diffusion limitation power law 
# Ghezzehei et al. 2018, MillennialV2
# APW: assumes matpot is abs(matpot), inconsistent with ELM which is expressed <1
# APW: inconsistent in using VWC and matpot as independent env vars, one should be calculated from other 
f_wcor_ghezzehei_biological <- function(., i ) {
  exp(.super$pars$wcor_bio_lambda*-.super$env$matpot) * 
    (.super$pars$wcor_bio_kamin + (1-.super$pars$wcor_bio_kamin) * 
    ((.super$env$porosity-.super$env$vwc)/.super$env$porosity) ^ .super$pars$wcor_bio_exp) * 
    (.super$env$vwc/.super$env$porosity) ^ .super$pars$wcor_bio_exp
}


# used in ELM for both CTC & CENTURY 
f_wcor_andren1987 <- function(., ... ) {
  soil_matpot <- .super$env$matpot
  if(soil_matpot >= .super$state_pars$matpot_sat) { 
    1
  } else if(soil_matpot <= .super$env$matpot_min) {
    0 
  } else 
    log(.super$env$matpot_min/soil_matpot) / log(.super$env$matpot_min/.super$state_pars$matpot_sat)  
}


# ELM soil saturated matric/water potential
f_matpot_sat_elm_cosby1984_tab5 <- function(., ... ) {
  #matpot_sat <- 10 * 10^(1.88 - 0.0131*.super$env$sand)
  # APW: sonvert sand to proportion
  matpot_sat <- 10 * 10^(1.88 - 1.31*.super$env$sand)

  # convert from mm to MPa
  matpot_sat * .super$pars$conv_mm_to_MPa 
}


# Unimodal function
# MEC: this function is not normalized (e.g. 0-1)
# APW: function gives 0-0.02 for values of theta 0-1  
# APW: I assume porosity is saturated VWC and in the same units as VWC?
f_wcor_sulman <- function(., i ) {
  theta <- .super$env$vwc/.super$env$porosity
  #theta^3 * (1-theta)^2.5
  theta^.super$pars$wcor_unimodal_exp * (1-theta)^(.super$pars$wcor_unimodal_exp+.super$pars$wcor_unimodal_exp_diff)
}

# APW: this puts above on scale 0-1
f_wcor_sulman_normalized <- function(.,C,t,i){
  # volumetric water content
  vwc <- .super$env$vwc
  # porosity could be considered as analogous to water-holding capacity I believe
  porosity <- .super$env$porosity
  theta <- vwc/porosity
  
  # calculate max value for normalization
  # default Vmax from CORPSE will also need to be normalized by this value
  theta_max = 3/(2.5*(1+3/2.5))  
  
  #calculate wcor at theta_max
  wcor_max = theta_max^3 * (1-theta_max)^2.5
  
  (theta^3 * (1-theta)^2.5)/wcor_max
}


# inverse exponential
# APW: I assume porosity is saturated VWC
f_wcor_abramoff <- function(., i )
  1 / (1 + .super$pars$wcor_century1*exp(-.super$pars$wcor_century2*.super$env$vwc/.super$env$porosity))




###########################################
# APW: functions below this line not currently used

f_wcor_skopp <- function(.,C,t,i){
  #sourced from SoilR 
  #SoilR::fW.Skopp
  vwc = .super$env$vwc
  porosity = .super$env$porosity
  rwc = vwc/porosity #relative water content
  alpha = 2 #empirical parameter
  beta = 2  #empirical parameter
  f = 1.3  #empirical parameter
  g = 0.8 #empirical parameter
  
  pmin(alpha * rwc^f, beta * (1 - rwc)^g)
}


f_wcor_daycent2 <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fW.Daycent2
  # APW: this doesn't seem quite right unless vwc is in units of proportion of plant available water, i.e. betweend WP and FC
  W = .super$env$vwc * 100    #volumetric water content (percentage)
  WP = 0       #WP - A scalar representing the wilting point in percentage.
  FC = 100     # FC - A scalar representing the field capacity in percentage.
  RWC = (W - WP) * 100/(FC - WP)
  fRWC = 5 * (0.287 + (atan(pi * 0.009 * (RWC - 17.47)))/pi)
  fRWC_max = 5 * (0.287 + (atan(pi * 0.009 * (100 - 17.47)))/pi)
  fRWC/fRWC_max
}


f_wcor_daycent1 <- function(.,C,t,i) {
  #sourced from SoilR
  #SoilR::fW.Daycent1
  swc = .super$env$vwc         #A scalar or vector with soil water content of a soil layer (cm).
  a=0.6       #Empirical coefficient. For fine textured soils a = 0.6. For coarse textured soils a = 0.55.
  b=1.27      #Empirical coefficient. For fine textured soils b = 1.27. For coarse textured soils b = 1.70.
  c=0.0012    #Empirical coefficient. For fine textured soils c = 0.0012. For coarse textured soils c = -0.007.
  d=2.84       #Empirical coefficient. For fine textured soils d = 2.84. For coarse textured soils d = 3.22.
  partd=2.65  #Particle density of soil layer.
  bulkd=1     #Bulk density of soil layer (g/cm^3).
  width=1      #Thickness of a soil layer (cm).
  porespace=1-(bulkd/partd)
  wfps=(swc/width)*(1/porespace)
  (((wfps-b)/(a-b))^(d*((b-a)/(a-c))))*((wfps-c)/(a-c))^d
}


f_wcor_moyano <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fW.Moyano
  theta = .super$env$vwc #theta - volumetric water content
  a = 3.11 #empirical parameter
  b = 2.42 #empirical parameter
  a * theta - b * theta^2
}


f_wcor_gompertz <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fW.Gompertz
  theta = .super$env$vwc   #volumetric soil water content
  a = 0.824 #empirical parameter
  b = 0.309 #empirical parameter
  exp(-exp(a-b*theta*100))
}


f_wcor_candy <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fW.Candy
  theta = .super$env$vwc #volumetric soil water content
  PV = .super$env$porosity #pore volume
  Mi=theta/PV 
  ifelse(Mi<=0.5, 4*Mi*(1-Mi),1)
}


f_wcor_standcarb <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fW.Candy
  #uses limitation due to water potential and limitation due to oxygen diffusion 
  #to calculate overall limitation of moisture on decomposition rates
  Moist = .super$env$vwc *100 # moisture content of a litter or soil pool (%)
  MatricShape = 5 # scalar that determines when matric limit is reduced to the point that decay can begin to occur
  MatricLag = 0 #scalar used to offset the curve to the left or right
  MoistMin = 15 #scalar determining the minimum moisture content (MC: switched default to 15 for soil as done in soilR)
  MoistMax = 100 #scalar determining the maximum moisture content without diffusion limitations (MC: switched default to 100 for soil as done in soilR)
  DiffuseShape = 15 #scalar that determines the range of moisture contents where diffusion is not limiting
  DiffuseLag=4  #scalar used to shift the point when moisture begins to limit diffusion
  IncreaseRate=3/MoistMin  
  MatricLimit=(1-exp(-IncreaseRate*(Moist+MatricLag)))^MatricShape
  DiffuseLimit=exp(-1*(Moist/(MoistMax+DiffuseLag))^DiffuseShape)
  MatricLimit*DiffuseLimit
}


f_wcor_century <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fW.Candy
  #Calculates the effects of precipitation and potential evapotranspiration on decomposition rates.
  PPT = .super$env$precip_monthly       #monthly precipitation
  PET = .super$env$pet_monthly          #potential evapotranspiration
  1/(1+30*exp(-8.5*(PPT/PET)))
}




#f_wcor_rothc
#need to add this one which is also in soilR
#tricky think about this function is that monthly moisture limitation
#depends on value from previous month.
# APW: just add that previous month water index as an env variable -- a higher level model that runs a whole ecosystem could calulcte in the future


# Calculate matric potential from volumetric soil water content
# van Genuchten 1980?
f_conv_water_vwc_to_mp_vanGenuchten <- function(.,C,t,i) {
  SWC0 = .super$env$vwc
  SWCres = 0.108 # MEC: probably should store these parameters elsewhere, but they are constant in MEND
  SWCsat = 0.6 
  alpha = 0.035 
  n = 1.445
  SWPmin = .super$env$SWP_min   # MEC: not sure why this is an argument in MEND code; It's not used in function

    const_cm2MPa = 98e-6
    rlim = 1.01
    m = 1-1/n
    if(SWC0<=SWCres*rlim) {
      SWC = SWCres*rlim
    } else {
      SWC = SWC0
    }
    if(SWC < SWCsat){
      eff_sat = (SWC - SWCres)/(SWCsat - SWCres)
      fSWC2SWP = (1/(eff_sat^(1/m)) - 1)^(1/n)/alpha
      fSWC2SWP = -1*fSWC2SWP*const_cm2MPa
    } else {
      fSWC2SWP = 0
    }
    fSWC2SWP
}
# fSWC2SWP <- function(SWC0,SWCres,SWCsat,alpha,n,SWPmin){
#   const_cm2MPa = 98e-6
#   rlim = 1.01
#   m = 1-1/n
#   if(SWC0<=SWCres*rlim) {
#     SWC = SWCres*rlim
#   } else {
#     SWC = SWC0
#   }
#   if(SWC < SWCsat){
#     eff_sat = (SWC - SWCres)/(SWCsat - SWCres)
#     fSWC2SWP = (1/(eff_sat^(1/m)) - 1)^(1/n)/alpha
#     fSWC2SWP = -1*fSWC2SWP*const_cm2MPa
#   } else {
#     fSWC2SWP = 0
#   }
#   fSWC2SWP
# }
  
  
  
### END ### 
