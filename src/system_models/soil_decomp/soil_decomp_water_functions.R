f_wcor_sulman_normalized <- function(.,C,t,i){
  #volumetric water content
  vwc <- .super$env$vwc
  #porosity could be considered as analogous to water-holding capacity I believe
  porosity <- .super$env$porosity
  theta <- vwc/porosity
  
  # calculate max value for normalization
  # default Vmax from CORPSE will also need to be normalized by this value
  theta_max = 3/(2.5*(1+3/2.5))  
  
  #calculate wcor at theta_max
  wcor_max = theta_max^3 * (1-theta_max)^2.5
  
  (theta^3 * (1-theta)^2.5)/wcor_max
}

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

f_wcor_abramoff <- function(.,C,t,i){
  w1 = 30
  w2 = 9
  (1 / (1+w1*exp(-w2*.super$env$vwc/.super$env$porosity)))
  #.35 is whc I think? this should be specified in env. Perhaps porosity as in sulman.
}

#f_wcor_rothc
#need to add this one which is also in soilR
#tricky think about this function is that monthly moisture limitation
#depends on value from previous month.



