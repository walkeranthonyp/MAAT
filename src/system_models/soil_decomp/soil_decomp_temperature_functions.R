f_tcor_century1 <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.Century1
  Tmax <- 45
  Topt <- 35
  ((((Tmax - .super$env$temp)/(Tmax - Topt))^0.2) * exp((0.2/2.63) * (1 - ((Tmax - .super$env$temp)/(Tmax - Topt))^2.63)))/
    ((((Tmax - .super$pars$reftemp)/(Tmax - Topt))^0.2) * exp((0.2/2.63) * (1 - ((Tmax - .super$pars$reftemp)/(Tmax - Topt))^2.63)))

}

f_tcor_century2 <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.Century2
  Tmax <- 45
  Topt <- 35
  (3.439423 * exp((0.2/2.63) * (1 - (((Tmax - .super$env$temp)/(Tmax - Topt))^2.63)) * ((Tmax - .super$env$temp)/(Tmax - Topt))^0.2))/
    (3.439423 * exp((0.2/2.63) * (1 - (((Tmax - .super$pars$reftemp)/(Tmax - Topt))^2.63)) * ((Tmax - .super$pars$reftemp)/(Tmax - Topt))^0.2))
}

f_tcor_daycent1 <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.Daycent1
  (0.08 * exp(0.095 * .super$env$temp))/
  (0.08 * exp(0.095 * .super$pars$reftemp))
}

f_tcor_daycent2 <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.Daycent2
  (0.56 + (1.46 * atan(pi * 0.0309 * (.super$env$temp - 15.7)))/pi)/
    (0.56 + (1.46 * atan(pi * 0.0309 * (.super$pars$reftemp - 15.7)))/pi)
}

f_tcor_daycent2_abramoff <- function(.,C,t,i){
  #daycent2 function as implemented in millennial model (i.e. different parameters)
  t1 = 15.4
  t2 = 11.75
  t3 = 29.7
  t4 = 0.031
  (t2 + (t3/pi)* atan(pi*t4*(.super$env$temp - t1))) / (t2 + (t3/pi)* atan(pi*t4*(.super$pars$reftemp - t1)))
}

f_tcor_q10 <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.Q10
  k_ref <- 1
  k_ref * .super$pars$q10^((.super$env$temp - .super$pars$reftemp)/10)
}

f_tcor_rothc <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.RothC
  (47.9/(1 + exp(106/(ifelse(.super$env$temp >= -18.3, .super$env$temp, NA) + 18.3))))/
  (47.9/(1 + exp(106/(ifelse(.super$pars$reftemp >= -18.3, .super$pars$reftemp, NA) + 18.3))))
}

f_tcor_lin_adair <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.linear
  a = 0.198306
  b = 0.036337
  (a + b * .super$env$temp)/ (a + b * .super$pars$reftemp)
}

f_tcor_landt <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.LandT
  (exp(308.56 * ((1/56.02) - (1/((.super$env$temp + 273) - 227.13)))))/
    (exp(308.56 * ((1/56.02) - (1/((.super$pars$reftemp + 273) - 227.13)))))
}

f_tcor_kb <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.KB
  (exp(-3.764 + 0.204 * .super$env$temp * (1 - 0.5 * .super$env$temp/36.9)))/
    (exp(-3.764 + 0.204 * .super$pars$reftemp * (1 - 0.5 * .super$pars$reftemp/36.9)))
}

f_tcor_demeter <- function(.,C,t,i){
  #sourced from SoilR
  #SoilR::fT.Demeter
  #this func is 1 at 20 without scaling
  (exp((log(.super$pars$q10)/10) * (.super$env$temp - 20)))
}

f_tcor_standcarb <- function(.,C,t,i){
  #sourced from SoilR and standardized to 20C
  #SoilR::fT.Standcarb
  Topt = 45
  Tlag = 4
  Tshape = 15
  (exp(-1 * (.super$env$temp/(Topt + Tlag))^Tshape) * .super$pars$q10^((.super$env$temp - 10)/10))/
    (exp(-1 * (.super$pars$reftemp/(Topt + Tlag))^Tshape) * .super$pars$q10^((.super$pars$reftemp - 10)/10))
}

f_tcor_arrhenius <- function(.,C,t,i) {
  #sourced from MAAT leaf model
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # Arrhenius equation
  
  # input parameters  
  # Ea     -- rate of increase to optimum  (J mol-1)
  # R      -- molar gas constant J mol-1 K-1
  
  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  #convert to Kelvin
  Trk <- .super$pars$reftemp + 273.15
  Tsk <- .super$env$temp + 273.15
  
  exp( .super$pars$ea[[i]]*(Tsk-Trk) / (.super$pars$R*Tsk*Trk) )
}

