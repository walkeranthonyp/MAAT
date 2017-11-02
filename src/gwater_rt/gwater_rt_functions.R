################################
#
# gwater_rt for MAAT object functions
# 
# AWalker March 2017
#
################################


f_none <- function(.) {
  NA
}

quad_sol <- function(a,b,c,out='lower') {
  # robust numerical solution to the quadratic
  # taken from Numerical Recipes

  q     <- -0.5 * ( b + sign(b)*(b^2 - 4*a*c)^0.5 )
  roots <- c( q/a , c/q )
  
  if(out=='lower')      min(roots,na.rm=T) 
  else if(out=='upper') max(roots,na.rm=T)
  else roots 
}





### PHOTOSYNTHESIS FUNCTIONS
################################

# Carboxylation limitation
f_wc_farquhar1980 <- function(.,cc=.$state$cc){   
  # Farquhar 1980 eq to calculate RuBisCO limited photosynthetic rate
  # umol m-2 s-1
  # Wc<-(Vcmax*(Ci-Gamma))/(Ci+Kc*(1+(O/Ko)))
  
  # dvars
  # Ci (or Cc but vcmax needs adjusting)   
  
  # params
  # Vcmax maximum carboxylation rate of RuBisCO  (umol m-2 s-1
  # Ko Michaelis constant for O2 mbar            (kPa)
  # Kc Michaelis constant for CO2                ( Pa) 
  
  # Ko and O must be in same units
  # Kc, Ci and Gamma must be in same units
  
  #gross of rd and G*
#   .$state_pars$vcmaxlt*cc /
#     (cc+.$state_pars$Kc*(1+(.$state$oi/.$state_pars$Ko)))
  .$state_pars$vcmaxlt /
    (cc + .$state_pars$Km)
}


# electron transport limitation
f_wj_generic <- function(.,cc=.$state$cc){
  # generic eq to calculate light limited photosynthetic rate from user selected electron transport and etoc function
  # umol m-2 s-1
  
  # calculate gross electron transport limited carboxylation 
  # i.e. rd and G* not included
  #   (.$state$J/4) * (cc/(cc+2*.$state_pars$gstar)) # currently no other formualtion exists (other than slighly different parameters in the denominator)
  .$state$J / (4*(cc+2*.$state_pars$gstar))        # currently no other formualtion exists (other than slighly different parameters in the denominator)
}

f_j_farquhar1980 <- function(.){
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per incident photon
  # Farquhar 1980 and others
  # this is the replacement for eq A2 added in the proofing stage (very last line of the manuscript after all refs etc)
  
  .$state_pars$jmaxlt * .$env$par * (1 - .$pars$f)/(.$env$par + 2.2*.$state_pars$jmaxlt)
}

f_j_harley1992 <- function(.){
  # Harley etal 1992 eq to calculate electron transport rate
  # umol e m-2 s-1
  harley_alpha <- .$pars$a * .$state_pars$alpha
  harley_alpha*.$env$par / ( 1.0 + (harley_alpha**2)*(.$env$par**2) / (.$state_pars$jmaxlt**2) )**0.5  
}


















