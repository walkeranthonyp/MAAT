################################
#
# General environmental functions
# 
# AWalker February 2018
#
################################



# VPD
################################
f_rh_from_vpd <- function(.) {
  rh <- (1 - .$env$vpd / f_sat_vp_allen(.))
  min(max(rh,0),1)
}

f_sat_vp_allen <- function(.) {
  # calculate saturation vapour pressure - Allen 1998
  # in kPa
  # t -- oC -- air temp
  
  0.6108 * exp( 17.27*.$state$leaf_temp /(.$state$leaf_temp+237.3) )
}

f_sat_vp_buck <- function(t){
  # calculate saturation vapour pressure - Buck 1981 J. App. Met.
  # in kPa
  # t -- oC -- air temp
  
  0.61121 * exp( 17.502*t / (240.97 + t) )
}



