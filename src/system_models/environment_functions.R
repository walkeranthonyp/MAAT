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
  rh <- (1 - .super$env$vpd / .$sat_vp())
  .super$env$rh <- min(max(rh,0),1)
}

f_vpd_from_rh <- function(.) {
  vpd <- (1 - .super$env$rh) * .$sat_vp()
  .super$env$vpd <- max(vpd,0)
}

f_sat_vp_allen1998 <- function(., t=.super$state$leaf_temp ) {
  # calculate saturation vapour pressure - Allen 1998
  # in kPa
  # t -- oC -- air temp
  
  0.6108 * exp( 17.27*t / (t+237.3) )
}

f_sat_vp_buck1981 <- function(., t=.super$state$leaf_temp ) {
  # calculate saturation vapour pressure - Buck 1981 J. App. Met.
  # in kPa
  # t -- oC -- air temp
  
  0.61121 * exp( 17.502*t / (240.97 + t) )
}



### END ###
