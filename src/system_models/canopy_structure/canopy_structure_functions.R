################################
#
# Canopy structure functions
# 
# AWalker October 2018
#
################################



# LAI
################################

f_lai_constant <- function(.) .super$pars$lai

f_lai_env      <- function(.) .super$env$lai

f_lai_leafdem_env <- function(.) {
   .super$state$lai_young  <- .super$env$lai_young 
   .super$state$lai_mature <- .super$env$lai_mature
   .super$state$lai_old    <- .super$env$lai_old

   .super$state$lai_young + .super$state$lai_mature + .super$state$lai_old
}

# calculate sphagnum 'lai' as a logistic function of water table depth
f_lai_sphagnum <- function(.) {
  b <- 0.1
    
  .super$pars$lai_max*b / (b + exp(.super$pars$lai_curve*(.super$env$water_td-.super$env$sphag_h)/10) )
}


# Other stuff 
################################

f_leafdem_vcmax0_constant <- function(.) unlist(.super$pars$vcmax0)


f_leafdem_upper_wu2017 <- function(.) {
  prop_young  <- .super$state$lai_young  * .super$pars$lai_ftop / .super$pars$lai_upper
  prop_mature <- .super$state$lai_mature * .super$pars$lai_ftop / .super$pars$lai_upper

  c(prop_young, prop_mature, 1 - prop_young - prop_mature )
}


f_leafdem_lower_wu2017 <- function(.) {
  prop_young  <- .super$state$lai_young  * (1-.super$pars$lai_ftop) / (.super$state$lai-.super$pars$lai_upper)
  prop_mature <- .super$state$lai_mature * (1-.super$pars$lai_ftop) / (.super$state$lai-.super$pars$lai_upper)

  c(prop_young, prop_mature, 1 - prop_young - prop_mature )
}


# Canopy/Leaf water status  
################################

f_water_status_none <- function(.) .super$leaf$state$fwdw_ratio <- NA

# set sphagnum water status
f_water_status_sphagnum <- function(.) {
  .super$leaf$state$fwdw_ratio <- .$fns$fwdw()
}


# Sphagnum functions
################################

# Calculates Sphagnum fresh weight : dry weight ratio as a function of water level (mm) - linear
f_fwdw_wtd_lin <- function(.) {
  
  if((.super$env$water_td - .super$env$sphag_h) > 0) .super$pars$fwdw_wl_sat 
  else .super$pars$fwdw_wl_sat + .super$pars$fwdw_wl_slope * -(.super$env$water_td - .super$env$sphag_h) 
}

# Calculates Sphagnum fresh weight : dry weight ratio as a function of water level (mm) - exponential
f_fwdw_wtd_exp <- function(.) {
  # Strack & Price 2009
  
  if((.super$env$water_td - .super$env$sphag_h) > 0) exp( .super$pars$fwdw_wl_exp_b + .super$pars$fwdw_wl_exp_a*0 )
  else exp( .super$pars$fwdw_wl_exp_b + .super$pars$fwdw_wl_exp_a * (-(.super$env$water_td - .super$env$sphag_h)/10) ) 
}



### END ###
