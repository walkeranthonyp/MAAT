################################
#
# Canopy structure functions
# 
# AWalker October 2018
#
################################



# LAI
################################

f_lai_constant <- function(.) .$pars$lai

f_lai_env      <- function(.) .$env$lai

f_lai_leafdem_env <- function(.) {
   .$state$lai_young  <- .$env$lai_young 
   .$state$lai_mature <- .$env$lai_mature
   .$state$lai_old    <- .$env$lai_old

   .$state$lai_young + .$state$lai_mature + .$state$lai_old
}

# calculate sphagnum 'lai' as a logistic function of water table depth
f_lai_sphagnum <- function(.) {
  b <- 0.1
    
  .$pars$lai_max*b / (b + exp(.$pars$lai_curve*(.$env$water_td-.$env$sphag_h)/10) )
}


# Other stuff 
################################

f_leafdem_vcmax0_constant <- function(.) unlist(.$pars$vcmax0)


f_leafdem_upper_wu2017 <- function(.) {
  prop_young  <- .$state$lai_young  * .$pars$lai_ftop / .$pars$lai_upper
  prop_mature <- .$state$lai_mature * .$pars$lai_ftop / .$pars$lai_upper

  c(prop_young, prop_mature, 1 - prop_young - prop_mature )
}


f_leafdem_lower_wu2017 <- function(.) {
  prop_young  <- .$state$lai_young  * (1-.$pars$lai_ftop) / (.$state$lai-.$pars$lai_upper)
  prop_mature <- .$state$lai_mature * (1-.$pars$lai_ftop) / (.$state$lai-.$pars$lai_upper)

  c(prop_young, prop_mature, 1 - prop_young - prop_mature )
}




### END ###
