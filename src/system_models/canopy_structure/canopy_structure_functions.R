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

f_lai_sphagnum <- function(.) {
  # calculate sphagnum 'lai' as a logistic function of water table depth
  b <- 0.1
    
  .$pars$lai_max*b / (b + exp(.$pars$lai_curve*(.$env$water_td-.$env$sphag_h)/10) )
}


# Other stuff 
################################

f_leafdem_vcmax0_constant <- function(.) unlist(.$pars$vcmax0)

f_leafdem_upper_env <- function(.) c(.$env$upper_can_prop_young, .$env$upper_can_prop_mature, .$env$upper_can_prop_old)
f_leafdem_lower_env <- function(.) c(.$env$lower_can_prop_young, .$env$lower_can_prop_mature, .$env$lower_can_prop_old)



### END ###
