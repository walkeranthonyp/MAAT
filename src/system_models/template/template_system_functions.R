################################
#
# template for MAAT object system functions
# 
# AWalker February 2018
#
################################



################################
# Template system function one
f_templatesys_1 <- function(.) {

  # calculate state parameters
  # photosynthetic parameters
  .$state_pars$vcmax   <- get(.$fnames$vcmax)(.)
  .$state_pars$jmax    <- get(.$fnames$jmax)(.)
  .$state_pars$tpu     <- get(.$fnames$tpu)(.)
  .$state_pars$rd      <- get(.$fnames$respiration)(.)
  .$state_pars$alpha   <- 0.5 * (1-.$pars$f)

  .$state$respiration  <- get(.$fnames$rl_rd_scalar)(.) * .$state$respiration
  # determine rate limiting step - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  .$state$A       <- get(.$fnames$solver)(.)      
  # assign the limitation state - assumes the minimum is the dominant limiting rate
  .$state$lim     <- c('wc','wj','wp')[which(c(.$state$wc,.$state$wj,.$state$wp)==min(c(.$state$wc,.$state$wj,.$state$wp),na.rm=T))]       

}
