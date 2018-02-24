################################
#
# General mathematical functions
# 
# AWalker February 2018
#
################################



# robust numerical solution to the quadratic
quad_sol <- function(a,b,c,out='lower') {
  # taken from Numerical Recipes

  q     <- -0.5 * ( b + sign(b)*(b^2 - 4*a*c)^0.5 )
  roots <- c( q/a , c/q )
  
  if(out=='lower')      min(roots,na.rm=T) 
  else if(out=='upper') max(roots,na.rm=T)
  else roots 
}




