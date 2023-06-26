################################
#
# General mathematical functions
# 
# AWalker February 2018
#
################################



# robust numerical solution to the quadratic
# taken from Numerical Recipes
quad_sol <- function(a,b,c,out='lower') {

  q     <- -0.5 * ( b + sign(b)*(b^2 - 4*a*c)^0.5 )
  roots <- c( q/a , c/q )
  
  if(out=='lower')      min(roots,na.rm=T) 
  else if(out=='upper') max(roots,na.rm=T)
  else roots 
}

# function that mimics some of the behaviour of the matlab 'diag' function
diag_m <- function(v,k=0) {
  if(k==0) diag(v)
  else if(k>0) {
    rbind(cbind(matrix(0,length(v),k),diag(v)), matrix(0,k,length(v)+k) )
  } else {
    k <- abs(k)
    rbind(matrix(0,k,length(v)+k), cbind(diag(v), matrix(0,length(v),k) ))
  }  
}



### END ###
