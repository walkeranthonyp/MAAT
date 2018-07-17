################################
#
# template for MAAT object system functions
# 
# AWalker February 2018
#
################################



################################
# Template system function one
# - this is a variant on the hello world program

f_templatesys_1 <- function(.) {

  # define state 
  .$state$text    <- get(.$fnames$text)(.)
  .$state$calcval <- get(.$fnames$calcval)(.)

  # call print function
  get(.$fnames$print)(.)  
  
}



### END ###
