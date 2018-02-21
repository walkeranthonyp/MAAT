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

  # print combined text and calculated state
  get(.$fnames$print)(.)  
  
}
