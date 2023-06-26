################################
#
# MAAT template system representation functions (SRFs) 
# 
# AWalker February 2018
#
################################



################################
# Template system function one
# - this is a variant on the hello world program

f_sys_1 <- function(.) {

  # define state 
  .super$state$text    <- .$text()
  .super$state$calcval <- .$calcval()

  # call print function
  .$print_out()  
  
}



### END ###
