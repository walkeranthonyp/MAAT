################################
#
# template for MAAT object functions
# 
# AWalker March 2017
#
################################


### FUNCTIONS
################################

# generate text 
f_text_combine <- function(.) {   
  paste(.super$pars$text1, .super$pars$text2)
}

f_text_one <- function(.) {   
  .super$pars$text3
}

# calculate a value
f_calcval_product <- function(.) {
  .super$pars$val1 * .super$pars$val2
}

f_calcval_power <- function(.) {
  .super$pars$val1 ^ .super$pars$val2
}

# print values
f_print_textonly <- function(.) {
  print(.super$state$text, quote=F )
}

f_print_textandval <- function(.) {
  print(paste(.super$state$text,.super$state$calcval), quote=F )
}



### END ###
