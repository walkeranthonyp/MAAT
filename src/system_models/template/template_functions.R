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
  paste(.$pars$text1, .$pars$text2)
}

f_text_one <- function(.) {   
  .$pars$text3
}

# calculate a value
f_calcval_product <- function(.) {
  .$pars$val1 * .$pars$val2
}

f_calcval_power <- function(.) {
  .$pars$val1 ^ .$pars$val2
}

# print values
f_print_textonly <- function(.) {
  print(.$state$text)
}

f_print_textandval <- function(.) {
  print(paste(.$state$text,.$state$calcval))
}



### END ###
