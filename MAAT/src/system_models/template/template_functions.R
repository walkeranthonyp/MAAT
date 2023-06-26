################################
#
# MAAT template process representation functions (PRFs)
# 
# AWalker March 2017
#
################################


### FUNCTIONS
################################

# generate text 
f_text_combine <- function(.) {   
  paste(.super$env$text1, .super$env$text2)
}

f_text_one <- function(.) {   
  .super$env$text3
}

# calculate a value
f_calcval_product <- function(.) {
  .super$pars$val1 * .super$pars$val2
}

f_calcval_power <- function(.) {
  .super$pars$val1 ^ .super$pars$val2
}

# print values
f_print_out_textonly <- function(.) {
  print(.super$state$text, quote=F )
  #print(.super$state$text)
}

f_print_out_textandval <- function(.) {
  print(paste(.super$state$text,.super$state$calcval), quote=F )
}



### END ###
