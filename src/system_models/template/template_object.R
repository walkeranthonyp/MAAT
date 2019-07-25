################################
#
# MAAT template model object 
#
# AWalker November 2017
#
################################

library(proto)
source('template_functions.R')
source('template_system_functions.R')



# TEMPLATE OBJECT
###############################################################################

# use generic template
setwd('..')
source('generic_model_object.R')
template_object <- as.proto( system_model_object$as.list() )
rm(system_model_object)
setwd('template')



# assign object functions
###########################################################################
template_object$name <- 'template'
# if this new object will contain nested child objects uncomment and edit the below lines of code and delete this comment 
# - otherwise delete all 5 lines
#template_object$child_list      <- list('child_name1') 
#template_object$build_child     <- build_child  
#template_object$configure_child <- configure_child  



# assign unique run function
###########################################################################
# if run function needs to be modified - add new function here



# functions unique to object that do not live in fnames/fns, i.e. do not vary ever
###########################################################################
# add structural functions (i.e not alternative process functions)  unique to model object here



# assign object variables 
###########################################################################

# function names
####################################
template_object$fnames <- list(
  sys       = 'f_sys_1',
  text      = 'f_text_combine',
  calcval   = 'f_calcval_product',
  print_out = 'f_print_textonly'
)


# environment
####################################
template_object$env <- list(
  text1 = 'hello',    
  text2 = 'world',    
  text3 = 'The answer is:'    
)


# state
####################################
template_object$state <- list(
  text    = character(0),     # text to print
  calcval = numeric(0)        # value to print
)


# state parameters (i.e. calculated parameters)
####################################
template_object$state_pars <- list(
  vcmax    = numeric(0)   
)


# parameters
####################################
template_object$pars   <- list(
  val1  = 6,           
  val2  = 7           
)


# run control parameters
####################################
template_object$cpars <- list(
  verbose  = F,          # write diagnostic output during runtime 
  cverbose = F,          # write diagnostic output from configure function 
  output   = 'run'       # type of output from run function
)



# output functions
#######################################################################        

f_object_template_run <- function(.) {
  unlist(.$state)
}

f_object_template_state <- function(.) {
  unlist(.$state)
}

f_object_template_full <- function(.) {
  c(unlist(.$state),unlist(.$statei_pars))
}



# test functions
#######################################################################        

template_object$.test <- function(., verbose=F) {
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))

  .$run()
}


template_object$.test_change_func <- function(., verbose=F,
                                              template.text='f_text_combine',
                                              template.calcval='f_calcval_product',
                                              template.print_out='f_print_out_textonly' ) {
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))
  
  .$fnames$text      <- template.text
  .$fnames$calcval   <- template.calcval
  .$fnames$print_out <- template.print_out

  # configure_test must be called if any variables in the fnames lista re reassigned
  .$configure_test()  
  .$run()
}


template_object$.test_change_pars <- function(., verbose=F,
                                              template.text1='hello',
                                              template.text2='world' ) {
  if(verbose) str(.)
  .$build(switches=c(F,verbose,F))

  .$pars$text1 <- template.text1
  .$run()
}



### END ###
