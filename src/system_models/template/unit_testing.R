###########################
#
# template for MAAT object unit testing
# 
# AWalker March 2017
#
###########################



### Load model scripts 
###############################

# template
source('template_object.R')
template_object$.test(verbose=F,verbose_loop=F)
template_object$.test_change_func(verbose=F,verbose_loop=F)
template_object$.test_change_func(verbose=F,verbose_loop=F,template.text='f_text_one')
template_object$.test_change_func(verbose=F,verbose_loop=F,template.text='f_text_one',template.print='f_print_textandval')
template_object$.test_change_func(verbose=F,verbose_loop=F,template.text='f_text_one',template.print='f_print_textandval',template.calcval='f_calcval_power')


template_object$fnames
template_object$state
template_object$state_pars
template_object$pars
template_object$env



### END ###
