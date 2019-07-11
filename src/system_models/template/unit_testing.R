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
template_object$.test()
template_object$.test(verbose=T)
template_object$.test_change_func()
template_object$.test_change_func(template.text='f_text_one')
template_object$.test_change_func(template.text='f_text_one',template.print_out='f_print_out_textandval')
template_object$.test_change_func(template.text='f_text_one',template.print_out='f_print_out_textandval',template.calcval='f_calcval_power')


template_object$fnames
template_object$fns
template_object$state
template_object$state_pars
template_object$pars
template_object$env



### END ###