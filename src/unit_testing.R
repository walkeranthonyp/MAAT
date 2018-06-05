###########################
#
# Unit testing of the MAAT wrapper
#
# AWalker
# May 2015
#
###########################



### Load model scripts 
###############################


# Wrapper
# Leaf
rm(list=ls())

source('wrapper_object.R')
out <- wrapper_object$.test_simple()
out <- wrapper_object$.test_simple(metd=T)

source('wrapper_object.R')
out <- wrapper_object$.test(mc=F,metd=F,oconf=F)
out <- wrapper_object$.test(mc=T,oconf=F)
out <- wrapper_object$.test(mc=F,oconf=F)

class(out)
head(out)

library(lattice)
xyplot(A~leaf.par |leaf.etrans*leaf.rs, out, groups=leaf.temp,type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs, out, groups=leaf.vpd, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs, out, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(ci~leaf.par |leaf.etrans*leaf.rs, out, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})

# test model mimic
source('wrapper_object.R')
out <- wrapper_object$.test_mimic()
out <- wrapper_object$.test_mimic(metd=T)
wrapper_object$model$output()

wrapper_object$model$env
wrapper_object$model$pars
wrapper_object$model$state_pars
wrapper_object$model$state
wrapper_object$model$fnames


wrapper_object$init_static
wrapper_object$init_dynamic


wrapper_object$dataf
wrapper_object$dataf$met

f_tpu_lin(wrapper_object$model)
f_tpu_tcor_independent(wrapper_object$model)
f_rs_ball1987(wrapper_object$model)
f_rs_ball1987_fe(wrapper_object$model)
f_rs_ball1987(wrapper_object$model,10,50)
f_R_Brent_solver(wrapper_object$model)
f_A_r_leaf(wrapper_object$model,100)
f_rd_lin_vcmax(wrapper_object$model)


# Leaf, Ye SA
source('wrapper_object.R')
out <- wrapper_object$.test_ye(mc=F,pr=6,oconf=F,n=10)
out <- wrapper_object$.test_ye(mc=T,pr=6,oconf=F,n=10)


# Leaf, Saltelli SA
source('wrapper_object.R')
out <- wrapper_object$.test_saltelli(mc=F, pr=6, oconf=F, n=10 )

source('wrapper_object.R')
out <- wrapper_object$.test_saltelli(mc=T, pr=6, oconf=F, n=10 )
wrapper_object$dataf$env
wrapper_object$dataf$fnames
wrapper_object$dataf$pars


# Canopy
source('wrapper_object.R')
out <- wrapper_object$.test_can(mc=F,verbose=F)
head(out[[1]])
print(out[[2]])

wrapper_object$model$state
wrapper_object$dataf$met



### END ###
