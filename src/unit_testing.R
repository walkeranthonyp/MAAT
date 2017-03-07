###########################
#
# Unit testing of the MAAT models and wrapper
#
# AWalker
# May 2015
#
###########################



### Load model scripts 
###############################



# Leaf MAAT
rm(list=ls())
library(lattice)
source('leaf_object.R')

# leaf
leaf_object$.test_leaf(verbose=F,verbose_loop=F)
leaf_object$fnames
leaf_object$state
leaf_object$state_pars
leaf_object$pars
leaf_object$env

leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=1100)
leaf_object$.test_aci(leaf.ca_conc=seq(0.1,2000,50))
leaf_object$.test_aci(diag=T)
leaf_object$.test_aci(diag=T,leaf.par=c(450,900) )
leaf_object$.test_aci(diag=T,leaf.par=c(450,900,1200,1500) )
leaf_object$.test_aci_light()
leaf_object$.test_aci_light(diag=T)

# out <- leaf_object$.test_aci_analytical()
out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011')
out <- leaf_object$.test_aci_analytical(rs='f_rs_ball1987')
out <- leaf_object$.test_aci_analytical(rs='f_rs_leuning1995')
out <- leaf_object$.test_aci_analytical(rs='f_rs_constantCiCa')
out <- leaf_object$.test_aci_analytical(rs='f_rs_cox1998')


# Wrapper
# Leaf
rm(list=ls())
source('wrapper_object.R')
library(lattice)
out <- wrapper_object$.test(mc=T,oconf=F)

head(out[[1]])
length(out[[1]][,1])

xyplot(A~leaf.par |leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs,out[[1]],groups=leaf.vpd, type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(cc~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(gs~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,1,0.02));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs,out[[1]],type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})


# Leaf, Ye SA
rm(list=ls())
source('wrapper_object.R')
library(lattice)
out <- wrapper_object$.test_ye(mc=T,pr=6,oconf=F,n=20)

out
out[[2]]


# Leaf, Saltelli SA
rm(list=ls())
source('wrapper_object.R')
library(lattice)
out <- wrapper_object$.test_saltelli(mc=T,pr=6,oconf=F,n=10)

out
out[[1]]
out[[2]]
dim(out[[1]])
length(out[[2]])
dim(out[[2]][[1]][[1]])
wrapper_object$dataf$env
wrapper_object$dataf$fnames
wrapper_object$dataf$pars

wrapper_object$model$env
wrapper_object$model$fnames
wrapper_object$model$pars










