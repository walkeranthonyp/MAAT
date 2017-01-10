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

leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=1100)
leaf_object$.test_aci()
leaf_object$.test_aci_light()



# Wrapper
# Leaf
rm(list=ls())
source('wrapper_object.R')
library(lattice)
out <- wrapper_object$.test(mc=T,oconf=F)

head(out[[1]])
length(out[[1]][,1])

xyplot(A~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.vpd,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(cc~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(gs~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,1,0.02));panel.xyplot(...)})
xyplot(A~leaf.par|leaf.etrans*leaf.rs,out[[1]],type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})


# Leaf, Ye SA
rm(list=ls())
source('wrapper_object.R')
library(lattice)
out <- wrapper_object$.test_ye(mc=F,oconf=F)

out
out[[1]]














