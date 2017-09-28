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

head(out[[1]])
length(out[[1]][,1])

library(lattice)
xyplot(A~leaf.par |leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs,out[[1]],groups=leaf.vpd, type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(cc~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(gs~leaf.par|leaf.etrans*leaf.rs,out[[1]],groups=leaf.temp,type='l',panel=function(...) { panel.abline(h=seq(0,1,0.02));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs,out[[1]],type='l',panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})


# Leaf, Ye SA
source('wrapper_object.R')
out <- wrapper_object$.test_ye(mc=T,pr=6,oconf=F,n=100)

out[[2]]
out <- out[[1]]
length(out)
names(out)
class(out)
class(out$fnames)
class(out$out)
class(out$env)
colnames(out$env)
out$env[,which(colnames(out$env)=='leaf.ca_conc')]


# Leaf, Saltelli SA
rm(list=ls())
source('wrapper_object.R')
out <- wrapper_object$.test_saltelli(mc=F,pr=6,oconf=F,n=10)
out <- wrapper_object$.test_saltelli(mc=T,pr=6,oconf=F,n=10)

out
length(out)
dim(out[[1]])
dim(out[[2]])

out[[1]][1,1,,]
out[[2]][1,1,,,]

wrapper_object$dataf$env
wrapper_object$dataf$fnames
wrapper_object$dataf$pars










