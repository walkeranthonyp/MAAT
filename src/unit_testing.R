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
if (! file.exists('~/tmp')) dir.create('~/tmp',recursive=TRUE)

source('wrapper_object.R')
out <- wrapper_object$.test_simple()
out <- wrapper_object$.test_simple(metd=T)

source('wrapper_object.R')
out <- wrapper_object$.test(mc=F, metd=F, oconf=F )
out <- wrapper_object$.test(mc=T, oconf=F )
out <- wrapper_object$.test(mc=F, oconf=F )

library(lattice)
xyplot(A~leaf.par |leaf.etrans*leaf.rs, out, groups=leaf.temp,type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs, out, groups=leaf.vpd, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(A~leaf.par |leaf.etrans*leaf.rs, out, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})
xyplot(ci~leaf.par |leaf.etrans*leaf.rs, out, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5));panel.xyplot(...)})

source('wrapper_object.R')
out <- wrapper_object$.test_con(mc=F, metd=F, oconf=F )
out <- wrapper_object$.test_con(mc=T, oconf=F )
out <- wrapper_object$.test_con(mc=F, oconf=F )
out <- wrapper_object$.test_con(mc=F, oconf=F,
                                  sfnames=list(leaf.vcmax_tcor_asc      = 'f_temp_scalar_Arrhenius',
                                               leaf.jmax_tcor_asc       = 'f_temp_scalar_Arrhenius'
                                  ),
                                  dpars=list(leaf.Ha.vcmax=c(4e4,7e4),leaf.Ha.jmax=c(1e4,7e4))
                                )

source('wrapper_object.R')
out <- wrapper_object$.test_init()
wrapper_object$init_static
wrapper_object$static

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


# MCMC Mixture test
source('wrapper_object.R')
out <- wrapper_object$.test_mcmc_mixture(mc=F, pr=6 )
out <- wrapper_object$.test_mcmc_mixture(mc=F, pr=6, mcmc_chains=8, mcmc_maxiter=10000 )
out[[3]]
df1 <- data.frame(lklihood=as.vector(t(out[[2]])), chain=rep(1:dim(out[[2]])[1],each=dim(out[[2]])[2]) )
xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) , df1, groups=chain, auto.key=T, type='l' )
xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) | chain , df1, auto.key=T, type='l' )

wrapper_object$dynamic
wrapper_object$dataf$pars
wrapper_object$model$pars
wrapper_object$dataf$pars_array
wrapper_object$dataf$pars_lklihood
wrapper_object$wpars


# MCMC linear regression test
source('wrapper_object.R')
out <- wrapper_object$.test_mcmc_linreg()
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=5000)
out[[3]]
dim(out[[2]])
df1 <- data.frame(lklihood=as.vector(t(out[[2]])), chain=rep(1:dim(out[[2]])[1],each=dim(out[[2]])[2]) )
xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) , df1, groups=chain, auto.key=T, type='l' )
xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) | chain , df1, auto.key=T, type='l' )


wrapper_object$dynamic
wrapper_object$dataf$pars
wrapper_object$model$pars
wrapper_object$dataf$pars_array
wrapper_object$dataf$pars_lklihood
wrapper_object$wpars


# Canopy
source('wrapper_object.R')
out <- wrapper_object$.test_can(mc=F,verbose=F)
head(out[[1]])
print(out[[2]])

wrapper_object$model$state
wrapper_object$dataf$met



### END ###