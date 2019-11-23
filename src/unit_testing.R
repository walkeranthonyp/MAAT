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
out <- wrapper_object$.test_simple(gen_metd=T)
wrapper_object$build
wrapper_object$run0
wrapper_object$run1
wrapper_object$run2
wrapper_object$run3
wrapper_object$run4

wrapper_object$model$fns$sys

wrapper_object$dataf
wrapper_object$dynamic


source('wrapper_object.R')
out <- wrapper_object$.test(mc=F, oconf=F, gen_metd=F )
out <- wrapper_object$.test(mc=F, oconf=F )
out <- wrapper_object$.test(mc=T, oconf=F )

library(lattice)
xyplot(A~leaf.par  | leaf.etrans*leaf.rs, out, groups=leaf.temp,type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5)); panel.xyplot(...)} )
xyplot(A~leaf.par  | leaf.etrans*leaf.rs, out, groups=leaf.vpd, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5)); panel.xyplot(...)} )
xyplot(A~leaf.par  | leaf.etrans*leaf.rs, out, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5)); panel.xyplot(...)} )
xyplot(ci~leaf.par | leaf.etrans*leaf.rs, out, type='l', panel=function(...) { panel.abline(h=seq(0,20,2.5)); panel.xyplot(...)} )


source('wrapper_object.R')
out <- wrapper_object$.test_con(mc=F, oconf=F, metd=NULL ) 
out <- wrapper_object$.test_con(mc=T, oconf=F )
out <- wrapper_object$.test_con(mc=F, oconf=F )
out <- wrapper_object$.test_con(mc=F, oconf=F,
                                  sfnames=list(leaf.tcor_asc.vcmax = 'f_tcor_asc_Arrhenius',
                                               leaf.tcor_asc.jmax  = 'f_tcor_asc_Arrhenius'
                                  ),
                                  dpars=list(leaf.Ha.vcmax=c(4e4,7e4),leaf.Ha.jmax=c(1e4,7e4))
                                )


source('wrapper_object.R')
out <- wrapper_object$.test_init()
wrapper_object$init_static
class(wrapper_object$static)
wrapper_object$static

source('wrapper_object.R')
out <- wrapper_object$.test_init_npools()


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
head(wrapper_object$dataf$out_saltelli)
head(wrapper_object$dataf$out)
wrapper_object$dataf$out_saltelli
wrapper_object$init_output_matrix
wrapper_object$init_output_matrix()
wrapper_object$write_output
wrapper_object$print_saltelli


source('wrapper_object.R')
out <- wrapper_object$.test_saltelli(mc=T, pr=6, oconf=F, n=10 )
wrapper_object$dataf$env
wrapper_object$dataf$fnames
wrapper_object$dataf$pars


# MCMC Mixture test
source('wrapper_object.R')
out <- wrapper_object$.test_mcmc_mixture(mcmc_maxiter=1000)
out <- wrapper_object$.test_mcmc_mixture(mcmc_type='dream',
                                         mcmc_maxiter=5000, mcmc_chains=8, mc=F, pr=4,
                                         mu_vector=c(-8, 0, 8), sd_vector=c(1, 1, 1),
                                         height_vector=c(0.1, 0.3, 0.6), mixture_scale=1e12
                                         )
names(out)
wrapper_object$dynamic
wrapper_object$dynamic$pars_eval
wrapper_object$dataf$pars
wrapper_object$model$pars
wrapper_object$dataf$pars_array[,,1]
wrapper_object$dataf$pars_lklihood
wrapper_object$dataf$pars_array
wrapper_object$wpars

out <- wrapper_object$.test_mcmc_mixture(mcmc_maxiter=20)
out <- wrapper_object$.test_mcmc_mixture(mcmc_maxiter=20, mc=T )

out <- wrapper_object$.test_mcmc_mixture(mcmc_maxiter=1000)
out <- wrapper_object$.test_mcmc_mixture(mcmc_maxiter=1000, mc=T )

wrapper_object$dataf$pars
wrapper_object$dataf$out
wrapper_object$dataf$pars_lklihood[,1]
wrapper_object$dataf$pars_array[,,1]
wrapper_object$dataf$pars_lklihood
wrapper_object$dataf$pars_array



# MCMC linear regression test
source('wrapper_object.R')
out <- wrapper_object$.test_mcmc_linreg()
out <- wrapper_object$.test_mcmc_linreg(cverbose=T)
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=10)
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=1000)
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=1000,
                                        mcmc_test.a='runif(n,-20,20)',
                                        mcmc_test.b='runif(n,-20,20)')
names(out)
plot(out[[3]])
wrapper_object$wpars
wrapper_object$mcmc
wrapper_object$dynamic

source('wrapper_object.R')
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=2, mc=F )
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=2, mc=T )

out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=1000, mc=F )
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=1000, mc=T )

wrapper_object$dataf$pars
wrapper_object$dataf$out
wrapper_object$dataf$pars_lklihood[,1]
wrapper_object$dataf$pars_array[,,1]
wrapper_object$dataf$pars_lklihood
wrapper_object$dataf$pars_array

wrapper_object$dataf$pars_array[,,10]
wrapper_object$dataf$out
wrapper_object$dataf$out_mcmc[,,1]
wrapper_object$dataf$out_mcmc[,,10]

wrapper_object$dataf$met
wrapper_object$dataf$obs
wrapper_object$dataf$obsse
wrapper_object$dataf$out

xyplot(wrapper_object$dataf$obs~wrapper_object$dataf$met)
lm(wrapper_object$dataf$obs~wrapper_object$dataf$met)

wrapper_object$mcmc$boundary_min
wrapper_object$mcmc$boundary_max


# Canopy
source('wrapper_object.R')
out <- wrapper_object$.test_can(mc=F,verbose=F)
head(out[[1]])
print(out[[2]])

wrapper_object$model$state
wrapper_object$model$leaf$fns$sys
wrapper_object$dataf$met



### END ###