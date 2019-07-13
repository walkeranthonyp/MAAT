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
wrapper_object$model$fns$sys
wrapper_object$run0
wrapper_object$run1
wrapper_object$run2
wrapper_object$run3
wrapper_object$run4


source('wrapper_object.R')
out <- wrapper_object$.test(mc=F, oconf=F, gen_metd=F )
out <- wrapper_object$.test(mc=T, oconf=F )
out <- wrapper_object$.test(mc=F, oconf=F )

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
wrapper_object$init_output_matrix_saltelli()
wrapper_object$write_output
wrapper_object$print_saltelli


source('wrapper_object.R')
out <- wrapper_object$.test_saltelli(mc=T, pr=6, oconf=F, n=10 )
wrapper_object$dataf$env
wrapper_object$dataf$fnames
wrapper_object$dataf$pars


# MCMC Mixture test
source('wrapper_object.R')
# ALJ: test DREAM algorithm with mixture model
out <- wrapper_object$.test_mcmc_mixture(mcmc_type = 'dream',
                                         mcmc_maxiter = 5000,
                                         mcmc_chains = 8,
                                         mc = F,
                                         pr = 4,
                                         mu_vector = c(-8, 0, 8),
                                         sd_vector = c(1, 1, 1),
                                         height_vector = c(0.1, 0.3, 0.6),
                                         mixture_scale = 1e12
                                         )
# ALJ: histogram of posterior parameter distributions
hist(out$pars_array,
     breaks = 200,
     col = 'darkmagenta',
     border = 'darkmagenta',
     xlab = 'Mixture Model Parameters',
     main = 'Posterior (Target) Parameter Distributions for Mixture Model')
# df1 <- data.frame(lklihood=as.vector(t(out[[2]])), chain=rep(1:dim(out[[2]])[1],each=dim(out[[2]])[2]) )
# xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) , df1, groups=chain, auto.key=T, type='l' )
# names(out)
# update(out[[3]],breaks=50)
# wrapper_object$dynamic
# wrapper_object$dynamic$pars_eval
# wrapper_object$dataf$pars
# wrapper_object$model$pars
# wrapper_object$dataf$pars_array[,,1]
# wrapper_object$dataf$pars_lklihood
# wrapper_object$dataf$pars_array
# wrapper_object$wpars


# MCMC linear regression test
source('wrapper_object.R')
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=150)
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=150, mcmc_type='dream' )
out <- wrapper_object$.test_mcmc_linreg(mcmc_maxiter=150,
                                        mcmc_test.a ='runif(n,-20,20)',
                                        mcmc_test.b ='runif(n,-20,20)')
#out <- wrapper_object$.test_mcmc_linreg(mc=T, pr=4, mcmc_chains=8, mcmc_maxiter=1000 )
update(out[[3]],breaks=50)
dim(out[[2]])
df1 <- data.frame(lklihood=as.vector(t(out[[2]])), chain=rep(1:dim(out[[2]])[1],each=dim(out[[2]])[2]) )
xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) , df1, groups=chain, auto.key=T, type='l' )
xyplot(lklihood ~ rep(1:dim(out[[2]])[2], dim(out[[2]])[1] ) | chain , df1, auto.key=T, type='l' )
names(out)
dim(out$pars_array)
histogram(as.numeric(out$pars_array[,1,]))
histogram(as.numeric(out$pars_array[,2,]))
class(out$pars_array[,2,])
as.numeric(out$pars_array[,2,])
wrapper_object$mcmc$boundary_min
wrapper_object$mcmc$boundary_max

wrapper_object$wpars
wrapper_object$mcmc
wrapper_object$dynamic
wrapper_object$dataf$pars
wrapper_object$model$pars
wrapper_object$dataf$pars_array
wrapper_object$dataf$pars_lklihood
wrapper_object$wpars
dim(wrapper_object$dataf$out_mcmc)
wrapper_object$dataf$out_mcmc[,,74]
wrapper_object$dataf$out_mcmc[,,75]


# Canopy
source('wrapper_object.R')
out <- wrapper_object$.test_can(mc=F,verbose=F)
head(out[[1]])
print(out[[2]])

wrapper_object$model$state
wrapper_object$model$leaf$fns$sys
wrapper_object$dataf$met



### END ###
