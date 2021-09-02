###########################
#
# soil_decomp for MAAT object unit testing
# 
# AWalker March 2017
#
###########################



### Load model scripts 
###############################

# soil_decomp
source('soil_decomp_object.R')
soil_decomp_object$.test()
soil_decomp_object$.test(litter=1.4)
soil_decomp_object$.test(litter=6)
soil_decomp_object$.test(verbose=T)

soil_decomp_object$state
soil_decomp_object$state$cpools
soil_decomp_object$state_pars
soil_decomp_object$pars
soil_decomp_object$env
soil_decomp_object$run


soil_decomp_object$.test(metdf=T)
soil_decomp_object$state_pars$solver_steadystate_out
soil_decomp_object$state_pars$solver_steadystate_out$y
soil_decomp_object$.test(metdf=T, ntimes=10 )

soil_decomp_object$.test(metdf=T, litter=6, ntimes=10 )
soil_decomp_object$.test(metdf=T, ntimes = 1 )

soil_decomp_object$fnames
soil_decomp_object$fns$transfer.t1_to_3
names(soil_decomp_object$fns)
soil_decomp_object$fns$as.list()
names(soil_decomp_object$fns$as.list())

soil_decomp_object$fns$decomp.d1
soil_decomp_object$fns$decomp.d1(soil_decomp_object$state$cpools,i=1)
soil_decomp_object$fns$decomp.d2
soil_decomp_object$fns$decomp.d3

soil_decomp_object$fns$transfer.t1_to_2
soil_decomp_object$fns[["transfer.t1_to_2"]]

soil_decomp_object$fns$input
soil_decomp_object$fns$input()
soil_decomp_object$fns$transfermatrix
soil_decomp_object$fns$transfermatrix(C=soil_decomp_object$state$cpools)
soil_decomp_object$fns$DotO
soil_decomp_object$fns$DotO(soil_decomp_object$state$cpools[,1])

source('soil_decomp_object.R')
soil_decomp_object$.test_changepool()

source('soil_decomp_object.R')
system.time(olist <- soil_decomp_object$.test_3pool(ntimes=2))
system.time(olist <- soil_decomp_object$.test_3pool(ntimes=365))
system.time(olist <- soil_decomp_object$.test_3pool(ntimes=3650))
system.time(olist <- soil_decomp_object$.test_3pool())
olist

source('soil_decomp_object.R')
#default MM (vmax1 = 88, vmax3 = 171, km1 = 144, km3 = 936) from MIMICS
soil_decomp_object$.test_var_kinetics_yearly(kinetics = 'mm')
#equivalent RMM
soil_decomp_object$.test_var_kinetics_yearly(kinetics = 'rmm', vmax1 = 4, vmax3 = 6, km1 = 6, km3 = 33)
#MM with low vmax and km
soil_decomp_object$.test_var_kinetics_yearly(kinetics = 'mm', vmax1 = 44, vmax3 = 86, km1 = 72, km3 = 468)
#RMM with low vmax and km
soil_decomp_object$.test_var_kinetics_yearly(kinetics = 'rmm', vmax1 = 2, vmax3 = 3, km1 = 3, km3 = 17)
#linear
soil_decomp_object$.test_var_kinetics_yearly(kinetics = 'lin', k1 = .3, k3 = .09)



### END ###