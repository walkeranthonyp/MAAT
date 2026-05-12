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

#soil_decomp_object$.test.mimics.ss()
soil_decomp_object$.test()
soil_decomp_object$.test(litter=172.8978/365/1000)
soil_decomp_object$.test(litter=.5/.15/365)
soil_decomp_object$state$cpools
soil_decomp_object$.test(litter=10)

soil_decomp_object$fns$outflux.of1
soil_decomp_object$fns$cue.cue3(i=3)
soil_decomp_object$fns$cue2.cue23(i=3)
names(soil_decomp_object$fns)
soil_decomp_object$fns$calc_state_pars
soil_decomp_object$fns[['outflux.of11']]
soil_decomp_object$fns$tcor
soil_decomp_object$fns$input()

unlist(soil_decomp_object$state)
soil_decomp_object$state$cpools
soil_decomp_object$state$outflux
soil_decomp_object$fnames
soil_decomp_object$fnames$cue
soil_decomp_object$fnames$outflux
soil_decomp_object$fnames$transfer
soil_decomp_object$pars
soil_decomp_object$fnames$k
soil_decomp_object$pars$k
soil_decomp_object$state_pars$k
soil_decomp_object$fnames$km
soil_decomp_object$pars$km
soil_decomp_object$state_pars$km
soil_decomp_object$pars$cue
soil_decomp_object$pars$cue2
soil_decomp_object$state_pars$cue
soil_decomp_object$state_pars$cue2
soil_decomp_object$pars$cat_pool
soil_decomp_object$pars$cstate0
soil_decomp_object$pars$input_coef
soil_decomp_object$pars$millennialV2

soil_decomp_object$env

names(soil_decomp_object$fns)
soil_decomp_object$fns$decomp.d1
soil_decomp_object$fns$decomp.d1(soil_decomp_object$state$cpools, i=1 )
soil_decomp_object$fns$decomp.d2(soil_decomp_object$state$cpools, i=2 )
soil_decomp_object$fns$transfer_fluxsum_prop(from=2)
soil_decomp_object$fns[["transfer_fluxsum_prop"]](from=2)
soil_decomp_object$fns[["transfer_fluxsum_prop"]](.=soil_decomp_object$fns,from=2)
  
soil_decomp_object$pars$millennialV2

soil_decomp_object$pars$decomp_outflux1
unlist(soil_decomp_object$pars$decomp_outflux1)[[2]]
soil_decomp_object$pars$decomp_outflux2
soil_decomp_object$pars$decomp_outflux3
soil_decomp_object$pars$decomp_outflux4
soil_decomp_object$pars$decomp_outflux5
soil_decomp_object$pars$decomp_outflux6
soil_decomp_object$pars$decomp_outflux7

soil_decomp_object$state$cpools[] <- 1
soil_decomp_object$fns$outflux.of1(C=soil_decomp_object$state$cpools, t=1, i=1, of=1 )
soil_decomp_object$fns$outflux.of2(C=soil_decomp_object$state$cpools, t=1, i=1, of=2 )
soil_decomp_object$fns$outflux.of3(C=soil_decomp_object$state$cpools, t=1, i=2, of=3 )
soil_decomp_object$fns$outflux.of4(C=soil_decomp_object$state$cpools, t=1, i=3, of=4 )
soil_decomp_object$fns$outflux.of5(C=soil_decomp_object$state$cpools, t=1, i=3, of=5 )
soil_decomp_object$fns$outflux.of6(C=soil_decomp_object$state$cpools, t=1, i=4, of=6 )
soil_decomp_object$fns$outflux.of7(C=soil_decomp_object$state$cpools, t=1, i=4, of=7 )
soil_decomp_object$fns$outflux.of8(C=soil_decomp_object$state$cpools, t=1, i=4, of=8 )
soil_decomp_object$fns$outflux.of9(C=soil_decomp_object$state$cpools, t=1, i=5, of=9 )
soil_decomp_object$fns$outflux.of4
soil_decomp_object$fnames$outflux

soil_decomp_object$fns$tcor.t2
soil_decomp_object$fns$tcor.t2(i=2)

soil_decomp_object$fnames$transfer
soil_decomp_object$state$outflux[] <- 1
soil_decomp_object$fns$transfer.t
soil_decomp_object$fns$transfer.t1_to_4
soil_decomp_object$fns$transfer.t1_to_4(from=1,to=4)
soil_decomp_object$fns$transfer.t1_to_5(from=1,to=5)
soil_decomp_object$fns$transfer.t2_to_3(from=2,to=3)
soil_decomp_object$fns$transfer.t2_to_4(from=2,to=4)
soil_decomp_object$fns$transfer.t3_to_4(from=3,to=4)
soil_decomp_object$fns$transfer.t3_to_5(from=3,to=5)
soil_decomp_object$fns$transfer.t4_to_2(from=4,to=2)
soil_decomp_object$fns$transfer.t4_to_3(from=4,to=3)
soil_decomp_object$fns$transfer.t5_to_1(from=5,to=1)
soil_decomp_object$fns$transfer.t5_to_3(from=5,to=3)

soil_decomp_object$init1

soil_decomp_object$fns$input()

soil_decomp_object$pars$decomp_outflux4

soil_decomp_object$state$cpools
soil_decomp_object$state$outflux
soil_decomp_object$state$outflux[[2]]
soil_decomp_object$state$outflux[[1:2]]
soil_decomp_object$state$outflux[] <- 1

soil_decomp_object$pars$k
soil_decomp_object$pars$k[[1]]
soil_decomp_object$pars$km
soil_decomp_object$pars$vmax
soil_decomp_object$pars$cat_pool
soil_decomp_object$pars$poolmax

soil_decomp_object$pars[[paste0('decomp_outflux',1)]]
unlist(soil_decomp_object$pars[[paste0('decomp_outflux',1)]])
soil_decomp_object$state$outflux[
  unlist(soil_decomp_object$pars[[paste0('decomp_outflux',1)]])
]
soil_decomp_object$fns$transfer.t1_to_3

a1 <- unlist(soil_decomp_object$pars[[paste0('decomp_outflux',1)]])
class(a1)

soil_decomp_object$.test(litter=172.8978/365, ntimes = 1, metdf = T)
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


source('soil_decomp_object.R')
soil_decomp_object$.test_ctc()
soil_decomp_object$pars$n_pools
soil_decomp_object$pars$cstate0
soil_decomp_object$state_pars$solver_out
soil_decomp_object$state_pars$solver_steadystate_out



### END ###