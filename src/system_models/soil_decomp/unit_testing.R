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
soil_decomp_object$.test()
soil_decomp_object$.test(litter=6)
soil_decomp_object$.test(verbose=T)

soil_decomp_object$state
soil_decomp_object$state$cpools
soil_decomp_object$state_pars
soil_decomp_object$pars
soil_decomp_object$env
soil_decomp_object$run

soil_decomp_object$.test(metdf=T)
soil_decomp_object$.test(metdf=T, ntimes=10 )
soil_decomp_object$.test(metdf=T, litter=6, ntimes=10 )
soil_decomp_object$.test(metdf=T, litter=1:10 )

soil_decomp_object$fnames
soil_decomp_object$fns
names(soil_decomp_object$fns)
soil_decomp_object$fns$as.list()

soil_decomp_object$fns$decomp.d1
soil_decomp_object$fns$decomp.d1(soil_decomp_object$state$cpools[,1],i=1)
soil_decomp_object$fns$decomp.d2
soil_decomp_object$fns$decomp.d3

soil_decomp_object$fns$transfer.t1_to_2
soil_decomp_object$fns[["transfer.t1_to_2"]]

soil_decomp_object$fns$inputrates
soil_decomp_object$fns$inputrates()
soil_decomp_object$fns$transfermatrix
soil_decomp_object$fns$transfermatrix()
soil_decomp_object$fns$DotO
soil_decomp_object$fns$DotO(soil_decomp_object$state$cpools[,1])


source('soil_decomp_object.R')
system.time(olist <- soil_decomp_object$.test_3pool(ntimes=365))
system.time(olist <- soil_decomp_object$.test_3pool())
olist



### END ###