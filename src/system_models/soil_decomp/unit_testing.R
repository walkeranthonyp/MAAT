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
soil_decomp_object$.test()
soil_decomp_object$.test()

soil_decomp_object$.test(verbose=T)

soil_decomp_object$fnames
soil_decomp_object$fns
names(soil_decomp_object$fns)
soil_decomp_object$fns$as.list()

grep('transfer\\.', names(soil_decomp_object$fns), value=T )


soil_decomp_object$fns$decomp.d1
soil_decomp_object$fns$decomp.d1(soil_decomp_object$state$c_pools[,1])
soil_decomp_object$fns$decomp.d2
soil_decomp_object$fns[["decomp.d2"]]()
soil_decomp_object$fns$transfer.t1_to_2
soil_decomp_object$fns[["transfer.t1_to_2"]]
soil_decomp_object$fns$transfer.t2_to_1

soil_decomp_object$state
soil_decomp_object$state$c_pools[,1]
soil_decomp_object$state_pars
soil_decomp_object$pars
soil_decomp_object$env
soil_decomp_object$run

soil_decomp_object$state$solver_out


soil_decomp_object$env$times
soil_decomp_object$fns$transfermatrix
soil_decomp_object$fns$transfermatrix()
soil_decomp_object$fns$DotO
soil_decomp_object$fns$DotO()
soil_decomp_object$fns$DotO(soil_decomp_object$state$c_pools[,1])
soil_decomp_object$fns$DotC
soil_decomp_object$fns$DotC(soil_decomp_object$state$c_pools[,1])
soil_decomp_object$fns$inputrates
soil_decomp_object$fns$inputrates()




### END ###
