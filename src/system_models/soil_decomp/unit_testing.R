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
soil_decomp_object$state$cpools
soil_decomp_object$.test(steadystate=T)
soil_decomp_object$.test(verbose=T)
soil_decomp_object$.test(litter=2)
soil_decomp_object$.test(metdf=T, ntimes=2 )
soil_decomp_object$.test(metdf=T, ntimes=100 )

# model mimics
source('soil_decomp_object.R')
soil_decomp_object$.test(mod_mimic='mimics')
soil_decomp_object$.test(mod_mimic='mimics', steadystate=T )
soil_decomp_object$.test(mod_mimic='corpse')
soil_decomp_object$.test(mod_mimic='corpse', steadystate=T )
soil_decomp_object$.test(mod_mimic='century')
soil_decomp_object$.test(mod_mimic='century', steadystate=T )
soil_decomp_object$.test(mod_mimic='elmcentury')
soil_decomp_object$.test(mod_mimic='elmcentury', steadystate=T )
soil_decomp_object$.test(mod_mimic='elmctc')
soil_decomp_object$.test(mod_mimic='elmctc', steadystate=T )
soil_decomp_object$.test(mod_mimic='mend2013')
soil_decomp_object$.test(mod_mimic='mend2013', steadystate=T )
soil_decomp_object$.test(mod_mimic='millennialv2')
soil_decomp_object$.test(mod_mimic='millennialv2', steadystate=T )
soil_decomp_object$.test(mod_mimic='millennialv2_original')
soil_decomp_object$.test(mod_mimic='millennialv2_original', steadystate=T )

# model output
soil_decomp_object$state
soil_decomp_object$state_pars

# model configuration
soil_decomp_object$env
soil_decomp_object$pars
soil_decomp_object$fnames

# outfluxes associated with each pool
soil_decomp_object$fnames$decomp
soil_decomp_object$pars$decomp_outflux1
soil_decomp_object$pars$decomp_outflux2
soil_decomp_object$pars$decomp_outflux3
soil_decomp_object$pars$decomp_outflux4
soil_decomp_object$pars$decomp_outflux5
soil_decomp_object$pars$decomp_outflux6
soil_decomp_object$pars$decomp_outflux7

# model initialisation
soil_decomp_object$env
soil_decomp_object$pars$cstate0
soil_decomp_object$pars$input_coef
soil_decomp_object$fns$input
soil_decomp_object$fns$input()

# transfer matrix
soil_decomp_object$state_pars$transfer_matrix
soil_decomp_object$fns$transfermatrix
soil_decomp_object$fns$transfer.t1_to_4
soil_decomp_object$fns$cue.cue3(i=3)
soil_decomp_object$fns$cue2.cue23(i=3)

# scalars
soil_decomp_object$state_pars$matpot_sat
soil_decomp_object$fns$matpot_sat()
soil_decomp_object$fns$wcor.wcor2(i=2)
soil_decomp_object$fns$tcor.tcor2(i=2)
soil_decomp_object$fns$wcor.wcor2(i=2) * 
  soil_decomp_object$fns$tcor.tcor2(i=2) *
  soil_decomp_object$state_pars$vmax$vmax2
soil_decomp_object$fns$tcor.t2
soil_decomp_object$fns$tcor.t2(i=2)
soil_decomp_object$state_pars$vmax[[1]] * soil_decomp_object$fns$tcor.tcor1(i=1)
soil_decomp_object$state_pars$vmax[[6]] * soil_decomp_object$fns$tcor.tcor6(i=6)

# outflux calculation
soil_decomp_object$fns$outflux.of1(C=soil_decomp_object$state$cpools, t=1, i=1, of=1 )
soil_decomp_object$fns$outflux.of2(C=soil_decomp_object$state$cpools, t=1, i=1, of=2 )
soil_decomp_object$fns$outflux.of3(C=soil_decomp_object$state$cpools, t=1, i=2, of=3 )
soil_decomp_object$fns$outflux.of4(C=soil_decomp_object$state$cpools, t=1, i=3, of=4 )
soil_decomp_object$fns$outflux.of5(C=soil_decomp_object$state$cpools, t=1, i=3, of=5 )



### END ###