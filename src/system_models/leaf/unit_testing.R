###########################
#
# Unit testing of the MAAT models and wrapper
#
# AWalker
# May 2015
#
###########################



# Leaf MAAT
rm(list=ls())
library(lattice)

# leaf
source('leaf_object.R')
leaf_object$.test_leaf(verbose=F,verbose_loop=F)
leaf_object$fnames
leaf_object$state
leaf_object$state_pars
leaf_object$pars
leaf_object$env

source('leaf_object.R')
leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=600,rs='f_rs_medlyn2011',gd='f_ficks_ci')
leaf_object$state
leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=600,rs='f_rs_medlyn2011',gd='f_ficks_ci_bound0')
leaf_object$state
leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=1000,rs='f_rs_constant')
leaf_object$state
leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=1000,rs='f_rs_constant',gd='f_ficks_ci_bound0')
leaf_object$state

source('leaf_object.R')
leaf_object$.test_leaf(verbose=F,leaf.par=1000,leaf.ca_conc=1100)
leaf_object$.test_aci(leaf.ca_conc=seq(0.1,2000,50))
leaf_object$.test_aci(diag=T)
leaf_object$.test_aci(diag=T,leaf.par=c(450,900) )
leaf_object$.test_aci(diag=T,leaf.par=c(450,900,1200,1500) )
leaf_object$.test_aci_light()
leaf_object$.test_aci_light(diag=T)

source('leaf_object.R')
odf <- leaf_object$.test_aci_lim(leaf.par=2000,leaf.ca_conc=seq(100,1500,10))
odf <- leaf_object$.test_aci_lim(et='f_j_collatz1991',leaf.par=2000,leaf.ca_conc=seq(100,1500,10))

source('leaf_object.R')
leaf_object$.test_tscalar()
leaf_object$.test_tscalar(tcor_des='f_temp_scalar_modArrhenius_des')
leaf_object$.test_tscalar(tcor_asc='f_scalar_none')
leaf_object$.test_tscalar(tcor_asc='f_temp_scalar_Q10')
leaf_object$.test_tscalar(tcor_asc='f_temp_scalar_Q10',tcor_des='f_temp_scalar_collatz1991_des')
leaf_object$.test_tscalar(tcor_asc='f_temp_scalar_Q10',tcor_des='f_temp_scalar_cox2001_des')

source('leaf_object.R')
out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011')
out <- leaf_object$.test_aci_analytical(rs='f_rs_ball1987')
out <- leaf_object$.test_aci_analytical(rs='f_rs_leuning1995')
out <- leaf_object$.test_aci_analytical(rs='f_rs_constantCiCa')
out <- leaf_object$.test_aci_analytical(rs='f_rs_cox1998')

out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011',leaf.rb=0.5)
out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011',leaf.rb=0.1)
out <- leaf_object$.test_aci_analytical(rs='f_rs_leuning1995',leaf.rb=0.0001)

# leaf solver function
source('leaf_object.R')
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F)
leaf_object$fnames$Alim <- 'f_lim_collatz1991'
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F)
leaf_object$fnames$jmax <- 'f_jmax_lin'
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F)
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_ball1987')
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_ball1987',leaf.ca_conc=50)

leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000)
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=-30:200) 
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-30,200) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=-10:1 )

leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,1) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,10) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,100) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(1,-10) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(10,-10) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(100,-10) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=-10:100 )

leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,10) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,100) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_medlyn2011',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,10) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_medlyn2011',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,100) )
leaf_object$.test_solverFunc(verbose=F,verbose_loop=F,rs='f_rs_medlyn2011',leaf.ca_conc=1,leaf.par=1000 )

source('leaf_object.R')
leaf_object$.test_solverFunc(verbose=T,verbose_loop=T,rs='f_rs_leuning1995',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,100) )
leaf_object$.test_solverFunc(verbose=T,verbose_loop=T,rs='f_rs_medlyn2011',leaf.ca_conc=1,leaf.par=1000, sinput=c(-10,100) )




