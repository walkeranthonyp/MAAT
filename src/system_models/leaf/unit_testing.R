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

source('leaf_object.R')
leaf_object$.test(verbose=F, verbose_loop=F )
leaf_object$fnames
leaf_object$state
leaf_object$state_pars
leaf_object$pars
leaf_object$env
leaf_object$.test(verbose=F, leaf.par=1000, leaf.ca_conc=1100 )
leaf_object$.test(verbose=F, leaf.par=1000, leaf.ca_conc=1 )

source('leaf_object.R')
leaf_object$.test_aci(leaf.ca_conc=seq(0.1,2000,50), verbose=F )
leaf_object$.test_aci(leaf.ca_conc=seq(0.1,2000,50), rs='f_rs_yin2009' )
leaf_object$.test_aci(diag=T)
leaf_object$.test_aci(diag=T, leaf.par=c(450,900) )
leaf_object$.test_aci(diag=T, leaf.par=c(450,900,1200,1500) )
leaf_object$.test_aci_light()
leaf_object$.test_aci_light(diag=T)

source('leaf_object.R')
out <- leaf_object$.test_aci_lim(leaf.par=2000, leaf.ca_conc=seq(100,1500,10) )
out <- leaf_object$.test_aci_lim(et='f_j_collatz1991', leaf.par=2000, leaf.ca_conc=seq(100,1500,10) )

source('leaf_object.R')
leaf_object$.test_tscalar()
leaf_object$.test_tscalar(Ha=40000)
leaf_object$.test_tscalar(tcor_des='f_tcor_des_modArrhenius')
leaf_object$.test_tscalar(tcor_des='f_tcor_des_modArrhenius', Ha=40000 )
leaf_object$.test_tscalar(tcor_asc='f_scalar_none')
leaf_object$.test_tscalar(tcor_asc='f_tcor_asc_Q10')
leaf_object$.test_tscalar(tcor_asc='f_tcor_asc_Q10', tcor_des='f_tcor_des_collatz1991' )
leaf_object$.test_tscalar(tcor_asc='f_tcor_asc_Q10', tcor_des='f_tcor_des_cox2001' )

source('leaf_object.R')
out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011')
out <- leaf_object$.test_aci_analytical(rs='f_rs_ball1987', leaf.rb=0.115293015 )
out <- leaf_object$.test_aci_analytical(rs='f_rs_leuning1995')
out <- leaf_object$.test_aci_analytical(rs='f_rs_constantCiCa')
out <- leaf_object$.test_aci_analytical(rs='f_rs_cox1998')
out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011', leaf.rb=0.5 )
out <- leaf_object$.test_aci_analytical(rs='f_rs_medlyn2011', leaf.rb=0.1 )
out <- leaf_object$.test_aci_analytical(rs='f_rs_leuning1995', leaf.rb=0.0001 )

# leaf solver function
source('leaf_object.R')
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_ball1987' )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_ball1987', leaf.ca_conc=50 )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_ball1987', leaf.ca_conc=50, sinput=seq(-3,10,0.1) )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_ball1987', leaf.ca_conc=300, leaf.par=100, sinput=seq(-2,5,0.01) ) 
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_leuning1995', leaf.ca_conc=1, leaf.par=1000 )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_leuning1995', leaf.ca_conc=1, leaf.par=1000, sinput=-30:200) 
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_medlyn2011', leaf.ca_conc=1, leaf.par=1000, sinput=c(-10,10) )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_medlyn2011', leaf.ca_conc=1, leaf.par=1000, sinput=c(-10,100) )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_medlyn2011', leaf.ca_conc=100, leaf.par=1000, seq(-1,5,0.01) )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_yin2009', leaf.ca_conc=50, leaf.par=1000 )
leaf_object$.test_solverFunc(verbose=F, verbose_loop=F, rs='f_rs_yin2009', leaf.ca_conc=500, leaf.par=1000 )



### END ###