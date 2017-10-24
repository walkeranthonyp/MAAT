################################
#
# Leaf level physiology functions
# 
# AWalker March 2014
#
################################


f_none <- function(.) {
  NA
}


# VPD
################################
f_rh_from_vpd <- function(.) {
  rh <- (1 - .$env$vpd / f_sat_vp_allen(.))
  max(rh,0)
}

f_sat_vp_allen <- function(.) {
  # calculate saturation vapour pressure - Allen 1998
  # in kPa
  # t -- oC -- air temp
  
  0.6108 * exp( 17.27*.$state$leaf_temp /(.$state$leaf_temp+237.3) )
}

sat_vp_buck <- function(t){
  # calculate saturation vapour pressure - Buck 1981 J. App. Met.
  # in kPa
  # t -- oC -- air temp
  
  0.61121 * exp( 17.502*t / (240.97 + t) )
}



### ANALYTICAL SOLUTIONS
################################

quad_sol <- function(a,b,c,out='lower') {
  # robust numerical solution to the quadratic
  # taken from Numerical Recipes

  q     <- -0.5 * ( b + sign(b)*(b^2 - 4*a*c)^0.5 )
  roots <- c( q/a , c/q )
  
  if(out=='lower')      min(roots,na.rm=T) 
  else if(out=='upper') max(roots,na.rm=T)
  else roots 
}

f_A_r_leaf_analytical <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes
  # solves A analytically by assuming g0 = 0 in the stomatal resistance function, and rb and ri are also assumed zero 
  
  .$state_pars$rb <- 0
  .$state_pars$ri <- 0
  
  # calculate cb, ci & cc
  .$state$cb      <- .$state$ca
  fe              <- get(paste(.$fnames$rs,'fe',sep='_') )(.) 
  .$state$ci      <- .$state$ca * (1 - (1.6 / fe) )
  .$state$cc      <- .$state$ci 
  
  # calculate Ag / cc for each limiting process
  .$state$Acg     <- get(.$fnames$Acg)(.)
  .$state$Ajg     <- get(.$fnames$Ajg)(.)
  .$state$Apg     <- get(.$fnames$Apg)(.)
  
  # calculate limiting cycle
  Amin            <- get(.$fnames$Alim)(.) 
  
  # calculate Ag (gross asimilation) for each limiting process
  .$state$Acg     <- .$state$Acg * .$state$cc 
  .$state$Ajg     <- .$state$Ajg * .$state$cc 
  .$state$Apg     <- .$state$Apg * .$state$cc
  
  # calculate net A
  Anet <- Amin*.$state$ci - Amin*.$state_pars$gstar - .$state$respiration
  
  # calculate rs
  .$state_pars$rs <- .$state$ca / (fe * Anet) 
  
  # return net A
  Anet
}

f_A_r_leaf_noR <- function(.,...) {
  # same as above function but with no resistance to CO2 diffusion
  # from the atmosphere to the site of carboxylation
  # These functions can be used to calculate the stomatal limitation to photosynthesis
  
  # combines all rate limiting processes
  
  # assume cc = ca
  cc <- .$state$ca
  
  # calculate Ag/cc
  # - Ag/cc and not Ag is calculated for numerical stability at low cc
  .$state$Acg <- get(.$fnames$Acg)(.,cc=cc)
  .$state$Ajg <- get(.$fnames$Ajg)(.,cc=cc)
  .$state$Apg <- get(.$fnames$Apg)(.,cc=cc)
  
  # calculate limiting cycle
  Amin <- get(.$fnames$Alim)(.) 
  
  # calculate actual w
  .$state$Acg <- .$state$Acg * cc
  .$state$Ajg <- .$state$Ajg * cc
  .$state$Apg <- .$state$Apg * cc
  
  # calculate net A
  Amin*cc - Amin*.$state_pars$gstar - .$state$respiration
}



### SOLVERS & RESIDUAL FUNCTIONS
################################

f_R_Brent_solver <- function(.,...) {
  # ... could be used to pass through different functions 'func(.$env$atm_press * 1.6e-6)' to the solver function, not currently necessary

  if(.$cpars$verbose_loop) print(.$env)  
  .$solver_out <- uniroot(get(.$fnames$solver_func),interval=c(-10,100),.=.,extendInt='no',...)
  .$solver_out$root
}

f_A_r_leaf <- function(A,.,...) {
  # combines A, rs, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  #  -- for use with uniroot solver
  #  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation

  # calculate cc from ca, rb, rs, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer and stomatal resistance terms are in h2o units
  # assumes mesophyll resistance is in co2 units
  # print(.$fnames$gas_diff)
  # print(.$fnames$rs)
  # print(.$state_pars$rb)
  # print(1.6*get(.$fnames$rs)(.,A=A,c=get(.$fnames$gas_diff)(.,A)))
  # print(.$state_pars$ri)
  .$state$cc <- get(.$fnames$gas_diff)( . , A , r=( 1.4*.$state_pars$rb + 1.6*get(.$fnames$rs)(.,A=A,c=get(.$fnames$gas_diff)(.,A)) + .$state_pars$ri ) )
  
  # calculate Ag/cc
  # - Ag/cc and not Ag is calculated for numerical stability at low cc
  .$state$Acg <- get(.$fnames$Acg)(.)
  .$state$Ajg <- get(.$fnames$Ajg)(.)
  .$state$Apg <- get(.$fnames$Apg)(.)

  # calculate limiting cycle
  Amin <- get(.$fnames$Alim)(.) 
  
  # calculate actual w
  .$state$Acg <- .$state$Acg * .$state$cc
  .$state$Ajg <- .$state$Ajg * .$state$cc
  .$state$Apg <- .$state$Apg * .$state$cc
  
  # calculate residual of net A
  Amin*.$state$cc - Amin*.$state_pars$gstar - .$state$respiration - A
} 

f_A_r_leaf_noRs <- function(A,.,...) {
  # same as above function but with no stomatal resistance 
  
  # combines A, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  #  -- for use with uniroot solver
  #  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation
  
  # calculate cc from ca, rb, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer resistance is in h2o units
  # assumes mesophyll resistance is in co2 units
  cc <- get(.$fnames$gas_diff)( . , A , r=( 1.4*.$state_pars$rb + .$state_pars$ri ) )
  
  # calculate Ag/cc
  # - Ag/cc and not Ag is calculated for numerical stability at low cc
  .$state$Acg <- get(.$fnames$Acg)(.,cc=cc)
  .$state$Ajg <- get(.$fnames$Ajg)(.,cc=cc)
  .$state$Apg <- get(.$fnames$Apg)(.,cc=cc)
  
  # calculate limiting cycle
  Amin <- get(.$fnames$Alim)(.) 
  
  # calculate actual w
  .$state$Acg <- .$state$Acg * cc
  .$state$Ajg <- .$state$Ajg * cc
  .$state$Apg <- .$state$Apg * cc
  
  # calculate residual of net A
  Amin*cc - Amin*.$state_pars$gstar - .$state$respiration - A
}



### PHOTOSYNTHESIS FUNCTIONS
################################

# Transition point function, calculates the cc at which Acg = Ajg, i.e. the transition cc
transition_cc <- function(.) {
  
  vcm_et_ratio <- .$state_pars$vcmaxlt/.$state$J 
  
  ( 8*.$state_pars$gstar*vcm_et_ratio - .$state_pars$Km ) /
    (1 - 4*vcm_et_ratio)
}


# Carboxylation limitation
f_Acg_farquhar1980 <- function(.,cc=.$state$cc){   
  # Farquhar 1980 eq to calculate RuBisCO limited photosynthetic rate
  # umol m-2 s-1
  # Wc<-(Vcmax*(Ci-Gamma))/(Ci+Kc*(1+(O/Ko)))
  
  # dvars
  # Ci (or Cc but vcmax needs adjusting)   
  
  # params
  # Vcmax maximum carboxylation rate of RuBisCO  (umol m-2 s-1
  # Ko Michaelis constant for O2 mbar            (kPa)
  # Kc Michaelis constant for CO2                ( Pa) 
  
  # Ko and O must be in same units
  # Kc, Ci and Gamma must be in same units
  
  #gross of rd and G*
  .$state_pars$vcmaxlt /
    (cc + .$state_pars$Km)
}


# electron transport limitation
f_Ajg_generic <- function(.,cc=.$state$cc){
  # generic eq to calculate light limited photosynthetic rate from user selected electron transport and etoc function
  # umol m-2 s-1
  
  # calculate gross electron transport limited carboxylation 
  # i.e. rd and G* not included
  .$state$J / (4*(cc+2*.$state_pars$gstar))        # currently no other formulation exists (other than slighly different parameters in the denominator)
}

f_j_farquhar1980 <- function(.){
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per incident photon
  # Farquhar 1980 and others
  # this is the replacement for eq A2 added in the proofing stage (very last line of the manuscript after all refs etc)
  
  .$state_pars$jmaxlt * .$env$par * (1 - .$pars$f)/(.$env$par + 2.2*.$state_pars$jmaxlt)
}

f_j_harley1992 <- function(.){
  # Harley etal 1992 eq to calculate electron transport rate
  # umol e m-2 s-1
  harley_alpha <- .$pars$a * .$state_pars$alpha
  harley_alpha*.$env$par / ( 1.0 + (harley_alpha**2)*(.$env$par**2) / (.$state_pars$jmaxlt**2) )**0.5  
}

f_j_farquharwong1984 <- function(.){
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per photon
  # von Caemmerer 2000 book and taken from Farquhar & Wong 1984
  
  # in von Caemmerer 2000 - I2 = I * a*(1-f)/2
  # where a is leaf light absorptance, 1-f corrects for spectral quality and light not absorbed by photosystems, and /2 accounts for the 2 photons needed to fully tranport 1 electron 
  # this is presumably the same as I * alpha
  # if so von C 2000's alpha is 0.36125
    
  I2 <- .$env$par * .$pars$a * .$state_pars$alpha    
  a  <- .$pars$theta_j
  b  <- -1 * (I2 + .$state_pars$jmaxlt)
  c  <- I2 * .$state_pars$jmaxlt
  
  quad_sol(a,b,c)
}

f_j_collatz1991 <- function(.){
  # Collatz etal 1991 eq to calculate electron transport rate
  # umol e m-2 s-1
  
  .$pars$a * .$state_pars$alpha * .$env$par  
}


# TPU limitation
f_Apg_vonc2000 <- function(.,cc=.$state$cc){
  # triose phosphate limitation from vonCaemmerer 2000 as corrected in Gu 2010
  # this is derived from Harley & Sharkey 1991
  # the difference is that Harley & Sharkey iteratively solved to find tpu, while here tpu is set independently
  
  # As described in Gu 2010, the collatz 1991 function for TPU limitation is a special case of this one where:
  # alpha = 0, tpu = vcmax/6
  # and tpu temperature scaling is identical to vcmax
  
  if( cc <= (1+3*.$pars$Apg_alpha)*.$state_pars$gstar) NA
  else 3*.$state_pars$tpu / ( cc-(1+3*.$pars$Apg_alpha)*.$state_pars$gstar )
}

f_Apg_foley1996 <- function(.,cc=.$state$cc){
  # triose phosphate limitation from Foley 1996, citing Harlkey & Sharkey 1991
  # its not clear to me how they derive this from Harley and Sharkey
  # Eq 5 from Foley 1996: Js = 3 * tpu * (1 - gstar/cc) + Jp * gstar / cc
  # where Jp is the limiting rate of wc or wj
  
  if( cc <= .$state_pars$gstar ) NA
  else {
    # calculate limiting cycle of wc and wj
    .$state$Apg <- NA
    wmin <- get(.$fnames$Alim)(.) 

    (3*.$state_pars$tpult + wmin) / cc
  }
}


# limiting factor
f_lim_farquhar1980 <- function(.){
  # simple minimum of all three possible limitating states
  
  min(c(.$state$Acg,.$state$Ajg,.$state$Apg),na.rm=T)
}

f_lim_collatz1991 <- function(.){
  # smoothed solution of all three possible limiting states

  a  <- .$pars$theta_col_cj
  b  <- -1 * (.$state$Acg + .$state$Ajg)
  c  <- .$state$Acg * .$state$Ajg
  
  sol1 <- quad_sol(a,b,c)

  if(!is.na(.$state$Apg)) {
    a  <- .$pars$theta_col_cjp
    b  <- -1 * (.$state$Apg + sol1)
    c  <- .$state$Apg * sol1
    
    quad_sol(a,b,c)
  } else sol1
}

 

#########################
### Respiration functions

f_rd_constant <- function(.) {
  .$pars$atref.rd
}

f_rd_lin_vcmax <- function(.) {
  .$pars$a_rdv_25 + .$state_pars$vcmax * .$pars$b_rdv_25    
}

f_rd_lin_vcmax_t <- function(.) {
  # rd25 is a proportion of vcmax 25 but that proportion changes with temperature
  .$pars$b_rdv_25 <- .$pars$a_rdv_25_t + .$state$leaf_temp * .$pars$b_rdv_25_t
  .$pars$a_rdv_25 + .$state_pars$vcmax * .$pars$b_rdv_25    
}

f_rd_lin_N <- function(.) {
  .$pars$a_rdn_25 + .$state$leafN_area * .$pars$b_rdn_25    
}

f_rd_tcor_independent <- function(.) {
  # TPU temperature scaling is independent

  get(.$fnames$rd_tcor_asc)(.,parlist=list(Tr=.$pars$reftemp.rd,Ha=.$pars$Ha.rd,q10_func=.$fnames$q10_func.rd,
                                           q10=.$pars$q10.rd,a_q10_t=.$pars$a_q10_t.rd,b_q10_t=.$pars$b_q10_t.rd)) *
  get(.$fnames$rd_tcor_des)(.,parlist=list(Tr=.$pars$reftemp.rd,Hd=.$pars$Hd.rd,Topt=.$pars$Topt.rd,
                                           tupp=.$pars$tupp_cox.rd,tlow=.$pars$tlow_cox.rd,exp=.$pars$exp_cox.rd,
                                           a_deltaS_t=.$pars$a_deltaS_t.rd,b_deltaS_t=.$pars$b_deltaS_t.rd,deltaS=.$pars$deltaS.rd))
}

f_rd_tcor_dependent <- function(.) {
  # TPU temperature scaling is identical to that of Vcmax
  
  .$state_pars$vcmaxlt / .$state_pars$vcmax
}

# light supression
f_rl_rd_fixed <- function(.) {
  .$pars$rl_rd_ratio
}

f_rl_rd_lloyd1995 <- function(.) {
  
  if(.$env$par <= 10) 1
  else                ( .$pars$rl_rd_lloyd_a - .$pars$rl_rd_lloyd_b * ln(.$env$par) )
}



### PHOTOSYNTHESIS PARAMETER FUNCTIONS
################################

# vcmax
f_vcmax_constant <- function(.) {
  .$pars$atref.vcmax
}

f_vcmax_lin <- function(.) {
  if(.$cpars$verbose) print('Vcmax_CLM_function')
  .$pars$avn_25 + .$state$leafN_area * .$pars$bvn_25    
}

f_vcmax_clm <- function(.) {
 
  .$state$leafN_area * .$pars$flnr * .$pars$fnr * .$pars$Rsa  
}

# jmax
f_jmax_constant <- function(.) {
  .$pars$atref.jmax
}

f_jmax_power <- function(.) {
  exp(.$pars$e_ajv_25) * .$state_pars$vcmax^.$pars$e_bjv_25    
}

f_jmax_lin <- function(.) {
  .$pars$ajv_25 + .$state_pars$vcmax * .$pars$bjv_25    
}

f_jmax_lin_t <- function(.) {
  .$pars$ajv_25 <- 0.0
  .$pars$bjv_25 <- .$pars$a_jvt_25 + .$state$leaf_temp * .$pars$b_jvt_25
  .$pars$ajv_25 + .$state_pars$vcmax * .$pars$bjv_25    
}


# TPU
f_tpu_constant <- function(.) {
  .$pars$atref.tpu
}

f_tpu_lin <- function(.) {
  .$pars$atv_25 + .$state_pars$vcmax * .$pars$btv_25    
}

f_tpu_tcor_independent <- function(.) {
  # TPU temperature scaling is independent
  
  get(.$fnames$tpu_tcor_asc)(.,parlist=list(Tr=.$pars$reftemp.tpu,Ha=.$pars$Ha.tpu,q10=.$pars$q10.tpu)) *
  get(.$fnames$tpu_tcor_des)(.,parlist=list(Tr=.$pars$reftemp.tpu,Hd=.$pars$Hd.tpu,Topt=.$pars$Topt.tpu,
                                            a_deltaS_t=.$pars$a_deltaS_t.tpu,b_deltaS_t=.$pars$b_deltaS_t.tpu,deltaS=.$pars$deltaS.tpu))
}

f_tpu_tcor_dependent <- function(.) {
  # TPU temperature scaling is identical to that of Vcmax
  .$state_pars$vcmaxlt / .$state_pars$vcmax
}



### STOMATAL & RELATED CONDUCTANCE / RESISTANCE FUNCTIONS
################################

# CO2 diffusion
f_ficks_ci <- function(.,A=.$state$A,r=1.4*.$state_pars$rb,c=.$state$ca) {
  # can be used to calculate cc, ci or cs (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1
  
  c-A*r*.$env$atm_press*1e-6
  #   c-A*r*.$env$atm_press*1e-6
}

f_ficks_ci_bound0 <- function(.,A=.$state$A,r=1.4*.$state_pars$rb,c=.$state$ca) {
  # can be used to calculate cc, ci or cs (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1
  
  # there is a catch 22 here, make this max of the function result or zero and it screws up the solver boundaries
  # remove the max term and this screws up at very low ca
  #   max( c-A*r*.$env$atm_press*1e-6 , 1e-6)
  c2 <- c-A*r*.$env$atm_press*1e-6
  if(c2>0) c2 else 0
}


# general
f_conv_ms_molm2s1 <- function(.){
  .$env$atm_press / ((.$state$leaf_temp + 273.15)*.$pars$R)
}

f_r_zero <- function(.,...){
  0
}

f_r_zero_fe <- function(.,...){
  0
}


# stomata
# stomatal resistances are all assumed to be in h2o units 
f_rs_constant <- function(.,...) {
  # output in m2s mol-1 h2o

  # this currently doesn't work with the solver.
  # When rs is not a function of A the solver interval doesn't span zero,
  # if this is a function to be used must reroute solver to an analytical method

  .$pars$rs
}


f_rs_medlyn2011 <- function(.,A=.$state$A,c=.$state$cb){
  # Medlyn et al 2011 eq for stomatal resistance
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  # if( A < 0 ) 1/.$pars$g0
  # else 
  1 / (.$pars$g0 + f_rs_medlyn2011_fe(.) * A * .$env$atm_press*1e-6 / c )
}

f_rs_medlyn2011_fe <- function(.) {
  # f(e) component of rs from Medlyn 2011   
  ( 1 + .$pars$g1_medlyn / .$env$vpd^0.5 )
}


f_rs_leuning1995 <- function(.,A=.$state$A,c=.$state$cb){
  # Leuning et al 1995 eq for stomatal resistance
  # expects c in Pa
  # output in m2s mol-1  h2o

  # if( A < 0 ) 1/.$pars$g0
  # else 
  1 / ( .$pars$g0 + f_rs_leuning1995_fe(.,c=c) * A * .$env$atm_press*1e-6 / c ) 
}

f_rs_leuning1995_fe <- function(.,c=.$state$cb) {
  # f(e) component of rs from Leuning 1995   
  .$pars$g1_leuning / ( (1 - .$state_pars$gamma/c) * (1 + .$env$vpd/.$pars$d0) )  
}


f_rs_ball1987 <- function(.,A=.$state$A,c=.$state$cb){
  # Ball et al 1987 eq for stomatal resistance
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  # if( A < 0 ) 1/.$pars$g0
  # else
  1 / ( .$pars$g0 + f_rs_ball1987_fe(.) * A * .$env$atm_press*1e-6 / c )
}

f_rs_ball1987_fe <- function(.) {
  # f(e) component of rs from Ball 1987   
  .$pars$g1_ball * .$env$rh
}


f_rs_constantCiCa <- function(.,A=.$state$A,c=.$state$cb) {
  # implied LPJ etc assumption for stomatal resistance that keeps Ci:Ca constant
  # expects c in Pa
  # output in m2s mol-1 h2o

  # set Ci:Ca ratio
  .$state_pars$cica_chi <- get(.$fnames$cica_ratio)(.)
  
  # if( A < 0 ) 1/1e-9
  # else
  1 / (f_rs_constantCiCa_fe(.) * A * .$env$atm_press*1e-6 / c )
}

f_rs_constantCiCa_fe <- function(.) {
  # f(e) component of rs for constant Ci:Ca   
  .$state_pars$cica_chi <- get(.$fnames$cica_ratio)(.)
  1.6 / (1 - .$state_pars$cica_chi)
}

f_cica_constant <- function(.) {
  .$pars$cica_chi
}


f_rs_cox1998 <- function(.,A=.$state$A,c=.$state$cb) {
  # implied JULES etc assumption for stomatal resistance that keeps a variant of Ci:Ca constant
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  f0    <- 1 - 1.6/.$pars$g1_leuning
  dstar <- (.$pars$g1_leuning/1.6 - 1) * .$pars$d0
  # CmCP  <- (c - .$state_pars$gamma)
  CmCP  <- (1 - .$state_pars$gamma/c)
  
  # if( A < 0 ) 1/1e-9
  # else
  1 / ( f_rs_cox1998_fe(.,c=c) * A * .$env$atm_press*1e-6 / c )
}

f_rs_cox1998_fe <- function(.,c=.$state$cb) {
  # f(e) component of rs for Cox 1998   
  
  f0    <- 1 - 1.6/.$pars$g1_leuning
  dstar <- (.$pars$g1_leuning/1.6 - 1) * .$pars$d0
  CmCP  <- (1 - .$state_pars$gamma/c)

  1.6 / ( CmCP - f0*CmCP * (1 - .$env$vpd/dstar))
}


# internal/mesophyll
# internal/mesophyll resistances are all assumed by the solver to be in co2 units 
f_ri_constant <- function(., ... ) {
  # output in m2s mol-1 co2
  
  .$pars$ri
}

# leaf boundary layer
# leaf boundary layer resistances are all assumed by the solver to be in h2o units  
f_rb_constant <- function(., ... ) {
  # output in m2s mol-1 h2o
  
  .$pars$rb
}

# f_rb_leafdim <- function(., ... ) {
#   # output in s m-1 h2o
#   
#   ( .$env$lwind / .$pars$leaf_width )^-0.5 / .$pars$can_ttc
# }




### TEMPERATURE DEPENDENCE FUNCTIONS
################################
# - all these functions should be scalars such that the current temperature and reference temperature can be specified and the correction scalar is returned

f_scalar_none <- function(...){
  1
}

# temperature dependence functions that cannot be separtaed into ascending and decending components
f_temp_scalar_bethy <- function(.,parlist,...) { 
  
  exp(-(.$state$leaf_temp-parlist$Tr)/10) * 
    ( (parlist$a_q10_t + parlist$b_q10_t*.$state$leaf_temp) ^ ((parlist$a_q10_t + parlist$b_q10_t*.$state$leaf_temp)/(10*parlist$b_q10_t))    /  
        ( (parlist$a_q10_t + parlist$b_q10_t*parlist$Tr) ^ ((parlist$a_q10_t + parlist$b_q10_t*parlist$Tr)/(10*parlist$b_q10_t)) ) )
  
}


# Ascending components of the temperature response function - can be run alone for an increasing repsonse only
f_temp_scalar_Arrhenius <- function(.,parlist,...){
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # Arrhenius equation
  
  # input parameters  
  # Ha     -- rate of increase to optimum  (J mol-1)
  # R      -- molar gas constant J mol-1 K-1
  
  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  #convert to Kelvin
  Trk <- parlist$Tr + 273.15
  Tsk <- .$state$leaf_temp + 273.15
  
  if(.$cpars$verbose) {
    print('Arrhenius_tcorr_function')
    print(parlist)
  }
  
  exp( parlist$Ha*(Tsk-Trk) / (.$pars$R*Tsk*Trk) )
}

f_temp_scalar_Q10 <- function(.,parlist,...){
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  if(.$cpars$verbose) {
    print('Q10_tcorr_function')
    print(.$state$leaf_temp)
    print(parlist)
  }
  
  q10 <- get(parlist$q10_func)(.,parlist)
  
  parlist$q10 ^ ((.$state$leaf_temp-parlist$Tr)/10)
  
}


# Descending components of the temperature response function 
f_temp_scalar_modArrhenius_des <- function(.,parlist,...){
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # descending component of modified Arrhenius temperature response function, Medlyn et al 2002
  
  # input parameters  
  # Hd     -- rate of decrease from optimum (J mol-1) often 200000 
  # deltaS -- entropy parameter      
  # R      -- molar gas constant J mol-1 K-1
  
  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  if(.$cpars$verbose) {
    print('Arrhenius descending function')
    print(.$state$leaf_temp)
    print(parlist)
  }
  
  # convert to Kelvin
  Trk <- parlist$Tr + 273.15
  Tsk <- .$state$leaf_temp + 273.15
  
  deltaS <- get(.$fnames$deltaS)(.,parlist)
  
  (1 + exp((Trk*deltaS-parlist$Hd) / (Trk*.$pars$R)) ) / 
    (1 + exp((Tsk*deltaS-parlist$Hd) / (Tsk*.$pars$R)) )  
  
}

f_temp_scalar_collatz1991_des <- function(.,parlist,...){
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # descending component of temperature scaling from Collatz etal 1991
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  # convert to Kelvin
  Tsk <- .$state$leaf_temp + 273.15
  
  # get deltaS
  deltaS <- get(.$fnames$deltaS)(.,parlist)
  
  1 / ( 1 + exp((Tsk*deltaS-parlist$Hd) / (Tsk*.$pars$R)) )
}

f_temp_scalar_cox2001_des <- function(.,parlist,...){
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # descending component of temperature scaling from Cox etal 2001
  
  # input parameters  
  # Tr     -- reference temperature (oC) 
  
  1 / ( (1 + exp(parlist$exp*(.$state$leaf_temp-parlist$tupp))) * (1 + exp(parlist$exp*(parlist$tlow-.$state$leaf_temp))) )
}



# functions that can allow for temperature acclimation of parameters, deltaS, Q10

# deltaS
f_deltaS_constant <- function(.,parlist,...) {
  #constant delta S

  parlist$deltaS
}

f_deltaS <- function(.,parlist,...) {
  #calculate delta S from T opt (temp where t scaling peaks) in oC
  #Medlyn 2002
  
  Toptk  <- parlist$Topt + 273.15
  
  parlist$Hd/Toptk + (.$pars$R*log(parlist$Ha/(parlist$Hd - parlist$Ha)))
}

f_deltaS_lin_t <- function(.,parlist,...) {
  #calculate delta S

  # CLM limits the range of growth temps 
  parlist$a_deltaS_t + .$state$leaf_temp * parlist$b_deltaS_t 
}

# Q10
f_q10_constant <- function(.,parlist,...) {
  # constant delta S
  
  parlist$q10
}

f_q10_lin_t <- function(.,parlist,...) {
  # calculate q10 as a function of T
  # Tjoelker 
  
  parlist$a_q10_t + .$state$leaf_temp * parlist$b_q10_t 
}



### Gamma star - CO2 compensation point in the absence of dark respiration
f_gstar_constant <- function(.,...){
  .$pars$atref.gstar
}

f_gstar_f1980 <- function(.,...) {
  # calculates Gamma star as a function of Kc & Ko, and thus their combined temperature dependence
  # Farquhar 1980 Eq 38
  # 0.21 = ko/kc
  
  .$pars$ko_kc_ratio * .$state_pars$Kc*.$state$oi/(2*.$state_pars$Ko)
}

f_gstar_constref <- function(.) {
  # takes a defined ref temperature value of gstar and scales to leaf temp
  # this will probably not give the correct response to a change in atmospheric pressure
  
  .$pars$atref.gstar * get(.$fnames$gstar_tcor)(.,parlist=list(Tr=.$pars$reftemp.gstar,Ha=.$pars$Ha.gstar)) 
}

f_gstar_c1991 <- function(.) {
  # takes a defined ref temperature value of tau and scales to leaf temp
  # calcualtes gstar at leaftemp from tau
  
  .$state_pars$tau <- .$pars$atref.tau * get(.$fnames$tau_tcor)(.,parlist=list(Tr=.$pars$reftemp.tau,q10=.$pars$q10.tau,Ha=.$pars$Ha.tau,q10_func='f_q10_constant'))
  .$state$oi/(2*.$state_pars$tau)   
}

f_temp_scalar_quadratic_bf1985 <- function(.,parlist,...) {
  # calculates Gamma star (umol mol-1) temperature scalar (K or oC)
  # could be expanded to incorporate parameters in parlist, but currently gstar specific
  # Brooks&Farquhar 1985
  # rearranged to give a scalar of value at 25oC 

  1 + (.$pars$gstar_bf_b*(.$state$leaf_temp-parlist$Tr) + .$pars$gstar_bf_a*(.$state$leaf_temp-parlist$Tr)^2) / .$pars$gstar_bf_c
}





