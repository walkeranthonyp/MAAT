################################
#
# Leaf functions
# 
# AWalker March 2014
#
################################

source('../generic_mathematical_functions.R')



f_none <- function(.) {
  NA
}

### PHOTOSYNTHESIS FUNCTIONS
################################

# Calculate assimilation for a given cc (.super$state$cc)
# - code block common to all C3 assimilation solvers
# - calculates Ag/cc, determines limiting rate, calculates and returns net A
f_assimilation_c3 <- function(.) {
 
  # calculate Ag / cc for each limiting process
  .super$state$Acg[] <- .$Acg()
  .super$state$Ajg[] <- .$Ajg()
  .super$state$Apg[] <- .$Apg()
  
  # determine rate limiting cycle - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  Amin <- .$Alim() 
  
  # calculate & return net A
  Amin*.super$state$cc - Amin*.super$state_pars$gstar - .super$state$rd
}
  
  
# C4 assimilation function   
f_assimilation_c4 <- function(.) {
 
  # calculate Ag for each limiting process
  .super$state$Acg[] <- .$Acg()
  .super$state$Ajg[] <- .$Ajg()
  .super$state$Apg[] <- .$Apg()
  
  # determine rate limiting cycle - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  Amin <- .$Alim() 
  
  # calculate & return net A
  Amin - .super$state$rd
}
  

# Transition point function, calculates cc at which Acg = Ajg, i.e. the transition cc
transition_cc <- function(.) {
  
  vcm_et_ratio <- .super$state_pars$vcmaxlt/.super$state$J 
  
  ( 8*.super$state_pars$gstar*vcm_et_ratio - .super$state_pars$Km ) /
    (1 - 4*vcm_et_ratio)
}


# Carboxylation limitation
f_Acg_farquhar1980 <- function(., cc=.super$state$cc ) {   
  # Farquhar 1980 eq to calculate RuBisCO limited photosynthetic rate
  # umol m-2 s-1
  
  # dvars
  # Ci (or Cc but vcmax needs adjusting)   
  
  # params
  # Vcmax maximum carboxylation rate of RuBisCO  (umol m-2 s-1)
  # Ko Michaelis constant for O2 mbar            (kPa)
  # Kc Michaelis constant for CO2                ( Pa) 
  
  # Ko and O must be in same units
  # Kc, Ci and Gamma must be in same units
  
  # calculate gross CO2 limited carboxylation rate / cc 
  .super$state_pars$vcmaxlt /
    (cc + .super$state_pars$Km)
}

# Carboxylation limitaion for C4 photosynthesis
f_Acg_c4_collatz1992 <- function(., ... ) {   
  .super$state_pars$vcmaxlt 
}


# electron transport limitation
f_Ajg_generic <- function(., cc=.super$state$cc ) { 
  # generic eq to calculate light limited photosynthetic rate 
  # currently no other formulation exists (other than slighly different parameters in the denominator)
  # umol m-2 s-1
  
  # calculate gross electron transport limited carboxylation rate / cc 
  .super$state$J / ( .super$pars$quantum_yield_to_eff * (cc+2*.super$state_pars$gstar) )     
}

# converts electron transport to gross carboxylation reactions
f_Ajg_c4_collatz1992 <- function(., ... ) { 
  .super$state$J / .super$pars$quantum_yield_to_eff 
}


# Farquhar 1980 and others
f_etrans_farquhar1980 <- function(.){
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per incident photon
  # this is the replacement for eq A2 added in the proofing stage (very last line of the manuscript after all refs etc)
  
  .super$state_pars$jmaxlt * .super$env$par * (1 - .super$pars$f)/(.super$env$par + 2.2*.super$state_pars$jmaxlt)
}

# Harley etal 1992 eq to calculate electron transport rate
f_etrans_harley1992 <- function(.) {
  # umol e m-2 s-1
  harley_alpha <- .super$pars$a * .super$state_pars$alpha
  harley_alpha*.super$env$par / ( 1.0 + (harley_alpha**2)*(.super$env$par**2) / (.super$state_pars$jmaxlt**2) )**0.5  
}

# von Caemmerer 2000 book and taken from Farquhar & Wong 1984
f_etrans_farquharwong1984 <- function(.) {
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per photon
  
  # in von Caemmerer 2000 - I2 = I * a*(1-f)/2
  # where a is leaf light absorptance, 1-f corrects for spectral quality and light not absorbed by photosystems, and /2 accounts for the 2 photons needed to fully tranport 1 electron 
  # this is presumably the same as I * alpha
  # if so von C 2000's alpha is 0.36125
    
  I2 <- .super$env$par * .super$pars$a * .super$state_pars$alpha    
  a  <- .super$pars$theta_j
  b  <- -1 * (I2 + .super$state_pars$jmaxlt)
  c  <- I2 * .super$state_pars$jmaxlt
  
  quad_sol(a,b,c)
}

# Collatz etal 1991 eq to calculate electron transport rate
f_etrans_collatz1991 <- function(.) {
  # umol e m-2 s-1
  
  .super$pars$a * .super$state_pars$alpha * .super$env$par  
}


# TPU limitation
# triose phosphate limitation from vonCaemmerer 2000 as corrected in Gu 2010
f_Apg_vonc2000 <- function(., cc=.super$state$cc ) {
  # this is derived from Harley & Sharkey 1991
  # the difference is that Harley & Sharkey iteratively solved to find tpu, while here tpu is set independently
  
  # As described in Gu 2010, the collatz 1991 function for TPU limitation is a special case of this one where:
  # alpha = 0, tpu = vcmax/6
  # and tpu temperature scaling is identical to vcmax
  
  # calculate gross TPU limited carboxylation rate / cc 
  ifelse( cc <= (1+3*.super$pars$Apg_alpha)*.super$state_pars$gstar, NA,
          3*.super$state_pars$tpu / ( cc-(1+3*.super$pars$Apg_alpha)*.super$state_pars$gstar ) 
          )
}

# triose phosphate limitation from Foley 1996, citing Harley & Sharkey 1991
f_Apg_foley1996 <- function(., cc=.super$state$cc ) {
  # its not clear to me how they derive this from Harley and Sharkey
  # Eq 5 from Foley 1996: Js = 3 * tpu * (1 - gstar/cc) + Jp * gstar / cc
  # where Jp is the limiting rate of wc or wj
  
  if( cc <= .super$state_pars$gstar ) NA
  else {
    # calculate limiting cycle of Acg and Ajg
    .super$state$Apg <- NA
    amin <- .$Alim(.) 

    (3*.super$state_pars$tpult + amin) / cc
  }
}

# no triose phosphate limitation
f_Apg_none <- function(.) {
  # returns NA which is ignored in limiting rate selection functions 
  NA
}

# PEPC limitation for C4 (i.e. not TPU limitation but it is the third possible limitation for C4) 
# - also, strictly speaking Cc should be Ci but this is to work with the solver
# - 
f_Apg_c4_pepc_collatz1992 <- function(., cc=.super$state$cc ) {
  cc * .super$state_pars$k_pepc_lt / .super$env$atm_press
}

# PEPC limitation for C4 - solution when rs = r0 
f_Apg_c4_pepc_collatz1992_r0soln <- function(.,r) {
  .super$state$ca * .super$state_pars$k_pepc_lt / (.super$env$atm_press*(1 + r*.super$state_pars$k_pepc_lt*1e-6))
}


# limiting rate selection 
# simple minimum of all three possible limitating states
f_Alim_farquhar1980 <- function(.) {
 
  min(c(.super$state$Acg,.super$state$Ajg,.super$state$Apg),na.rm=T)
}

# smoothed solution of all three possible limiting states
f_Alim_collatz1991 <- function(.) {

  a  <- .super$pars$theta_col_cj
  b  <- -1 * (.super$state$Acg + .super$state$Ajg)
  c  <- .super$state$Acg * .super$state$Ajg
  
  sol1 <- quad_sol(a,b,c)

  ifelse(!is.na(.super$state$Apg), {
    a  <- .super$pars$theta_col_cjp
    b  <- -1 * (.super$state$Apg + sol1)
    c  <- .super$state$Apg * sol1
    quad_sol(a,b,c)

  }, sol1)
}

 

#########################
### Respiration functions

f_rd_constant <- function(.) {
  .super$pars$atref$rd
}

f_rd_lin_vcmax <- function(.) {
  .super$pars$a_rdv_25 + .super$state_pars$vcmax * .super$pars$b_rdv_25    
}

# rd25 is a proportion of vcmax 25 but that proportion changes with temperature
f_rd_lin_vcmax_t <- function(.) {
  .super$pars$b_rdv_25 <- .super$pars$a_rdv_25_t + .super$state$leaf_temp * .super$pars$b_rdv_25_t
  .super$pars$a_rdv_25 + .super$state_pars$vcmax * .super$pars$b_rdv_25    
}

f_rd_lin_N <- function(.) {
  .super$pars$a_rdn_25 + .super$state$leafN_area * .super$pars$b_rdn_25    
}

# light supression of respiration
# - return scalars of respiration in light : respiration in dark
f_rl_rd_fixed <- function(.) {
  .super$pars$rl_rd_ratio
}

f_rl_rd_lloyd1995 <- function(.) {
  
  if(.super$env$par <= 10) 1
  else                ( .super$pars$rl_rd_lloyd_a - .super$pars$rl_rd_lloyd_b * ln(.super$env$par) )
}



### PHOTOSYNTHESIS PARAMETER FUNCTIONS
################################

# vcmax
f_vcmax_constant <- function(.) {
  .super$pars$atref$vcmax
}

f_vcmax_lin <- function(.) {
  if(.super$cpars$verbose) print('Vcmax_CLM_function')
  .super$pars$avn_25 + .super$state$leafN_area * .super$pars$bvn_25    
}

f_vcmax_clm <- function(.) {
 
  .super$state$leafN_area * .super$pars$flnr * .super$pars$fnr * .super$pars$Rsa  
}


# jmax
f_jmax_constant <- function(.) {
  .super$pars$atref$jmax
}

f_jmax_power <- function(.) {
  exp(.super$pars$e_ajv_25) * .super$state_pars$vcmax^.super$pars$e_bjv_25    
}

f_jmax_lin <- function(.) {
  .super$pars$ajv_25 + .super$state_pars$vcmax * .super$pars$bjv_25 * .$tcor_jmax(.)    
}

# scale JV ratio by temperature, returns a scalar
f_tcor_jmax_lin <- function(.) {
  (.super$pars$a_jvt_25 + .super$pars$b_jvt_25 * .super$state$leaf_temp ) /
    (.super$pars$a_jvt_25 + .super$pars$b_jvt_25 * 25 )
}

# Kumarathunge et al 2019 New Phytologist
f_tcor_jmax_lin_thome <- function(., var, ... ) {

  # currently just setting home temp to also be leaf temp 
  #.super$pars$home_temp[] <- .super$state$leaf_temp 
 
  print('')
  print(var)
  
  # CLM limits the range of growth temps 
  ( .super$pars$a_jvt_25 + .super$pars$home_temp * .super$pars$b_jvt_25 + 
    (.super$state$leaf_temp - .super$pars$home_temp ) * .super$pars$c_jvt_25)  / 
    (.super$pars$a_jvt_25 + 25*.super$pars$b_jvt_25 ) 
 
}


# TPU
f_tpu_constant <- function(.) {
  .super$pars$atref$tpu
}

f_tpu_lin <- function(.) {
  .super$pars$atv_25 + .super$state_pars$vcmax * .super$pars$btv_25    
}


# k PEPC  
f_k_pepc_constant <- function(.) {
  .super$pars$atref$k_pepc
}

f_k_pepc_lin <- function(.) {
  .super$pars$akv_25 + .super$state_pars$vcmax * .super$pars$bkv_25    
}



### STOMATAL & RELATED CONDUCTANCE / RESISTANCE FUNCTIONS
################################

# CO2 diffusion
f_gas_diff_ficks_ci <- function(., A=.super$state$A, r=1.4*.super$state_pars$rb, c=.super$state$ca ) {
  # can be used to calculate cc, ci or cb (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1

  c-A*r*.super$env$atm_press*1e-6
}

f_gas_diff_ficks_ci_bound0 <- function(., A=.super$state$A, r=1.4*.super$state_pars$rb, c=.super$state$ca ) {
  # can be used to calculate cc, ci or cs (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1
  
  # there is a catch 22 here, make this max of the function result or zero and it screws up the solver boundaries
  # remove the max term and this screws up at very low ca
  c2 <- c-A*r*.super$env$atm_press*1e-6
  if(c2>0) c2 else 0
}


# general
f_conv_ms_molm2s1 <- function(.) {
  .super$env$atm_press / ((.super$state$leaf_temp + 273.15)*.super$pars$R)
}

f_r_zero <- function(., ... ) {
  0
}

f_rs_constant_fe <- f_r_zero_fe <- function(., ... ) {
  0
}

f_r_zero_r0 <- function(., ... ) {
  0
}

# 13C discrimination
f_d13c_lin <- function(.) {
  #ci_conc <- 1e6 * .super$state$ci / .super$env$atm_press
  #.super$pars$d13c_a + (.super$pars$d13c_b_prime - .super$pars$d13c_a) * ci_conc / .super$env$ca_conc
  .super$pars$d13c_a + (.super$pars$d13c_b_prime - .super$pars$d13c_a) * (.super$state$ci / .super$state$ca)
}

f_d13c_classical <- function(.) {
#  ci_conc <- 1e6 * .super$state$ci / .super$env$atm_press
#  cc_conc <- 1e6 * .super$state$cc / .super$env$atm_press
#  gstar_conc <- 1e6 * .super$state_pars$gstar / .super$env$atm_press
#
#  .super$pars$d13c_a * ( .super$env$ca_conc - ci_conc ) / .super$env$ca_conc +
#    .super$pars$d13c_am * ( ci_conc - cc_conc ) / .super$env$ca_conc +
#    .super$pars$d13c_b  * cc_conc / .super$env$ca_conc -
#    .super$pars$d13c_f  * gstar_conc / .super$env$ca_conc
  .super$pars$d13c_a * (( .super$state$ca - .super$state$ci ) / .super$state$ca ) +
    .super$pars$d13c_am * (( .super$state$ci - .super$state$cc ) / .super$state$ca ) +
    .super$pars$d13c_b  * ( .super$state$cc / .super$state$ca ) -
    .super$pars$d13c_f  * ( .super$state_pars$gstar / .super$state$ca )
}


# stomata
# stomatal resistances are all assumed to be in h2o units 
f_rs_constant_r0 <- f_rs_constant <- function(., ... ) {
  # output in m2s mol-1 h2o

  # this currently doesn't work with the solver.
  # When rs is not a function of A the solver interval doesn't span zero,
  # if this is a function to be used must reroute solver to an analytical method

  .super$pars$rs
}


# Medlyn et al 2011 eq for stomatal resistance
# - with Medlyn 2012 correction - 1.6 in the f(e) term
f_rs_medlyn2011 <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / ( .super$pars$g0 + .$rs_fe()*A*.super$env$atm_press*1e-6/c ) 
}

f_rs_medlyn2011_fe <- function(.) {
  # f(e) component of rs from Medlyn 2011   
  1.6 * ( 1 + .super$pars$g1_medlyn/(.$env$vpd^0.5) )
}

f_rs_medlyn2011_r0 <- function(.) {
  # g0 component of rs from Medlyn 2011   
  1 / max(.super$pars$g0,1e-6)
}

# as Medlyn 2011 but with no g0
f_rs_medlyn_lin2015 <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / ( .$rs_fe()*A*.super$env$atm_press*1e-6/c ) 
}

f_rs_medlyn_lin2015_fe <- f_rs_medlyn2011_fe
f_rs_medlyn_lin2015_r0 <- f_rs_medlyn2011_r0

f_rs_dewar2019 <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / ( .super$pars$g0 + .$rs_fe()*A*.super$env$atm_press*1e-6/(c-.super$state_pars$gstar) ) 
}

f_rs_dewar2019_gamma <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / ( .super$pars$g0 + .$rs_fe()*A*.super$env$atm_press*1e-6/(c-.super$state_pars$gamma) ) 
}


f_rs_dewar2019_fe  <- f_rs_dewar2019_gamma_fe  <- f_rs_medlyn2011_fe
f_rs_dewar2019_r0  <- f_rs_dewar2019_gamma_r0  <- f_rs_medlyn2011_r0


# Leuning et al 1995 eq for stomatal resistance
f_rs_leuning1995 <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1  h2o

  1 / (.super$pars$g0 + .$rs_fe(c=c)*A*.super$env$atm_press*1e-6/c )  
}

f_rs_leuning1995_fe <- function(., c=.super$state$cb ) {
  # f(e) component of rs from Leuning 1995  

  .super$pars$g1_leuning / ( (1 - .super$state_pars$gamma/c) * (1 + .super$env$vpd/.super$pars$d0) )  
}

f_rs_leuning1995_r0 <- function(.) {
  # g0 component of rs from Leuning 1995   
  1 / max(.super$pars$g0,1e-6)
}


# Ball et al 1987 eq for stomatal resistance
f_rs_ball1987 <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / ( .super$pars$g0 + .$rs_fe()*A*.super$env$atm_press*1e-6/c ) 
}

f_rs_ball1987_fe <- function(.) {
  # f(e) component of rs from Ball 1987   
  .super$pars$g1_ball * .super$env$rh
}

f_rs_ball1987_r0 <- function(.) {
  # g0 component of rs from Ball 1987   
  1 / max(.super$pars$g0,1e-6)
}


# implied LPJ etc assumption for stomatal resistance that keeps Ci:Ca constant
f_rs_constantCiCa <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o

  # set Ci:Ca ratio
  .super$state_pars$cica_chi <- .$cica_ratio()
  
  1 / ( .$rs_fe()*A*.super$env$atm_press*1e-6/c ) 
}

f_rs_constantCiCa_fe <- function(.) {
  # f(e) component of rs for constant Ci:Ca   
  .super$state_pars$cica_chi <- .$cica_ratio()
  1.6 / (1 - .super$state_pars$cica_chi)
}

f_rs_constantCiCa_r0 <- function(.) {
  # g0 component of rs for constant Ci:Ca
  # there is no g0 so an arbitrarily low number is used to avoid / 0   
  1 / 1e-6 
}

f_cica_ratio_constant <- function(.) {
  .super$pars$cica_chi
}


# implied JULES etc assumption for stomatal resistance that keeps a variant of Ci:Ca constant
f_rs_cox1998 <- function(., A=.super$state$A, c=.super$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / ( .$rs_fe(c=c)*A*.super$env$atm_press*1e-6/c )
}

f_rs_cox1998_fe <- function(., c=.super$state$cb ) {
  # f(e) component of rs for Cox 1998   
  
  f0    <- 1 - 1.6/.super$pars$g1_leuning
  dstar <- (.super$pars$g1_leuning/1.6 - 1) * .super$pars$d0
  CmCP  <- (1 - .super$state_pars$gamma/c)

  1.6 / ( CmCP - f0*CmCP * (1 - .super$env$vpd/dstar))
}

f_rs_cox1998_r0 <- function(.) {
  # g0 component of rs for cox1998
  # there is no g0 so an arbitrarily low number is used to avoid / 0   
  1 / 1e-6 
}


# ORCHIDEE assumption for stomatal resistance, from Yin & Struik 2009
f_rs_yin2009 <- function(., A=.super$state$A, c=.super$state$ci ) {
  # This will not work with the either analytical solution (as they are currently coded) due to the A + Rd in the denominator
  # this also prevents negative values when A is negative 
  
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / (1.6*( .super$pars$g0 + ( (A + .super$state$rd)  / (1e6/.super$env$atm_press*(.super$state$ci-.super$state_pars$gstar)) ) * (1/(1/(.super$pars$g_a1_yin - .super$pars$g_b1_yin*.super$env$vpd) - 1)) ))
}

f_rs_yin2009_fe <- function(.) {
  # This will not work with the either analytical solution (as they are currently coded) due to the A + Rd in the denominator
  stop('Yin & Struik 2009 rs function not compatible with analytical solution')
}

f_rs_yin2009_r0 <- function(.) {
  # This function is designed to avoid the problem of negative rs when A is negative, could just call the native function
  1 / max(.super$pars$g0,1e-6)
}


# internal/mesophyll
# internal/mesophyll resistances are all assumed by the solver to be in co2 units 
f_ri_constant <- function(., ... ) {
  # output in m2s mol-1 co2
  
  .super$pars$ri
}

f_ri_sphag_wf1998 <- function(., W=.super$state$fwdw_ratio ) {
  # Flannagan & Williams 1998 internal resistance of Sphagnum
  # W is water content ratio (fresh weight / dry weight)
  # output in  m2s mol-1 
  # restricts max fwdw is 17 as higher than this would extrapolate beyond the calibration range of the empirical function and all goes weird 
  
  if(W>17) W <- 17
  1 / max( 1e-6 , (-0.195 + 0.1345*W - 0.02565*W^2 + 0.002276*W^3 - 0.00009842*W^4 + 0.000001677*W^5)  )
}



# leaf boundary layer
# leaf boundary layer resistances are all assumed by the solver to be in h2o units  
f_rb_constant <- function(., ... ) {
  # output in m2s mol-1 h2o
  
  .super$pars$rb
}

f_rb_leafdim <- function(., ... ) {
  # output in s mol-1m-2 h2o
  
  cf <- (.super$pars$R * (.super$state$leaf_temp+273.15) ) / .super$env$atm_press
  .super$pars$can_ttc^-1 * ( .super$env$wind / .super$pars$leaf_width )^-0.5 * cf  
}

f_rb_water_e2009 <- function(.){
  # calculates rb in m2s mol-1 h20 
  # where water level is below 0, returns zero resistance 
  # assumes water level is in mm
  
  # calculates rb in s m-1
  rb <- ((.super$env$water_l-.super$env$sphag_l)/1000) / (.super$pars$co2_diff * (1 + .super$pars$hco_co2_ratio * .super$pars$hco_co2_diff_ratio) )   # both of these parameters are temp dependent    
  if((.super$env$water_l-.super$env$sphag_l)>0) {
    # convert to resistance m2s mol-1 h2o
    #  (not sure h2o is correct for this diffusion though, but it converts to an atmospheric h2o transport equivalent relative to co2 so that it will work in the solver) 
    rb / f_conv_ms_molm2s1(.) / 1.4
  } else f_r_zero(.)  
}



### TEMPERATURE DEPENDENCE FUNCTIONS
################################
# - all these functions should be scalars such that the current temperature and reference temperature can be specified and the correction scalar is returned

f_scalar_none <- function(...) {
  1
}


# temperature scaling is identical to that of Vcmax
 f_tcor_dep_dependent<- function(., ... ) {
  .super$state_pars$vcmaxlt / .super$state_pars$vcmax
}

# temperature scaling is independent of vcmax
f_tcor_dep_independent <- function(., var ) {

  .[[paste('tcor_asc',var,sep='.')]](., var=var ) * .[[paste('tcor_des',var,sep='.')]](., var=var )
}


# temperature dependence functions that cannot be separated into ascending and decending components
f_tcor_asc_bethy <- function(., var, ... ) { 
  #tcor_des <- .super$fnames$tcor_des[[var]]
  if(.super$fnames$tcor_des[[var]]!='f_scalar_none') print('Warning: temp scaling will not work correctly, f_tcor_asc_bethy specified with a descending temperature scaling other than f_scalar_none')
 
  exp(-(.super$state$leaf_temp-.super$pars$reftemp[[var]])/10) * 
    ( (.super$pars$a_q10_t[[var]] + .super$pars$b_q10_t[[var]]*.super$state$leaf_temp) ^ ((.super$pars$a_q10_t[[var]] + .super$pars$b_q10_t[[var]]*.super$state$leaf_temp)/(10*.super$pars$b_q10_t[[var]]))    /  
        ( (.super$pars$a_q10_t[[var]] + .super$pars$b_q10_t[[var]]*.super$pars$reftemp[[var]]) ^ ((.super$pars$a_q10_t[[var]] + .super$pars$b_q10_t[[var]]*.super$pars$reftemp[[var]])/(10*.super$pars$b_q10_t[[var]])) ) )
  
}

# calculates Gamma star (umol mol-1) temperature scalar (K or oC)
f_tcor_asc_quadratic_bf1985 <- function(., var, ... ) {
  # could be expanded to allow additional parameters to use this T scaling method but currently gstar specific 
  # Brooks&Farquhar 1985
  # rearranged to give a scalar of value at 25oC 

  1 + (.super$pars$gstar_bf_b*(.super$state$leaf_temp-.super$pars$reftemp[[var]]) + .super$pars$gstar_bf_a*(.super$state$leaf_temp-.super$pars$reftemp[[var]])^2) / .super$pars$gstar_bf_c
}


# Ascending components of the temperature response function - can be run alone for an increasing repsonse only
f_tcor_asc_Arrhenius <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # Arrhenius equation
  
  # input parameters  
  # Ha     -- rate of increase to optimum  (J mol-1)
  # R      -- molar gas constant J mol-1 K-1
  
  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  #convert to Kelvin
  Trk <- .super$pars$reftemp[[var]] + 273.15
  Tsk <- .super$state$leaf_temp + 273.15
  
  exp( .super$pars$Ha[[var]]*(Tsk-Trk) / (.super$pars$R*Tsk*Trk) )
}


# Q10 temperature scaling
f_tcor_asc_Q10 <- function(., var, ... ) {
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
 
  q10 <- .[[paste('q10',var,sep='.')]](., var=var )
 
  q10 ^ ((.super$state$leaf_temp - .super$pars$reftemp[[var]])/10)
  
}


# Descending components of the temperature response function 
f_tcor_des_modArrhenius <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # descending component of modified Arrhenius temperature response function, Medlyn et al 2002
  
  # input parameters  
  # Hd     -- rate of decrease from optimum (J mol-1) often 200000 
  # deltaS -- entropy parameter      
  # R      -- molar gas constant J mol-1 K-1
  
  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  if(.super$cpars$verbose) {
    print('Arrhenius descending function')
    print(.super$state$leaf_temp)
    print(var)
  }
  
  # convert to Kelvin
  Trk <- .super$pars$reftemp[[var]] + 273.15
  Tsk <- .super$state$leaf_temp + 273.15
 
  # calculate deltaS 
  deltaS <- .[[paste('deltaS',var,sep='.')]](var=var)
  
  (1 + exp((Trk*deltaS-.super$pars$Hd[[var]]) / (Trk*.super$pars$R)) ) / 
    (1 + exp((Tsk*deltaS-.super$pars$Hd[[var]]) / (Tsk*.super$pars$R)) )  
  
}


# descending component of temperature scaling from Collatz etal 1991
f_tcor_des_collatz1991 <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  # convert to Kelvin
  Tsk <- .super$state$leaf_temp + 273.15
  
  # calculate deltaS 
  deltaS <- .[[paste('deltaS',var,sep='.')]](., var=var )
  
  1 / ( 1 + exp((Tsk*deltaS-.super$pars$Hd[[var]]) / (Tsk*.super$pars$R)) )
}


# descending component of temperature scaling from Collatz etal 1991
f_tcor_des_collatz1992 <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 

  1 / (1 + exp(.super$pars$exp_cox[[var]]*(.super$state$leaf_temp-.super$pars$tupp_cox[[var]])))  
}


# descending component of temperature scaling from Cox etal 2001
f_tcor_des_cox2001 <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Tr     -- reference temperature (oC) 
  
  1 / ( (1 + exp(.super$pars$exp_cox[[var]]*(.super$state$leaf_temp-.super$pars$tupp_cox[[var]]))) * 
        (1 + exp(.super$pars$exp_cox[[var]]*(.super$pars$tlow_cox[[var]]-.super$state$leaf_temp))) )
}


# functions that can allow for temperature acclimation of parameters, deltaS, Q10
# deltaS
f_deltaS_constant <- function(., var, ... ) {
  #constant delta S

  .super$pars$deltaS[[var]]
}

# calculate delta S from T opt (temp where t scaling peaks) in oC
f_deltaS <- function(., var, ... ) {
  #Medlyn 2002
  
  Toptk  <- .super$pars$Topt[[var]] + 273.15
  
  .super$pars$Hd[[var]]/Toptk + (.super$pars$R*log(.super$pars$Ha[[var]]/(.super$pars$Hd[[var]] - .super$pars$Ha[[var]])))
}

# calculate deltaS as a function of growth temp
f_Ea_lin_t <- function(., var, ... ) {

  # CLM limits the range of growth temps 
  .super$pars$a_Ea_t[[var]] + .super$state$leaf_temp * .super$pars$b_Ea_t[[var]] 
}

# calculate deltaS as a function of growth temp
f_deltaS_lin_t <- function(., var, ... ) {

  # CLM limits the range of growth temps 
  .super$pars$a_deltaS_t[[var]] + .super$state$leaf_temp * .super$pars$b_deltaS_t[[var]] 
}

# calculate deltaS as a function of growth temp and 'home' temp
# Kumarathunge et al 2019 New Phytologist
f_deltaS_lin_thome <- function(., var, ... ) {

  # currently just setting home temp to also be leaf temp 
  #.super$pars$home_temp[] <- .super$state$leaf_temp 
  
  # CLM limits the range of growth temps 
  .super$pars$a_deltaS_t[[var]] + .super$pars$home_temp * .super$pars$b_deltaS_t[[var]] + 
    (.super$state$leaf_temp - .super$pars$home_temp ) * .super$pars$c_deltaS_t[[var]]  
}


# Q10
f_q10_constant <- function(., var, ... ) {
  # constant delta S
  
  .super$pars$q10[[var]]
}


# calculate q10 as a function of T
f_q10_lin_t <- function(., var, ... ) {
  # Tjoelker 
  
  .super$pars$a_q10_t[[var]] + .super$state$leaf_temp * .super$pars$b_q10_t[[var]] 
}


### Gamma star - CO2 compensation point in the absence of dark respiration
f_gstar_constant <- function(., ... ) {
  .super$pars$atref$gstar
}


# calculates Gamma star as a function of Kc & Ko, and thus their combined temperature dependence
f_gstar_f1980 <- function(., ... ) {
  # Farquhar 1980 Eq 38
  # 0.21 = ko/kc
  
  .super$pars$ko_kc_ratio * .super$state_pars$Kc*.super$state$oi/(2*.super$state_pars$Ko)
}

# calculates Gamma star as a function of Kc & Ko, and thus their combined temperature dependence
f_photorespiration_f1980 <- function(.) {

  # based on last line of C3 assimilation function, Farquhar 1980 Eq 1
  Amin <- .$Alim() 
  Amin - .super$state$A - .super$state$rd
}


# takes a defined ref temperature value of gstar and scales to leaf temp
f_gstar_constref <- function(.) {
  # this will probably not give the correct response to a change in atmospheric pressure
  
  .super$pars$atref$gstar * .[['tcor_asc.gstar']](var='gstar') 
}


# calculates gstar at leaftemp from tau
f_gstar_c1991 <- function(.) {
  # takes a defined ref temperature value of tau and scales to leaf temp
  
  .super$state_pars$tau <- .super$pars$atref$tau * .[['tcor_asc.tau']](., var='tau', q10='f_q10_constant' )
  .super$state$oi/(2*.super$state_pars$tau)   
}



### END ###
