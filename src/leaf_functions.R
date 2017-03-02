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

quad_sol <- function(a,b,c,out='lower') {
  # robust numerical solution to the quadratic
  # taken from Numerical Recipes

  q     <- -0.5 * ( b + sign(b)*(b^2 - 4*a*c)^0.5 )
  roots <- c( q/a , c/q )
  
  if(out=='lower')      min(roots,na.rm=T) 
  else if(out=='upper') max(roots,na.rm=T)
  else roots 
}



### SOLVERS & SOLVER FUNCTIONS
################################

f_R_analytical <- function(.,v,k,r) {
  
  # num       <-  v * (.$state$cb - r - .$state_pars$gstar)
  # denom     <- .$state$cb - r + k
  # carboxylation gross of rd
  # num/denom
  # carboxylation net of rd
  # num/denom - .$state$respiration
  
  v / (.$state$cb - r + k)
}

f_A_r_leaf_analytical <- function(.) {
  # combines A, rs, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  # passes cc to the assimilation function
  
  # calculate cc from ca, rb, rs, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer and stomatal resistance terms are in h2o units
  # assumes mesophyll resistance is in co2 units

  # currently this function does not account for leaf boundary layer and internal resitance
  # this is because these terms sum to get overall resistance of CO2 to the chloroplast
  # and therefore the simplification in f_R_analytical does not apply
  # to include rb and ri a general analytical solution needs work, may not be possible
  rs_simple <- get(paste(.$fnames$rs,'fg1',sep='_') )(.) 
  r_simple  <- rs_simple * .$env$atm_press * 1.6e-6
  # r_simple <- rs_simple + .$state_pars$rb + .$state_pars$ri
  
  # calculate electron transport rate
  get(.$fnames$wj)(.)

  # calculate w
  # - this calculation of w does not use the same functions as the numerical solver function 
  # - wc and wj do not have multiple respective equations for their calculation
  # - wc and wj equations also have the same form, the f_R_analytical function exploits this similarity in form
  .$state$wc <- f_R_analytical(.,v=.$state_pars$vcmaxlt, k=.$state_pars$Km,      r=r_simple)
  .$state$wj <- f_R_analytical(.,v=.$state$J/4,          k=2*.$state_pars$gstar, r=r_simple)
  .$state$wp <- NULL
  
  # calculate limiting cycle
  wmin <- get(.$fnames$Alim)(.) 
  
  # calculate ci
  .$state$ci <- (.$state$cb - r_simple)
  
  # calculate actual w
  .$state$wc <- .$state$wc * .$state$ci 
  .$state$wj <- .$state$wj * .$state$ci 
  .$state$wp <- .$state$wp * .$state$ci
  
  # calculate net A
  Anet <- wmin*.$state$ci - wmin*.$state_pars$gstar - .$state$respiration

  # calculate rs
  .$state_pars$rs <- rs_simple / Anet

  # return net A
  Anet
}

f_R_Brent_solver <- function(.,...) {
  # ... could be used to pass through different functions 'func' to the solver function, not currently necessary

  if(.$cpars$verbose_loop) print(.$env)  
  .$solver_out <- uniroot(get(.$fnames$solver_func),interval=c(-100,100),.=.,extendInt='no',...)
  .$solver_out$root
}

f_A_r_leaf <- function(A,.,...) {
  # combines A, rs, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  # passes cc to the assimilation functions
  #  -- for use with uniroot solver
  #  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation

  # calculate cc from ca, rb, rs, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer and stomatal resistance terms are in h2o units
  # assumes mesophyll resistance is in co2 units
  .$state$cc <- get(.$fnames$gas_diff)( . , A , r=( 1.4*.$state_pars$rb + 1.6*get(.$fnames$rs)(.,A=A,c=get(.$fnames$gas_diff)(.,A)) + .$state_pars$ri ) )
  
  # calculate w/cc
  # - w/cc and not w is calculated for numerical stability at low cc
  .$state$wc <- get(.$fnames$wc)(.)
  .$state$wj <- get(.$fnames$wj)(.)
  .$state$wp <- get(.$fnames$wp)(.)

  # calculate limiting cycle
  wmin <- get(.$fnames$Alim)(.) 
  
  # calculate actual w
  .$state$wc <- .$state$wc * .$state$cc
  .$state$wj <- .$state$wj * .$state$cc
  .$state$wp <- .$state$wp * .$state$cc
  
  # calculate net A
  wmin*.$state$cc - wmin*.$state_pars$gstar - .$state$respiration - A
} 

f_A_r_leaf_noRs <- function(A,.,...) {
  # same as above function but with no stomatal resistance 
  
  # combines A, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  # passes cc to the assimilation functions
  #  -- for use with uniroot solver
  #  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation
  
  # calculate cc from ca, rb, rs, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer and stomatal resistance terms are in h2o units
  # assumes mesophyll resistance is in co2 units
  cc <- get(.$fnames$gas_diff)( . , A , r=( 1.4*.$state_pars$rb + .$state_pars$ri ) )
  
  # calculate w/cc
  # - w/cc and not w is calculated for numerical stability at low cc
  .$state$wc <- get(.$fnames$wc)(.,cc=cc)
  .$state$wj <- get(.$fnames$wj)(.,cc=cc)
  .$state$wp <- get(.$fnames$wp)(.,cc=cc)
  
  # calculate limiting cycle
  wmin <- get(.$fnames$Alim)(.) 
  
  # calculate actual w
  .$state$wc <- .$state$wc * cc
  .$state$wj <- .$state$wj * cc
  .$state$wp <- .$state$wp * cc
  
  # calculate net A
  wmin*cc - wmin*.$state_pars$gstar - .$state$respiration - A
}

f_A_r_leaf_noR <- function(.,...) {
  # same as above function but with no resistance to CO2 diffusion
  # - from the atmosphere to the site of carboxylation
  # These functions can be used to calculate the stomatal limitation to photosynthesis  
  
  # combines all rate limiting processes  
  # cc is initialised at ca by the leaf run function 
  # so if this is run prior to solver for A then cc is already set to ca
  
  # calculate w/cc
  # - w/cc and not w is calculated for numerical stability at low cc
  .$state$wc <- get(.$fnames$wc)(.)
  .$state$wj <- get(.$fnames$wj)(.)
  .$state$wp <- get(.$fnames$wp)(.)
  
  # calculate limiting cycle
  wmin <- get(.$fnames$Alim)(.) 
  
  # calculate net A
  wmin*.$state$cc - wmin*.$state_pars$gstar - .$state$respiration
}

transition_cc <- function(.) {
  # calculates the cc at which wc = wj, i.e. the transition cc
  
  vcm_et_ratio <- .$state_pars$vcmaxlt/.$state$J 
  
  (8*.$state_pars$gstar*vcm_et_ratio - (.$state_pars$Kc*(1+(.$state$oi/.$state_pars$Ko))) ) /
    (1 - 4*vcm_et_ratio)
}


### PHOTOSYNTHESIS FUNCTIONS
################################

# CO2 diffusion
f_ficks_ci <- function(.,A=.$state$A,r=.$state_pars$rb,c=.$state$ca) {
  # can be used to calculate cc, ci or cs (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1
  
  c-A*r*.$env$atm_press*1e-6
#   c-A*r*.$env$atm_press*1e-6
}

f_ficks_ci_bound0 <- function(.,A=.$state$A,r=.$state_pars$rb,c=.$state$ca) {
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


# Carboxylation limitation
f_wc_farquhar1980 <- function(.,cc=.$state$cc){   
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
#   .$state_pars$vcmaxlt*cc /
#     (cc+.$state_pars$Kc*(1+(.$state$oi/.$state_pars$Ko)))
  .$state_pars$vcmaxlt /
    (cc + .$state_pars$Km)
}


# electron transport limitation
f_wj_generic <- function(.,cc=.$state$cc){
  # generic eq to calculate light limited photosynthetic rate from user selected electron transport and etoc function
  # umol m-2 s-1
  
  # calculate electron transport rate
  .$state$J <- get(.$fnames$etrans)(.)
  
  # calculate gross electron transport limited carboxylation 
  # i.e. rd and G* not included
  #   (.$state$J/4) * (cc/(cc+2*.$state_pars$gstar)) # currently no other formualtion exists (other than slighly different parameters in the denominator)
  .$state$J / (4*(cc+2*.$state_pars$gstar))        # currently no other formualtion exists (other than slighly different parameters in the denominator)
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
  a  <- .$pars$theta
  b  <- -1 * (I2 + .$state_pars$jmaxlt)
  c  <- I2 * .$state_pars$jmaxlt
  
  #(b - (b^2 - 4*a*c)^0.5) / (2*a) # this probably needs to be more robust
  quad_sol(a,b,c)
}

f_j_collatz1991 <- function(.){
  # Collatz etal 1991 eq to calculate electron transport rate
  # umol e m-2 s-1
  
  .$pars$a * .$state_pars$alpha * .$env$par  
}


# TPU limitation
# f_wp_collatz1991 <- function(.,cc=.$state$cc){
#   # triose phosphate limitation from collatz 1991
#   # the /cc is to make this function compatible with the solver and the form of these equations
#   
# #   .$state_pars$vcmaxlt/2
#   # could also read .$state_pars$TPU/cc and set TPU when setting Vcmax etc
#   .$state_pars$vcmaxlt/2/cc
# }

f_wp_vonc2000 <- function(.,cc=.$state$cc){
  # triose phosphate limitation from vonCaemmerer 2000 as corrected in Gu 2010
  # this is derived from Harley & Sharkey 1991
  # the difference is that Harley & Sharkey iteratively solved to find tpu, while here tpu is set independently
  
  # As described in Gu 2010, the collatz 1991 function for TPU limitation is a special case of this one where:
  # alpha = 0, tpu = vcmax/6
  # and tpu temperature scaling is identical to vcmax
  
  if( .$state$cc <= (1+3*.$pars$wp_alpha)*.$state_pars$gstar) NA
#   3*.$state_pars$tpu*cc / ( cc-(1+3*.$pars$wp_alpha)*.$state_pars$gstar )
  else 3*.$state_pars$tpu / ( cc-(1+3*.$pars$wp_alpha)*.$state_pars$gstar )
}

f_wp_foley1996 <- function(.,cc=.$state$cc){
  # triose phosphate limitation from Foley 1996, citing Harlkey & Sharkey 1991
  # its not clear to me how they derive this from Harley and Sharkey
  # Eq 5 from Foley 1996: Js = 3 * tpu * (1 - gstar/cc) + Jp * gstar / cc
  # where Jp is the limiting rate of wc or wj
  
  if( .$state$cc <= .$state_pars$gstar ) NA
  else {
    # calculate limiting cycle of wc and wj
    .$state$wp <- NA
    wmin <- get(.$fnames$Alim)(.) 

    (3*.$state_pars$tpu + wmin) / cc
  }
}


# limiting factor
f_lim_farquhar1980 <- function(.){
  # simple minimum of all three possible limitating states
  
  min(c(.$state$wc,.$state$wj,.$state$wp),na.rm=T)
}

f_lim_collatz1991 <- function(.){
  # smoothed solution of all three possible limiting states

  a  <- .$pars$theta_collatz
  b  <- -1 * (.$state$wc + .$state$wj)
  c  <- .$state$wc * .$state$wj
  
  #sol1 <- (b - (b^2 - 4*a*c)^0.5) / (2*a) # this probably needs to be more robust
  sol1 <- quad_sol(a,b,c)

  if(!is.na(.$state$wp)) {
    a  <- .$pars$beta_collatz
    b  <- -1 * (.$state$wp + sol1)
    c  <- .$state$wp * sol1
    
    #(b - (b^2 - 4*a*c)^0.5) / (2*a) # this probably needs to be more robust    
    quad_sol(a,b,c)
  } else sol1
}

 

#########################
### Respiration functions

f_rd_constant <- function(.) {
  .$pars$atref.rd
}

f_rd_farquhar1980 <- function(.) {
  # Farquhar etal 1980 dark respiration rate umol m-2 s-1
  # @ 25oC 23 Pa Ci & 1000 umol m-2 s-1 PAR
  
  # could make this fixed but have rd_prop_vcmax multiplied by 100 say 
  # to be able to vary the constant with a parameter 
  1.1
}

f_rd_collatz1991 <- function(.) {
  .$pars$rd_prop_vcmax * .$state_pars$vcmaxlt
}

f_rd_fN <- function(.) {
  .$pars$rd_prop_N * .$state$leafN_area
}

# light supression



### PHOTOSYNTHESIS PARAMETER FUNCTIONS
################################

# vcmax
f_constant_vcmax <- function(.) {
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
f_constant_jmax <- function(.) {
  .$pars$atref.jmax
}

f_jmax_walker2014 <- function(.) {
  exp(.$pars$e_ajv_25) * .$state_pars$vcmax^.$pars$e_bjv_25    
}

f_jmax_lin <- function(.) {
  .$pars$ajv_25 + .$state_pars$vcmax25 * .$pars$bjv_25    
}

# TPU
f_constant_tpu <- function(.) {
  .$pars$atref.tpu
}

f_tpu_lin <- function(.) {
  .$pars$atv_25 + .$state_pars$vcmax25 * .$pars$btv_25    
}



### STOMATAL & RELATED CONDUCTANCE / RESISTANCE FUNCTIONS
################################

# general
f_conv_ms_molm2s1 <- function(.){
  .$env$atm_press / ((.$state$leaf_temp + 273.15)*.$pars$R)
}

f_r_zero <- function(.,...){
  0
}

f_r_zero_fcb <- function(.,...){
  0
}

f_r_zero_fg1 <- function(.,...){
  0
}


# stomata
# stomatal resistances are all assumed by the solver to be in h2o units 
f_rs_constant <- function(.,...) {
  # output in m2s mol-1 h2o 
  
  .$pars$rs
}

f_rs_medlyn2011 <- function(.,A=.$state$A,c=.$state$cb){
  # Medlyn et al 2011 eq for stomatal resistance
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  # if( A < 0 ) 1/.$pars$g0
  # else 1 / (.$pars$g0 + (1 + .$pars$g1_medlyn/.$env$vpd^0.5) * A / (c/(.$env$atm_press * 1e-6)) )
  1 / (.$pars$g0 + (1 + .$pars$g1_medlyn/.$env$vpd^0.5) * A / (c/(.$env$atm_press * 1e-6)) )
}

f_rs_medlyn2011_fg1 <- function(.,c=.$state$cb) {
  # simplified version of rs   
  c / (.$env$atm_press * 1e-6) / ( 1 + .$pars$g1_medlyn / .$env$vpd^0.5 )
}


f_rs_leuning1995 <- function(.,A=.$state$A,c=.$state$cb){
  # Leuning et al 1995 eq for stomatal resistance
  # expects c in Pa
  # output in m2s mol-1  h2o

  # if( A < 0 ) 1/.$pars$g0
  # else 1 / (.$pars$g0 + .$pars$g1_leuning/(1+.$env$vpd/.$pars$d0) * A / ( (c-.$state_pars$gamma)/(.$env$atm_press * 1e-6) - .$state$respiration) )
  1 / ( .$pars$g0 + .$pars$g1_leuning * A / ( (1+.$env$vpd/.$pars$d0) * (c-.$state_pars$gamma)/(.$env$atm_press * 1e-6)) ) 
}

f_rs_leuning1995_fg1 <- function(.,c=.$state$cb) {
  # simplified version of rs   
  ( (c - .$state_pars$gamma)/(.$env$atm_press * 1e-6) * (1+.$env$vpd/.$pars$d0) ) / .$pars$g1_leuning 
}


f_rs_ball1987 <- function(.,A=.$state$A,c=.$state$cb){
  # Ball et al 1987 eq for stomatal resistance
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  # if( A < 0 ) 1/.$pars$g0
  # else 1 / (.$pars$g0 + .$pars$g1_ball*.$env$rh * A / (c/(.$env$atm_press * 1e-6)) )
  1 / ( .$pars$g0 + .$pars$g1_ball*.$env$rh*A / (c/(.$env$atm_press * 1e-6)) )
}

f_rs_ball1987_fg1 <- function(.,c=.$state$cb) {
  # simplified version of rs   
  c / (.$env$atm_press * 1e-6) / (.$pars$g1_ball*.$env$rh)
}


f_rs_constantCiCa <- function(.,A=.$state$A,c=.$state$cb) {
  # implied LPJ etc assumption for stomatal resistance that keeps Ci:Ca constant
  # expects c in Pa
  # output in m2s mol-1 h2o

  # set Ci:Ca ratio
  .$state_pars$cica_chi <- get(.$fnames$cica_ratio)(.)
  
  # if( A < 0 ) 1/1e-9
  # else ( c * (1 - .$state_pars$cica_chi) ) / ( .$env$atm_press*1.6e-6*A ) 
  ( c * (1 - .$state_pars$cica_chi) ) / ( .$env$atm_press*1.6e-6*A )
}

f_rs_constantCiCa_fg1 <- function(.,c=.$state$cb) {
  # simplified version of rs   
  .$state_pars$cica_chi <- get(.$fnames$cica_ratio)(.)
  ( c * (1 - .$state_pars$cica_chi) ) / ( .$env$atm_press*1.6e-6 )
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
  CmCP  <- (c - .$state_pars$gamma)
  
  # if( A < 0 ) 1/1e-9
  # else ( CmCP - f0*CmCP * (1 - .$env$vpd/dstar)) / ( .$env$atm_press*1.6e-6*A )
  ( CmCP - f0*CmCP * (1 - .$env$vpd/dstar)) / ( .$env$atm_press*1.6e-6*A )
}

f_rs_cox1998_fg1 <- function(.,c=.$state$cb) {
  # simplified version of rs   
  f0    <- 1 - 1.6/.$pars$g1_leuning
  dstar <- (.$pars$g1_leuning/1.6 - 1) * .$pars$d0
  CmCP  <- (c - .$state_pars$gamma)
  
  ( CmCP - f0*CmCP * (1 - .$env$vpd/dstar)) / (.$env$atm_press * 1.6e-6)
}



# internal/mesophyll
# internal/mesophyll resistances are all assumed by the solver to be in co2 units 

f_ri_constant <- function(.,...){
  # output in m2s mol-1 co2
  
  .$pars$ri
}

# leaf boundary layer
# leaf boundary layer resistances are all assumed by the solver to be in h2o units  



### TEMPERATURE DEPENDENCE FUNCTIONS
################################
# - all these functions should be scalars such that the current temperature and reference temperature can be specified and the correction scalar is returned
#   (this has not yet been completed)

f_temp_scalar_no_response <- function(...){
  1
}


# Arrhenius temperature response function 
f_temp_scalar_Arrhenius <- function(.,parlist,Tr=25,...){
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  # Arrhenius equation
  
  # input parameters  
  # Ha     -- rate of increase to optimum  (J mol-1)
  # R      -- molar gas constant J mol-1 K-1
  
  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  #convert to Kelvin
  Trk <- Tr + 273.15
  Tsk <- .$state$leaf_temp + 273.15
  
  if(.$cpars$verbose) {
    print('Arrhenius_tcorr_function')
    print(paste(Trk,Tsk))
    print(parlist)
  }
  
  exp( parlist$Ha*(Tsk-Trk) / (.$pars$R*Tsk*Trk) )
}


# modified Arrhenius temperature response function 
f_temp_scalar_modArrhenius <- function(.,parlist,Tr=25,...){
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  #Medlyn et al 2002
  
  # input parameters  
  # Ha     -- rate of increase to optimum  (J mol-1)
  # Hd     -- rate of decrease from optimum (J mol-1) often 200000 
  # deltaS --      
  # R      -- molar gas constant J mol-1 K-1

  # Tr     -- reference temperature (oC) 
  # Trk    -- reference temperature (K) 
  # Tsk    -- temperature to adjust parameter to (K) 
  
  #convert to Kelvin
  Trk <- Tr + 273.15
  Tsk <- .$state$leaf_temp + 273.15
  
  if(.$cpars$verbose) {
    print('Medlyn_tcorr_function')
    print(paste(Trk,Tsk))
    print(parlist)
  }
  
  deltaS <- f_deltaS(.,parlist)
  
  exp(parlist$Ha*(Tsk-Trk) / (.$pars$R*Tsk*Trk)) * ( 
    (1 + exp((Trk*deltaS-parlist$Hd) / (Trk*.$pars$R)) ) 
    / (1 + exp((Tsk*deltaS-parlist$Hd) / (Tsk*.$pars$R)) ) )  

}


# Maximum rate parameters - Vcmax & Jmax 
f_deltaS <- function(.,parlist,...){
  #calculate delta S
  #Medlyn 2002
  
  Toptk  <- parlist$Topt + 273.15
  
  parlist$Hd/Toptk + (.$pars$R*log(parlist$Ha/(parlist$Hd - parlist$Ha)))
}

f_temp_scalar_Q10 <- function(.,parlist,Tr=25,...){
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  #convert to Kelvin
  #   Trk <- Tr + 273.15
  #   Tsk <- .$state$leaf_temp + 273.15
  Ts <- .$state$leaf_temp
  
  if(.$cpars$verbose) {
    print('Q10_tcorr_function')
    print(paste(Tr,Ts))
    print(parlist)
  }
  
  parlist$q10 ^ ((Ts-Tr)/10)
  
}

f_temp_scalar_Q10_collatz1991 <- function(.,parlist,Tr=25,...){
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  #convert to Kelvin
  Tsk <- .$state$leaf_temp + 273.15
  
  # if(.$cpars$verbose) {
  #   print('Q10_collatz1991_function')
  #   print(paste(Tr,Ts))
  #   print(parlist)
  # }
  
  f_temp_scalar_Q10(.,parlist,Tr=25,...) + exp((Tsk*parlist$deltaS-parlist$Hd) / (Tsk*.$pars$R))
  
}

f_temp_scalar_Q10_cox1991 <- function(.,parlist,Tr=25,...){
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  # if(.$cpars$verbose) {
  #   print('Q10_cox2001_function')
  #   print(paste(Tr,Ts))
  #   print(parlist)
  # }
  
  f_temp_scalar_Q10(.,parlist,Tr=25,...) * 1 / ( (1 + exp(0.3*(.$state$leaf_temp-parlist$tupp))) * (1 + exp(0.3*(parlist$tlow-.$state$leaf_temp))) ) 
}


### Gamma star - CO2 compensation point in the absence of dark respiration
f_constant_gstar <- function(.,...){
  .$pars$atref.gstar
}

f_gstar_f1980 <- function(.,Tr=25,...){
  # calculates Gamma star as a function of Kc & Ko, and thus their combined temperature dependence
  # Farquhar 1980 Eq 38
  # 0.21 = ko/kc
  
  0.21 * .$state_pars$Kc*.$state$oi/(2*.$state_pars$Ko)
}

f_gstar_constref <- function(.) {
  # takes a defined ref temperature value of gstar and scales to leaf temp
  # this will probably not give the correct response to a change in atmospheric pressure
  
  .$state_pars$gstar <- .$pars$atref.gstar * get(.$fnames$gstar_tcor)(.,parlist=list(Ha=.$pars$Ha.gstar)) 
}

f_gstar_c1991 <- function(.) {
  # takes a defined ref temperature value of tau and scales to leaf temp
  # calcualtes gstar at leaftemp from tau
  
  .$state_pars$tau   <- .$pars$atref.tau * get(.$fnames$tau_tcor)(.,parlist=list(q10=.$pars$q10.tau))
  .$state_pars$gstar <- .$state$oi/(2*.$state_pars$tau)   
}

f_temp_scalar_quadratic_bf1985 <- function(.,Tr=25,...){
  # calculates Gamma star (umol mol-1) as a function of temperature difference (K or oC)
  # Brooks&Farquhar 1985
  # rearranged to give a scalar of value at 25oC 
  # calculated by B&F1985 from Jordan & Ogren 1984
  #   1 + (1.88*(.$state$leaf_temp-Tr) + 0.036*(.$state$leaf_temp-Tr)^2) / 44.7
  
  # B&F1985  
#   (42.7 + 1.68*(.$state$leaf_temp-Tr) + 0.012*(.$state$leaf_temp-Tr)^2)
#   1 + (1.68*(.$state$leaf_temp-Tr) + 0.012*(.$state$leaf_temp-Tr)^2) / 42.7
  1 + (.$pars$gstar_bf_b*(.$state$leaf_temp-Tr) + .$pars$gstar_bf_a*(.$state$leaf_temp-Tr)^2) / .$pars$gstar_bf_c
}

# f_temp_scalar_quadratic_jo1984 <- function(.,Tr=25,...){
#   # calculates Gamma star (umol mol-1) as a function of temperature difference (K or oC)
#   # Brooks&Farquhar 1985
#   # rearranged to give a scalar of value at 25oC 
#   
#   # calculated by B&F1985 from Jordan & Ogren 1984
# #   (44.7 + 1.88*(.$state$leaf_temp-Tr) + 0.036*(.$state$leaf_temp-Tr)^2) 
# #   1 + (1.88*(.$state$leaf_temp-Tr) + 0.036*(.$state$leaf_temp-Tr)^2) / 44.7
#   1 + (.$pars$gstar_jo_b*(.$state$leaf_temp-Tr) + .$pars$gstar_jo_a*(.$state$leaf_temp-Tr)^2) / .$pars$gstar_jo_c
#   
# }




