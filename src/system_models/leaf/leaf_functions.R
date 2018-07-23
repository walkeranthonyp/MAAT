################################
#
# Leaf level physiology functions
# 
# AWalker March 2014
#
################################

source('../generic_mathematical_functions.R')



f_none <- function(.) {
  NA
}

### SOLVERS & RESIDUAL FUNCTIONS
################################

# Solver to find root of .$fnames$solver_func
f_R_Brent_solver <- function(.) {
  if(.$cpars$verbose_loop) print(.$env) 

  .$solver_out <- uniroot(get(.$fnames$solver_func),interval=c(-0.002765326,50.1234),.=.,extendInt='downX')
  .$solver_out$root
}


# Calculate assimilation for a given cc (.$state$cc)
# - code block common to all assimilation solvers
# - calculates Ag/cc, determines limiting rate, calculates and returns net A
f_assimilation <- function(.) {
 
  # calculate Ag / cc for each limiting process
  .$state$Acg     <- get(.$fnames$Acg)(.)
  .$state$Ajg     <- get(.$fnames$Ajg)(.)
  .$state$Apg     <- get(.$fnames$Apg)(.)
  
  # determine rate limiting cycle - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  Amin            <- get(.$fnames$Alim)(.) 
  
  # calculate & return net A
  Amin*.$state$cc - Amin*.$state_pars$gstar - .$state$rd
}
  

# Residual function for solver to calculate assimilation
f_A_r_leaf <- function(., A ) {
  # combines A, rs, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  #  -- for use with uniroot solver
  #  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation

  # calculate cc from ca, rb, rs, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer and stomatal resistance terms are in h2o units
  # assumes mesophyll resistance is in co2 units
  .$state$cc <- get(.$fnames$gas_diff)( . , A , r=( 1.4*.$state_pars$rb + 1.6*get(.$fnames$rs)(.,A=A,c=get(.$fnames$gas_diff)(.,A)) + .$state_pars$ri ) )
  
  #print(c(A,f_assimilation(.)))
  
  # calculate residual of net A
  f_assimilation(.) - A
} 


# same as above function but with no stomatal resistance 
f_A_r_leaf_noRs <- function(.,A) {
  
  # calculate cc from ca, rb, and ri
  # assumes boundary layer resistance is in h2o units
  # assumes mesophyll resistance is in co2 units
  .$state$cc <- get(.$fnames$gas_diff)( . , A , r=( 1.4*.$state_pars$rb + .$state_pars$ri ) )
  
  # calculate residual of net A
  f_assimilation(.) - A
}


# Semi-analytical solution 
f_A_r_leaf_semiana <- function(.) {
  # - finds the analytical solution assuming rb and ri are zero to use as first guess (a0)
  # - this should always be larger than the full solution (unless rb and ri are zero)
  # - make a second (a1) and third (a2) guess smaller than a0
  # - calculate solver function value for these three guesses (fa0, fa1, fa2) 
  # - check to make sure these span the root
  # - fit a quadratic through the threee sets of co-ordinates to find the root 
 
  # save values that may get reset by analytical function 
  rb_hold <- .$state_pars$rb
  g0_hold <- .$pars$g0
  ri_hold <- .$state_pars$ri
  
  # find the analytical solution assuming rb and ri are zero to use as first guess (a0)
  a0 <- .$state$A_ana_rbzero <- f_A_r_leaf_analytical_quad(.)

  # return values potentially reset by analytical function to the data structure
  .$state_pars$rb <- rb_hold
  .$state_pars$ri <- ri_hold
  .$pars$g0       <- g0_hold 

  # currently no method to handle a0 < 0 - this should be fine 
  #if(a0<0) stop('a0 < 0 ,',a0)

  # a0 should mostly be larger than the full solution (unless rb and ri are zero) meaning that fa0 should always be < 0
  fa0 <- f_A_r_leaf(., a0 ) 
 
  #  
  if(abs(fa0) < 1e-6 ) return(a0)
  else {
 
    # in rare case when initial guess is smaller than solution and fa0 > 0
#    if(fa0 > 1 ) { 
#      print(paste('initial fa > 0,',fa0,'; a,',a0))
#      a0  <- 2 * a0 
#      fa0 <- f_A_r_leaf(., a0 )
#      print(paste('new initial fa > 0,',fa0,'; a,',a0))
#    } else if(fa0 > 1e-6 ) { 
#      print(paste('initial fa > 0,',fa0,'; a,',a0))
#      # second guess, also analytical
#      .$state$cc <- get(.$fnames$gas_diff)( . , a0 , r=( 1.4*.$state_pars$rb + 1.6*get(.$fnames$rs)(.,A=a0,c=get(.$fnames$gas_diff)(.,a0)) + .$state_pars$ri ) )
#      a1         <- f_assimilation(.) 
#      fa1        <- get(.$fnames$solver_func)(., a1 ) 
#      print(c('analytical fa1,',fa1,'a,',a1))
#      sinput <- seq(a0-2,a0+2,0.1 )
#      out    <- numeric(length(sinput))
#      for( i in 1:length(sinput) ) {
#        out[i] <- f_A_r_leaf(., A=sinput[i] )
#      }
#      print(cbind(sinput,out))
#      print(unlist(.$env)); print(unlist(.$state_pars))  
#      a0  <- a0 + .$pars$deltaA_prop * a0
#      fa0 <- f_A_r_leaf(., a0 )
#      print(paste('new initial fa > 0,',fa0,'; a,',a0))
#      print('')
#    } 

#    # second guess
#    deltaA1 <- .$pars$deltaA_prop * a0 
#    a1      <- a0 - deltaA1
#    fa1     <- get(.$fnames$solver_func)(., a1 ) 
#
#    # third guess
#    deltaA2 <- if(fa1*fa0 < 0) 0.5 * deltaA1 else 2 * deltaA1   
#    a2      <- a0 - deltaA2
#    fa2     <- get(.$fnames$solver_func)(., a2 ) 

    # second guess, also analytical
    .$state$cc <- get(.$fnames$gas_diff)( . , a0 , r=( 1.4*.$state_pars$rb + 1.6*get(.$fnames$rs)(.,A=a0,c=get(.$fnames$gas_diff)(.,a0)) + .$state_pars$ri ) )
    a1         <- f_assimilation(.) 
    fa1        <- get(.$fnames$solver_func)(., a1 ) 

    # third guess
    a2 <- mean(c(a0,a1))  
    #if(fa1*fa0 < 0) a2 <- mean(c(a0,a1))  
    #else            a2 <- a1 - (a0-a1)   # works if initial guess is too high but what if initial guess is too small? need to change sign
    fa2 <- get(.$fnames$solver_func)(., a2 ) 

    # check f(guesses) span zero i.e. that guesses span the root
    guesses  <- c(a0,a1,a2)
    fguesses <- c(fa0,fa1,fa2)
    if( min(fguesses) * max(fguesses) >= 0 ) {
      print(unlist(.$fnames)); print(unlist(.$pars)); print(unlist(.$state_pars)); print(unlist(.$env)) 
      #print(paste('Solver error: fa0-2 all > or < 0,',fa0,fa1,fa2,'; guesses,',a0,a1,a2))
      sinput <- seq(a0-2,a0+2,0.1 )
      out    <- numeric(length(sinput))
      for( i in 1:length(sinput) ) {
        out[i] <- f_A_r_leaf(., A=sinput[i] )
      }
      print(cbind(sinput,out))
      stop(paste('Solver error: fa0-2 all > or < 0,',fa0,fa1,fa2,'; guesses,',a0,a1,a2))      
#      a2  <- a2 - deltaA1
#      fa2 <- get(.$fnames$solver_func)(., a2 ) 
#      guesses  <- c(a0,a1,a2)
#      fguesses <- c(fa0,fa1,fa2)
#      if( min(fguesses) * max(fguesses) >= 0 ) stop(paste('Solver error: fa0-2 all > or < 0,',fa0,fa1,fa2,'; guesses,',a0,a1,a2))
    } 
    .$state$aguess  <- guesses 
    .$state$faguess <- fguesses
   
    # fit a quadratic through the three sets of co-ordinates a la Lomas (one step of Muller's method) 
    # Muller, David E., "A Method for Solving Algebraic Equations Using an Automatic Computer," Mathematical Tables and Other Aids to Computation, 10 (1956)    
    bx <- ((a1+a0)*(fa2-fa1)/(a2-a1)-(a2+a1)*(fa1-fa0)/(a1-a0))/(a0-a2)
    ax <- (fa2-fa1-bx*(a2-a1))/(a2**2-a1**2)
    cx <- fa1-bx*a1-ax*a1**2
 
    # find the root of the quadratic that lies between guesses 
    assim <- .$state$assim <- quad_sol(ax,bx,cx,'both')
    .$state$fA_ana_final[1] <- get(.$fnames$solver_func)(., assim[1] ) 
    .$state$fA_ana_final[2] <- get(.$fnames$solver_func)(., assim[2] ) 
    ss    <- which(assim>min(guesses)&assim<max(guesses))
    
    # catch potential errors
    if(length(ss)==0)      stop('no solution within initial 3 guesses') 
    else if(length(ss)==2) stop('both solutions within initial 3 guesses') 
    else return(assim[ss])
  }
}



### ANALYTICAL SOLUTIONS
################################

# solves A analytically by assuming g0 = 0 in the stomatal resistance function, and rb and ri also = zero 
f_A_r_leaf_analytical <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes
  
  .$state_pars$rb <- 0
  .$state_pars$ri <- 0
  
  # calculate cb, ci & cc
  .$state$cb      <- .$state$ca
  fe              <- get(paste0(.$fnames$rs,'_fe') )(.) 
  .$state$ci      <- .$state$ca * (1 - (1.6 / fe) )
  .$state$cc      <- .$state$ci 
  
  # calculate net A
  Anet <- f_assimilation(.)
  
  # calculate rs
  .$state_pars$rs <- .$state$ca / (fe * Anet * .$env$atm_press*1e-6) 
  
  # return net A
  Anet
}


# solves A analytically by assuming rb and ri are zero 
f_A_r_leaf_analytical_quad <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes

  # set cb, rb & ri
  .$state$cb      <- .$state$ca
  .$state_pars$rb <- 0
  .$state_pars$ri <- 0

  # correctly assign g0 if rs functions assume g0 = 0
  g0_hold <- .$pars$g0
  if(.$fnames$rs=='f_rs_cox1998'|.$fnames$rs=='f_rs_constantCiCa') .$pars$g0[] <- 0 

  # calculate coefficients of quadratic to solve A
  assim_quad_soln <- function(., V, K ) {
    gsd <- get(paste0(.$fnames$rs,'_fe'))(.) / .$state$ca
    p   <- .$env$atm_press*1e-6
    a   <- p*( 1.6 - gsd*(.$state$ca + K) )
    b   <- p*gsd*( .$state$ca*(V - .$state$rd) - .$state$rd*K - V*.$state_pars$gstar ) - .$pars$g0*(.$state$ca + K) + 1.6*p*(.$state$rd - V)
    c   <- .$pars$g0*( V*(.$state$ca - .$state_pars$gstar) - .$state$rd*(K + .$state$ca) )
 
    # return A 
    quad_sol(a,b,c,'upper')
  }

  .$state$Acg <- assim_quad_soln(., V=.$state_pars$vcmaxlt, K=.$state_pars$Km )
  .$state$Ajg <- assim_quad_soln(., V=(.$state$J/4),        K=(2*.$state_pars$gstar) )
  .$state$Apg <- assim_quad_soln(., V=(3*.$state_pars$tpu), K=(-(1+3*.$pars$Apg_alpha)*.$state_pars$gstar) )

  # determine rate limiting cycle - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  Amin        <- get(.$fnames$Alim)(.) 
  
  # determine cc/ci based on Amin
  # calculate rs
  .$pars$g0[]     <- g0_hold 
  .$state_pars$rs <- get(.$fnames$rs)(.,A=Amin)
  .$state$cc <-.$state$ci <- f_ficks_ci(., A=Amin, r=1.6*.$state_pars$rs )
    
  # recalculate Ag for each limiting process
  # necessary if Alim is Collatz smoothing as it reduces A, decoupling A from cc calculated in the quadratic solution 
  .$state$Acg <- get(.$fnames$Acg)(.) * .$state$cc
  .$state$Ajg <- get(.$fnames$Ajg)(.) * .$state$cc
  .$state$Apg <- get(.$fnames$Apg)(.) * .$state$cc

  # return net A
  Amin
}


# solves A analytically by assuming rs is equal to 1/g0 
f_A_r0_leaf_analytical_quad <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes

  r0 <- get(paste0(.$fnames$rs,'_r0'))(.) 

  # calculate coefficients of quadratic to solve A
  assim_quad_soln <- function(.,V,K) {
    r   <- 1.4*.$state_pars$rb + 1.6*r0 + .$state_pars$ri
    p   <- .$env$atm_press*1e-6
    a   <- -p*r
    b   <- .$state$ca + K - .$state$rd*p*r + V*p*r
    c   <- .$state$ca*(.$state$rd-V) + .$state$rd*K + V*.$state_pars$gstar 
    
    # return A 
    quad_sol(a,b,c,'lower')
  }

  .$state$Acg <- assim_quad_soln(., V=.$state_pars$vcmaxlt, K=.$state_pars$Km )
  .$state$Ajg <- assim_quad_soln(., V=(.$state$J/4),        K=(2*.$state_pars$gstar) )
  .$state$Apg <- assim_quad_soln(., V=(3*.$state_pars$tpu), K=(-(1+3*.$pars$Apg_alpha)*.$state_pars$gstar) )

  # determine rate limiting cycle - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  Amin        <- get(.$fnames$Alim)(.) 
  
  # determine cc/ci based on Amin
  .$state_pars$rs <- r0 
  .$state$cb <- f_ficks_ci(., A=Amin )
  .$state$cc <-.$state$ci <- f_ficks_ci(., A=Amin, r=1.6*.$state_pars$rs )
    
  # recalculate Ag for each limiting process
  # necessary if Alim is Collatz smoothing as it reduces A, decoupling A from cc calculated in the quadratic solution 
  .$state$Acg <- get(.$fnames$Acg)(.) * .$state$cc
  .$state$Ajg <- get(.$fnames$Ajg)(.) * .$state$cc
  .$state$Apg <- get(.$fnames$Apg)(.) * .$state$cc

  # return net A
  Amin
}


# Calculate assimilation assuming zero resistance to CO2 diffusion from the atmosphere to the site of carboxylation
f_A_r_leaf_noR <- function(.,...) {
  # This function can be used to calculate the stomatal limitation to photosynthesis when rb and ri are assumed zero 
  # (when ri and rb are non-zero, use below function 'f_A_r_leaf_noRs' )
  
  # combines all rate limiting processes
  
  # assume cc = ca
  .$state$cc <- .$state$ca
  
  # calculate net A
  f_assimilation(.)
}



### PHOTOSYNTHESIS FUNCTIONS
################################

# Transition point function, calculates cc at which Acg = Ajg, i.e. the transition cc
transition_cc <- function(.) {
  
  vcm_et_ratio <- .$state_pars$vcmaxlt/.$state$J 
  
  ( 8*.$state_pars$gstar*vcm_et_ratio - .$state_pars$Km ) /
    (1 - 4*vcm_et_ratio)
}


# Carboxylation limitation
f_Acg_farquhar1980 <- function(., cc=.$state$cc ) {   
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
  .$state_pars$vcmaxlt /
    (cc + .$state_pars$Km)
}

# electron transport limitation
f_Ajg_generic <- function(., cc=.$state$cc ) { 
  # generic eq to calculate light limited photosynthetic rate 
  # currently no other formulation exists (other than slighly different parameters in the denominator)
  # umol m-2 s-1
  
  # calculate gross electron transport limited carboxylation rate / cc 
  .$state$J / (4*(cc+2*.$state_pars$gstar))     
}

# Farquhar 1980 and others
f_j_farquhar1980 <- function(.){
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per incident photon
  # this is the replacement for eq A2 added in the proofing stage (very last line of the manuscript after all refs etc)
  
  .$state_pars$jmaxlt * .$env$par * (1 - .$pars$f)/(.$env$par + 2.2*.$state_pars$jmaxlt)
}

# Harley etal 1992 eq to calculate electron transport rate
f_j_harley1992 <- function(.) {
  # umol e m-2 s-1
  harley_alpha <- .$pars$a * .$state_pars$alpha
  harley_alpha*.$env$par / ( 1.0 + (harley_alpha**2)*(.$env$par**2) / (.$state_pars$jmaxlt**2) )**0.5  
}

# von Caemmerer 2000 book and taken from Farquhar & Wong 1984
f_j_farquharwong1984 <- function(.) {
  # calculates J given Jmax, I - irradiance & alpha - electrons transported per photon
  
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

# Collatz etal 1991 eq to calculate electron transport rate
f_j_collatz1991 <- function(.) {
  # umol e m-2 s-1
  
  .$pars$a * .$state_pars$alpha * .$env$par  
}


# TPU limitation
# triose phosphate limitation from vonCaemmerer 2000 as corrected in Gu 2010
f_Apg_vonc2000 <- function(., cc=.$state$cc ) {
  # this is derived from Harley & Sharkey 1991
  # the difference is that Harley & Sharkey iteratively solved to find tpu, while here tpu is set independently
  
  # As described in Gu 2010, the collatz 1991 function for TPU limitation is a special case of this one where:
  # alpha = 0, tpu = vcmax/6
  # and tpu temperature scaling is identical to vcmax
  
  # calculate gross TPU limited carboxylation rate / cc 
  ifelse( cc <= (1+3*.$pars$Apg_alpha)*.$state_pars$gstar, NA,
          3*.$state_pars$tpu / ( cc-(1+3*.$pars$Apg_alpha)*.$state_pars$gstar ) 
          )
}

# triose phosphate limitation from Foley 1996, citing Harley & Sharkey 1991
f_Apg_foley1996 <- function(., cc=.$state$cc ) {
  # its not clear to me how they derive this from Harley and Sharkey
  # Eq 5 from Foley 1996: Js = 3 * tpu * (1 - gstar/cc) + Jp * gstar / cc
  # where Jp is the limiting rate of wc or wj
  
  if( cc <= .$state_pars$gstar ) NA
  else {
    # calculate limiting cycle of Acg and Ajg
    .$state$Apg <- NA
    amin <- get(.$fnames$Alim)(.) 

    (3*.$state_pars$tpult + amin) / cc
  }
}

# no triose phosphate limitation
f_Apg_none <- function(.) {
  # returns a high value so that TPU is never limiting
  
  9e3
}


# limiting rate selection 
# simple minimum of all three possible limitating states
f_lim_farquhar1980 <- function(.) {
  
  min(c(.$state$Acg,.$state$Ajg,.$state$Apg),na.rm=T)
}

# smoothed solution of all three possible limiting states
f_lim_collatz1991 <- function(.) {

  a  <- .$pars$theta_col_cj
  b  <- -1 * (.$state$Acg + .$state$Ajg)
  c  <- .$state$Acg * .$state$Ajg
  
  sol1 <- quad_sol(a,b,c)

  ifelse(!is.na(.$state$Apg), {
    a  <- .$pars$theta_col_cjp
    b  <- -1 * (.$state$Apg + sol1)
    c  <- .$state$Apg * sol1
    quad_sol(a,b,c)

  }, sol1)
}

 

#########################
### Respiration functions

f_rd_constant <- function(.) {
  .$pars$atref[['rd']]
}

f_rd_lin_vcmax <- function(.) {
  .$pars$a_rdv_25 + .$state_pars$vcmax * .$pars$b_rdv_25    
}

# rd25 is a proportion of vcmax 25 but that proportion changes with temperature
f_rd_lin_vcmax_t <- function(.) {
  .$pars$b_rdv_25 <- .$pars$a_rdv_25_t + .$state$leaf_temp * .$pars$b_rdv_25_t
  .$pars$a_rdv_25 + .$state_pars$vcmax * .$pars$b_rdv_25    
}

f_rd_lin_N <- function(.) {
  .$pars$a_rdn_25 + .$state$leafN_area * .$pars$b_rdn_25    
}

# light supression of respiration
# - return scalars of respiration in light : respiration in dark
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
  .$pars$atref[['vcmax']]
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
  .$pars$atref[['jmax']]
}

f_jmax_power <- function(.) {
  exp(.$pars$e_ajv_25) * .$state_pars$vcmax^.$pars$e_bjv_25    
}

#f_jmax_lin <- function(.) {
#  .$pars$ajv_25 + .$state_pars$vcmax * .$pars$bjv_25    
#}

#f_jmax_lin_t <- function(.) {
#  .$pars$ajv_25 <- 0.0
#  .$pars$bjv_25 <- .$pars$a_jvt_25 + .$state$leaf_temp * .$pars$b_jvt_25
#  .$pars$ajv_25 + .$state_pars$vcmax * .$pars$bjv_25    
#}

f_jmax_lin <- function(.) {
  .$pars$ajv_25 + .$state_pars$vcmax * .$pars$bjv_25 * get(.$fnames$tcor_jmax)(.)    
}

f_tcor_jmax_lin <- function(.) {
  (.$pars$a_jvt_25 + .$pars$b_jvt_25 * .$state$leaf_temp ) /
    (.$pars$a_jvt_25 + .$pars$b_jvt_25 * 25 )
}


# TPU
f_tpu_constant <- function(.) {
  .$pars$atref[['tpu']]
}

f_tpu_lin <- function(.) {
  .$pars$atv_25 + .$state_pars$vcmax * .$pars$btv_25    
}



### STOMATAL & RELATED CONDUCTANCE / RESISTANCE FUNCTIONS
################################

# CO2 diffusion
f_ficks_ci <- function(., A=.$state$A, r=1.4*.$state_pars$rb, c=.$state$ca ) {
  # can be used to calculate cc, ci or cb (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1
  
  c-A*r*.$env$atm_press*1e-6
}

f_ficks_ci_bound0 <- function(., ... ) {
  # can be used to calculate cc, ci or cs (boundary CO2 conc) from either ri, rs or rb respectively
  # by default calculates cb from ca and rb
  # c units in Pa
  # A units umol m-2 s-1
  # r units  m2 s mol-1
  
  # there is a catch 22 here, make this max of the function result or zero and it screws up the solver boundaries
  # remove the max term and this screws up at very low ca
  #   max( c-A*r*.$env$atm_press*1e-6 , 1e-6)
  c2 <- f_ficks_ci(.,...) 
  if(c2>0) c2 else 0
}


# general
f_conv_ms_molm2s1 <- function(.) {
  .$env$atm_press / ((.$state$leaf_temp + 273.15)*.$pars$R)
}

f_r_zero <- function(., ... ) {
  0
}

f_r_zero_fe <- function(., ... ) {
  0
}

f_r_zero_r0 <- function(., ... ) {
  0
}


# stomata
# stomatal resistances are all assumed to be in h2o units 
f_rs_constant <- function(., ... ) {
  # output in m2s mol-1 h2o

  # this currently doesn't work with the solver.
  # When rs is not a function of A the solver interval doesn't span zero,
  # if this is a function to be used must reroute solver to an analytical method

  .$pars$rs
}


# Medlyn et al 2011 eq for stomatal resistance
f_rs_medlyn2011 <- function(., A=.$state$A, c=.$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / (.$pars$g0 + f_rs_medlyn2011_fe(.) * A * .$env$atm_press*1e-6 / c) 
}

f_rs_medlyn2011_fe <- function(.) {
  # f(e) component of rs from Medlyn 2011   
  ( 1 + .$pars$g1_medlyn / .$env$vpd^0.5 )
}

f_rs_medlyn2011_r0 <- function(.) {
  # g0 component of rs from Medlyn 2011   
  1 / max(.$pars$g0,1e-6)
}


# Leuning et al 1995 eq for stomatal resistance
f_rs_leuning1995 <- function(., A=.$state$A, c=.$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1  h2o

  1 / (.$pars$g0 + f_rs_leuning1995_fe(.,c=c) * A * .$env$atm_press*1e-6 / c)  
}

f_rs_leuning1995_fe <- function(., c=.$state$cb ) {
  # f(e) component of rs from Leuning 1995   
  .$pars$g1_leuning / ( (1 - .$state_pars$gamma/c) * (1 + .$env$vpd/.$pars$d0) )  
}

f_rs_leuning1995_r0 <- function(.) {
  # g0 component of rs from Leuning 1995   
  1 / max(.$pars$g0,1e-6)
}


# Ball et al 1987 eq for stomatal resistance
f_rs_ball1987 <- function(., A=.$state$A, c=.$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / (.$pars$g0 + f_rs_ball1987_fe(.) * A * .$env$atm_press*1e-6 / c) 
}

f_rs_ball1987_fe <- function(.) {
  # f(e) component of rs from Ball 1987   
  .$pars$g1_ball * .$env$rh
}

f_rs_ball1987_r0 <- function(.) {
  # g0 component of rs from Ball 1987   
  1 / max(.$pars$g0,1e-6)
}


# implied LPJ etc assumption for stomatal resistance that keeps Ci:Ca constant
f_rs_constantCiCa <- function(., A=.$state$A, c=.$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o

  # set Ci:Ca ratio
  .$state_pars$cica_chi <- get(.$fnames$cica_ratio)(.)
  
  1 / (f_rs_constantCiCa_fe(.) * A * .$env$atm_press*1e-6 / c) 
}

f_rs_constantCiCa_fe <- function(.) {
  # f(e) component of rs for constant Ci:Ca   
  .$state_pars$cica_chi <- get(.$fnames$cica_ratio)(.)
  1.6 / (1 - .$state_pars$cica_chi)
}

f_rs_constantCiCa_r0 <- function(.) {
  # g0 component of rs for constant Ci:Ca
  # there is no g0 so an arbitrarily low number is used to avoid / 0   
  1 / 1e-6 
}

f_cica_constant <- function(.) {
  .$pars$cica_chi
}


# implied JULES etc assumption for stomatal resistance that keeps a variant of Ci:Ca constant
f_rs_cox1998 <- function(., A=.$state$A, c=.$state$cb ) {
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / (f_rs_cox1998_fe(.,c=c) * A * .$env$atm_press*1e-6 / c)
}

f_rs_cox1998_fe <- function(., c=.$state$cb ) {
  # f(e) component of rs for Cox 1998   
  
  f0    <- 1 - 1.6/.$pars$g1_leuning
  dstar <- (.$pars$g1_leuning/1.6 - 1) * .$pars$d0
  CmCP  <- (1 - .$state_pars$gamma/c)

  1.6 / ( CmCP - f0*CmCP * (1 - .$env$vpd/dstar))
}

f_rs_cox1998_r0 <- function(.) {
  # g0 component of rs for cox1998
  # there is no g0 so an arbitrarily low number is used to avoid / 0   
  1 / 1e-6 
}


# ORCHIDEE assumption for stomatal resistance, from Yin & Struik 2009
f_rs_yin2009 <- function(., A=.$state$A, c=.$state$ci ) {
  # This will not work with the either analytical solution (as they are currently coded) due to the A + Rd in the denominator
  # this also prevents negative values when A is negative 
  
  # expects c in Pa
  # output in m2s mol-1 h2o
  
  1 / (1.6*( .$pars$g0 + ( (A + .$state$rd)  / (1e6/.$env$atm_press*(.$state$ci-.$state_pars$gstar)) ) * (1/(1/(.$pars$g_a1_yin - .$pars$g_b1_yin*.$env$vpd) - 1)) ))
}

f_rs_yin2009_fe <- function(.) {
  # This will not work with the either analytical solution (as they are currently coded) due to the A + Rd in the denominator
  stop('Yin & Struik 2009 rs function not compatible with analytical solution')
}

f_rs_yin2009_r0 <- function(.) {
  # This function is designed to avoid the problem of negative rs when A is negative, could just call the native function
  1 / max(.$pars$g0,1e-6)
}


# internal/mesophyll
# internal/mesophyll resistances are all assumed by the solver to be in co2 units 
f_ri_constant <- function(., ... ) {
  # output in m2s mol-1 co2
  
  .$pars$ri
}

f_ri_sphag_wf1998 <- function(., W=.$state$fwdw_ratio ) {
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
  
  .$pars$rb
}

f_rb_leafdim <- function(., ... ) {
  # output in s mol-1m-2 h2o
  
  cf <- (.$pars$R * (.$state$leaf_temp+273.15) ) / .$env$atm_press
  .$pars$can_ttc^-1 * ( .$env$wind / .$pars$leaf_width )^-0.5 * cf  
}

f_rb_water_e2009 <- function(.){
  # calculates rb in m2s mol-1 h20 
  # where water level is below 0, returns zero resistance 
  # assumes water level is in mm
  
  # calculates rb in s m-1
  rb <- ((.$env$water_l-.$env$sphag_l)/1000) / (.$pars$co2_diff * (1 + .$pars$hco_co2_ratio * .$pars$hco_co2_diff_ratio) )   # both of these parameters are temp dependent    
  if((.$env$water_l-.$env$sphag_l)>0) {
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
  .$state_pars$vcmaxlt / .$state_pars$vcmax
}

# temperature scaling is independent of vcmax
f_tcor_dep_independent <- function(., var ) {

  get(.$fnames$tcor_asc[[var]])(., var=var ) * get(.$fnames$tcor_des[[var]])(., var=var )
}

# temperature dependence functions that cannot be separated into ascending and decending components
f_tcor_asc_bethy <- function(., var, ... ) { 
  #tcor_des <- .$fnames$tcor_des[[var]]
  if(.$fnames$tcor_des[[var]]!='f_scalar_none') print('Warning: temp scaling will not work correctly, f_tcor_asc_bethy specified with a descending temperature scaling other than f_scalar_none')
 
  exp(-(.$state$leaf_temp-.$pars$reftemp[[var]])/10) * 
    ( (.$pars$a_q10_t[[var]] + .$pars$b_q10_t[[var]]*.$state$leaf_temp) ^ ((.$pars$a_q10_t[[var]] + .$pars$b_q10_t[[var]]*.$state$leaf_temp)/(10*.$pars$b_q10_t[[var]]))    /  
        ( (.$pars$a_q10_t[[var]] + .$pars$b_q10_t[[var]]*.$pars$reftemp[[var]]) ^ ((.$pars$a_q10_t[[var]] + .$pars$b_q10_t[[var]]*.$pars$reftemp[[var]])/(10*.$pars$b_q10_t[[var]])) ) )
  
}

# calculates Gamma star (umol mol-1) temperature scalar (K or oC)
f_tcor_asc_quadratic_bf1985 <- function(., var, ... ) {
  # could be expanded to allow additional parameters to use this T scaling method but currently gstar specific 
  # Brooks&Farquhar 1985
  # rearranged to give a scalar of value at 25oC 

  1 + (.$pars$gstar_bf_b*(.$state$leaf_temp-.$pars$reftemp[[var]]) + .$pars$gstar_bf_a*(.$state$leaf_temp-.$pars$reftemp[[var]])^2) / .$pars$gstar_bf_c
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
  Trk <- .$pars$reftemp[[var]] + 273.15
  Tsk <- .$state$leaf_temp + 273.15
  
  exp( .$pars$Ha[[var]]*(Tsk-Trk) / (.$pars$R*Tsk*Trk) )
}

# Q10 temperature scaling
f_tcor_asc_Q10 <- function(., var, ... ) {
  #returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
 
  q10 <- get(.$fnames$q10_func[[var]])(., var=var )
 
  q10 ^ ((.$state$leaf_temp - .$pars$reftemp[[var]])/10)
  
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
  
  if(.$cpars$verbose) {
    print('Arrhenius descending function')
    print(.$state$leaf_temp)
    print(var)
  }
  
  # convert to Kelvin
  Trk <- .$pars$reftemp[[var]] + 273.15
  Tsk <- .$state$leaf_temp + 273.15
  
  deltaS <- get(.$fnames$deltaS[[var]])(., var )
  
  (1 + exp((Trk*deltaS-.$pars$Hd[[var]]) / (Trk*.$pars$R)) ) / 
    (1 + exp((Tsk*deltaS-.$pars$Hd[[var]]) / (Tsk*.$pars$R)) )  
  
}

# descending component of temperature scaling from Collatz etal 1991
f_tcor_des_collatz1991 <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Q10    -- factor by which rate increases for every 10 oC of temp increase  
  # Tr     -- reference temperature (oC) 
  
  # convert to Kelvin
  Tsk <- .$state$leaf_temp + 273.15
  
  # get deltaS
  deltaS <- get(.$fnames$deltaS[[var]])(., var )
  
  1 / ( 1 + exp((Tsk*deltaS-.$pars$Hd[[var]]) / (Tsk*.$pars$R)) )
}

# descending component of temperature scaling from Cox etal 2001
f_tcor_des_cox2001 <- function(., var, ... ) {
  # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
  
  # input parameters  
  # Tr     -- reference temperature (oC) 
  
  1 / ( (1 + exp(.$pars$exp_cox[[var]]*(.$state$leaf_temp-.$pars$tupp_cox[[var]]))) * 
        (1 + exp(.$pars$exp_cox[[var]]*(.$pars$tlow_cox[[var]]-.$state$leaf_temp))) )
}


# functions that can allow for temperature acclimation of parameters, deltaS, Q10
# deltaS
f_deltaS_constant <- function(., var, ... ) {
  #constant delta S

  .$pars$deltaS[[var]]
}

#calculate delta S from T opt (temp where t scaling peaks) in oC
f_deltaS <- function(., var, ... ) {
  #Medlyn 2002
  
  Toptk  <- .$pars$Topt[[var]] + 273.15
  
  .$pars$Hd[[var]]/Toptk + (.$pars$R*log(.$pars$Ha[[var]]/(.$pars$Hd[[var]] - .$pars$Ha[[var]])))
}

#calculate delta S as a function of growth temp
f_deltaS_lin_t <- function(., var, ... ) {

  # CLM limits the range of growth temps 
  .$pars$a_deltaS_t[[var]] + .$state$leaf_temp * .$pars$b_deltaS_t[[var]] 
}

# Q10
f_q10_constant <- function(., var, ... ) {
  # constant delta S
  
  .$pars$q10[[var]]
}

# calculate q10 as a function of T
f_q10_lin_t <- function(., var, ... ) {
  # Tjoelker 
  
  .$pars$a_q10_t[[var]] + .$state$leaf_temp * .$pars$b_q10_t[[var]] 
}


### Gamma star - CO2 compensation point in the absence of dark respiration
f_gstar_constant <- function(., ... ) {
  .$pars$atref[['gstar']]
}

# calculates Gamma star as a function of Kc & Ko, and thus their combined temperature dependence
f_gstar_f1980 <- function(., ... ) {
  # Farquhar 1980 Eq 38
  # 0.21 = ko/kc
  
  .$pars$ko_kc_ratio * .$state_pars$Kc*.$state$oi/(2*.$state_pars$Ko)
}

# takes a defined ref temperature value of gstar and scales to leaf temp
f_gstar_constref <- function(.) {
  # this will probably not give the correct response to a change in atmospheric pressure
  
  .$pars$atref[['gstar']] * get(.$fnames$tcor_asc[['gstar']])(., var='gstar' ) 
}

# calcualtes gstar at leaftemp from tau
f_gstar_c1991 <- function(.) {
  # takes a defined ref temperature value of tau and scales to leaf temp
  
  .$state_pars$tau <- .$pars$atref[['tau']] * get(.$fnames$tcor_asc[['tau']])(., var='tau', q10_func='f_q10_constant' )
  .$state$oi/(2*.$state_pars$tau)   
}



### END ###
