################################
#
# MAAT soil_decomp model functions
# - calculates parameters and specifies ODE system that is unique to each soil decomposition model 
# 
# Matt Craig, Anthony Walker, October 2020
#
################################

#function for converting ANPP to litter inputs from MIMICS
f_input_mimics <- function(., t ) {
  EST_LIT_in = .super$env$anpp / (365*24) # gC/m2/h (from gC/m2/y)
  EST_LIT    = EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h)
  EST_LIT/.super$env$depth 
  
}

############## CORPSE #################
# lsoda style function to solve CORPSE
#######################################
# - parms is a dummy argument to work with lsoda
f_solver_func_corpse <- function(., t, y, parms) {
  dCs <- .$input(t)[[1]] + .$desorp.ds5(t=t,C=y,i=5) - .$sorp.s1(t=t,C=y,i=1)*.$scor(.) - .$decomp.d1(t = t, C=y, i=1)*.$wcor(.)*.$tcor.t1(i=1)
  dCr <- .$input(t)[[2]] + .$desorp.ds6(t=t,C=y,i=6) - .$sorp.s2(t=t,C=y,i=2)*.$scor(.) - .$decomp.d2(t = t, C=y, i=2)*.$wcor(.)*.$tcor.t2(i=2)
  dCn <- .$decomp.d4(t = t, C=y, i=4)*.super$pars$cue[[4]] + .$desorp.ds7(t=t,C=y,i=7) - .$sorp.s3(t=t,C=y,i=3)*.$scor(.) - .$decomp.d3(t = t, C=y, i=3)*.$wcor(.)*.$tcor.t3(i=3)
  dM <- .$decomp.d1(t = t, C=y, i=1)*.$wcor(.)*.$tcor.t1(i=1)*.super$pars$cue[[1]]+.$decomp.d2(t = t, C=y, i=2)*.$wcor(.)*.$tcor.t2(i=2)*.super$pars$cue[[2]]+.$decomp.d3(t = t, C=y, i=3)*.$wcor(.)*.$tcor.t3(i=3)*.super$pars$cue[[3]] - .$decomp.d4(t = t, C=y, i=4)
  dPs <- .$sorp.s1(t=t,C=y,i=1)*.$scor(.) - .$desorp.ds5(t=t,C=y,i=5)
  dPr <- .$sorp.s2(t=t,C=y,i=2)*.$scor(.) - .$desorp.ds6(t=t,C=y,i=6)
  dPn <- .$sorp.s3(t=t,C=y,i=3)*.$scor(.) - .$desorp.ds7(t=t,C=y,i=7)
  list(c(dCs, dCr, dCn, dM, dPs, dPr, dPn))
}


############## MIMICS #################
# lsoda style function to solve MIMICS
#######################################
# - parms is a dummy argument to work with lsoda
f_solver_func_mimics <- function(., t, y, parms){
  
  .super$pars$vmax[[6]] = .super$pars$vmax[[2]]
  .super$pars$vmax2[[6]] = .super$pars$vmax2[[2]]
  
  #dynamic parameters
  fmet = .super$pars$mimics[['fmet_p1']] * (.super$pars$mimics[['fmet_p2']] - .super$pars$mimics[['fmet_p3']]*(.super$env$lignin/.super$env$N))
  
  
  # ensures that tau_mod1 is between two values
  # in Will's script, anpp is multipled by 0 in the manipulation scirpt... which would imply that 
  # tau_mod1 might always be set to 0.6 in the simulations in Ben's paper
  tau_mod1 = min(max(sqrt(.super$env$anpp/.super$pars$mimics[['tau_mod1_p1']]),.super$pars$mimics[['tau_mod1_p2']]),.super$pars$mimics[['tau_mod1_p3']]) 
  tau_r = .super$pars$mimics[['tau_r_p1']] * exp(.super$pars$mimics[['tau_r_p2']] * fmet) *tau_mod1 * .super$pars$mimics[['tau_mod2']]
  tau_k = .super$pars$mimics[['tau_k_p1']] * exp(.super$pars$mimics[['tau_k_p2']] * fmet) *tau_mod1 * .super$pars$mimics[['tau_mod2']]
  
  .super$pars$k[[3]] = tau_r
  .super$pars$k[[4]] = tau_k
  
  #this could be added to decomp functions
  desorb = .super$pars$mimics[['desorb_p1']] * exp(.super$pars$mimics[['desorb_p2']] * .super$env$clay) * 0.1
  
  .super$pars$k[[5]] = desorb
  
  #total decomp fluxes (for pools with more than one decomp output)
  # d1 = .$decomp.d1(t=t,C=y,i=1) + .decomp.d1(t=t,C=y,i=1, cat = 'cat2', cat_pool = 4)
  d1 = .$decomp.d1(t=t,C=y,i=1) + .$decomp.d1(t=t,C=y,i=1, cat_pool = 4)
  # d2 = .$decomp.d2(t=t,C=y,i=2) + .decomp.d2(t=t,C=y,i=2, cat = 'cat2', cat_pool = 4)
  d2 = .$decomp.d2(t=t,C=y,i=2) + .$decomp.d2(t=t,C=y,i=2, cat_pool = 4)
  # d6 = .$decomp.d6(t=t,C=y,i=6) + .decomp.d6(t=t,C=y,i=6, cat = 'cat2', cat_pool = 4)
  d6 = .$decomp.d6(t=t,C=y,i=6) + .$decomp.d6(t=t,C=y,i=6, cat_pool = 4)
  # d7 = .$decomp.d7(t=t,C=y,i=7) + .decomp.d7(t=t,C=y,i=7, cat = 'cat2', cat_pool = 4)
  d7 = .$decomp.d7(t=t,C=y,i=7) + .$decomp.d7(t=t,C=y,i=7, cat_pool = 4)
  
  #transfers
  t1_to_3 = (.$decomp.d1(t=t,C=y,i=1)/d1)*.super$pars$cue[[1]]
  t1_to_4 = (.$decomp.d1(t=t,C=y,i=1, cat_pool = 4)/d1)*.super$pars$cue2[[1]]
  t2_to_3 = (.$decomp.d2(t=t,C=y,i=2)/d2)*.super$pars$cue[[2]]
  t2_to_4 = (.$decomp.d2(t=t,C=y,i=2, cat_pool = 4)/d2)*.super$pars$cue2[[2]]
  t3_to_5 = .super$pars$mimics[['fSOMp_r_p1']] * exp(.super$pars$mimics[['fSOMp_r_p2']]*.super$env$clay) * 0.1 #0.1 manual calibration from Will's script #fSOMp_r
  t3_to_6 = .super$pars$mimics[['fSOMc_r_p1']] * exp(.super$pars$mimics[['fSOMc_r_p2']]*fmet)*.super$pars$mimics[['fSOMc_r_p3']]  #fSOMc_r
  t3_to_7 = 1 - (t3_to_5 + t3_to_6) #fSOMa_r
  t4_to_5 = .super$pars$mimics[['fSOMp_k_p1']] * exp(.super$pars$mimics[['fSOMp_k_p2']]*.super$env$clay) * 0.1 #fSOMp_k
  t4_to_6 = .super$pars$mimics[['fSOMc_k_p1']] * exp(.super$pars$mimics[['fSOMc_k_p2']]*fmet)*.super$pars$mimics[['fSOMc_k_p3']]  #fSOMc_k
  t4_to_7 = 1 - (t4_to_5 + t4_to_6)  #fSOMa_k
  t7_to_3 = (.$decomp.d7(t=t,C=y,i=7)/d7)*.super$pars$cue[[7]]
  t7_to_4 = (.$decomp.d7(t=t,C=y,i=7, cat_pool = 4)/d7)*.super$pars$cue2[[7]]
  
  #partitioning inputs
  i_LITm = .$input(t) * fmet * (1-.super$pars$mimics[['fi_LITm']])
  i_LITs = .$input(t) * (1-fmet) * (1-.super$pars$mimics[['fi_LITs']])
  i_SOMp = .$input(t) * fmet * .super$pars$mimics[['fi_LITm']]
  i_SOMc = .$input(t) * (1-fmet) * .super$pars$mimics[['fi_LITs']]
  
  #ODE system
  dLITm = i_LITm - d1*.$tcor(.)
  dLITs = i_LITs - d2*.$tcor(.)
  dMICr = d1*.$tcor(.) *t1_to_3 + d2*.$tcor(.) *t2_to_3 + d7*.$tcor(.) *t7_to_3 - .$decomp.d3(t=t,C=y,i=3)
  dMICk = d1*.$tcor(.) *t1_to_4 + d2*.$tcor(.) *t2_to_4 + d7*.$tcor(.) *t7_to_4 - .$decomp.d4(t=t,C=y,i=4)
  dSOMp = i_SOMp + .$decomp.d3(t=t,C=y,i=3)*t3_to_5 + .$decomp.d4(t=t,C=y,i=4)*t4_to_5 - .$desorp.ds5(t=t, C=y,i=5)
  dSOMc = i_SOMc + .$decomp.d3(t=t,C=y,i=3)*t3_to_6 + .$decomp.d4(t=t,C=y,i=4)*t4_to_6 - d6*.$tcor(.)
  dSOMa = .$decomp.d3(t=t,C=y,i=3)*t3_to_7 + .$decomp.d4(t=t,C=y,i=4)*t4_to_7 + .$desorp.ds5(t=t,C=y,i=5) + d6*.$tcor(.) - d7*.$tcor(.)
  list(c(dLITm, dLITs, dMICr, dMICk, dSOMp, dSOMc, dSOMa))
  #print(fmet)
}

f_solver_func_mimics_sat <- function(., t, y, parms){
  
  .super$pars$vmax[[6]] = .super$pars$vmax[[2]]
  .super$pars$vmax2[[6]] = .super$pars$vmax2[[2]]
  
  #dynamic parameters
  fmet = .super$pars$mimics[['fmet_p1']] * (.super$pars$mimics[['fmet_p2']] - .super$pars$mimics[['fmet_p3']]*(.super$env$lignin/.super$env$N))
  
  
  # ensures that tau_mod1 is between two values
  # in Will's script, anpp is multipled by 0 in the manipulation scirpt... which would imply that 
  # tau_mod1 might always be set to 0.6 in the simulations in Ben's paper
  tau_mod1 = min(max(sqrt(.super$env$anpp/.super$pars$mimics[['tau_mod1_p1']]),.super$pars$mimics[['tau_mod1_p2']]),.super$pars$mimics[['tau_mod1_p3']]) 
  tau_r = .super$pars$mimics[['tau_r_p1']] * exp(.super$pars$mimics[['tau_r_p2']] * fmet) *tau_mod1 * .super$pars$mimics[['tau_mod2']]
  tau_k = .super$pars$mimics[['tau_k_p1']] * exp(.super$pars$mimics[['tau_k_p2']] * fmet) *tau_mod1 * .super$pars$mimics[['tau_mod2']]
  
  .super$pars$k[[3]] = tau_r
  .super$pars$k[[4]] = tau_k
  
  #this could be added to decomp functions
  desorb = .super$pars$mimics[['desorb_p1']] * exp(.super$pars$mimics[['desorb_p2']] * .super$env$clay) * 0.1
  
  .super$pars$k[[5]] = desorb
  
  #total decomp fluxes (for pools with more than one decomp output)
  # d1 = .$decomp.d1(t=t,C=y,i=1) + .decomp.d1(t=t,C=y,i=1, cat = 'cat2', cat_pool = 4)
  d1 = .$decomp.d1(t=t,C=y,i=1) + .$decomp.d1(t=t,C=y,i=1, cat_pool = 4)
  # d2 = .$decomp.d2(t=t,C=y,i=2) + .decomp.d2(t=t,C=y,i=2, cat = 'cat2', cat_pool = 4)
  d2 = .$decomp.d2(t=t,C=y,i=2) + .$decomp.d2(t=t,C=y,i=2, cat_pool = 4)
  # d6 = .$decomp.d6(t=t,C=y,i=6) + .decomp.d6(t=t,C=y,i=6, cat = 'cat2', cat_pool = 4)
  d6 = .$decomp.d6(t=t,C=y,i=6) + .$decomp.d6(t=t,C=y,i=6, cat_pool = 4)
  # d7 = .$decomp.d7(t=t,C=y,i=7) + .decomp.d7(t=t,C=y,i=7, cat = 'cat2', cat_pool = 4)
  d7 = .$decomp.d7(t=t,C=y,i=7) + .$decomp.d7(t=t,C=y,i=7, cat_pool = 4)
  
  #transfers
  t1_to_3 = (.$decomp.d1(t=t,C=y,i=1)/d1)*.super$pars$cue[[1]]
  t1_to_4 = (.$decomp.d1(t=t,C=y,i=1, cat_pool = 4)/d1)*.super$pars$cue2[[1]]
  t2_to_3 = (.$decomp.d2(t=t,C=y,i=2)/d2)*.super$pars$cue[[2]]
  t2_to_4 = (.$decomp.d2(t=t,C=y,i=2, cat_pool = 4)/d2)*.super$pars$cue2[[2]]
  t3_to_5 = .super$pars$mimics[['fSOMp_r_p1']] * exp(.super$pars$mimics[['fSOMp_r_p2']]*.super$env$clay) * 0.1 #0.1 manual calibration from Will's script #fSOMp_r
  t3_to_6 = .super$pars$mimics[['fSOMc_r_p1']] * exp(.super$pars$mimics[['fSOMc_r_p2']]*fmet)*.super$pars$mimics[['fSOMc_r_p3']]  #fSOMc_r
  t3_to_7 = 1 - (t3_to_5 + t3_to_6) #fSOMa_r
  t4_to_5 = .super$pars$mimics[['fSOMp_k_p1']] * exp(.super$pars$mimics[['fSOMp_k_p2']]*.super$env$clay) * 0.1 #fSOMp_k
  t4_to_6 = .super$pars$mimics[['fSOMc_k_p1']] * exp(.super$pars$mimics[['fSOMc_k_p2']]*fmet)*.super$pars$mimics[['fSOMc_k_p3']]  #fSOMc_k
  t4_to_7 = 1 - (t4_to_5 + t4_to_6)  #fSOMa_k
  t7_to_3 = (.$decomp.d7(t=t,C=y,i=7)/d7)*.super$pars$cue[[7]]
  t7_to_4 = (.$decomp.d7(t=t,C=y,i=7, cat_pool = 4)/d7)*.super$pars$cue2[[7]]
  
  #partitioning inputs
  i_LITm = .$input(t) * fmet * (1-.super$pars$mimics[['fi_LITm']])
  i_LITs = .$input(t) * (1-fmet) * (1-.super$pars$mimics[['fi_LITs']])
  i_SOMp = .$input(t) * fmet * .super$pars$mimics[['fi_LITm']]
  i_SOMc = .$input(t) * (1-fmet) * .super$pars$mimics[['fi_LITs']]
  
  sat_frac = (1-.super$state$cpools[[5]]/.super$pars$poolmax[[5]])
  SOMp_input = (i_SOMp + .$decomp.d3(t=t,C=y,i=3)*t3_to_5 + .$decomp.d4(t=t,C=y,i=4)*t4_to_5)
  
  #ODE system
  dLITm = i_LITm - d1*.$tcor(.)
  dLITs = i_LITs - d2*.$tcor(.)
  dMICr = d1*.$tcor(.) *t1_to_3 + d2*.$tcor(.) *t2_to_3 + d7*.$tcor(.) *t7_to_3 - .$decomp.d3(t=t,C=y,i=3)
  dMICk = d1*.$tcor(.) *t1_to_4 + d2*.$tcor(.) *t2_to_4 + d7*.$tcor(.) *t7_to_4 - .$decomp.d4(t=t,C=y,i=4)
  dSOMp = SOMp_input*sat_frac - .$desorp.ds5(t=t, C=y,i=5) #this now saturates
  dSOMc = i_SOMc + .$decomp.d3(t=t,C=y,i=3)*t3_to_6 + .$decomp.d4(t=t,C=y,i=4)*t4_to_6 - d6*.$tcor(.) + SOMp_input*(1-sat_frac) #routs somp excess to somc (maybe this should go to soma or co2?)
  dSOMa = .$decomp.d3(t=t,C=y,i=3)*t3_to_7 + .$decomp.d4(t=t,C=y,i=4)*t4_to_7 + .$desorp.ds5(t=t,C=y,i=5) + d6*.$tcor(.) - d7*.$tcor(.)
  list(c(dLITm, dLITs, dMICr, dMICk, dSOMp, dSOMc, dSOMa))
  #print(fmet)
}

############## MILLENNIALv1 #################
# lsoda style function to solve MILLENNIALv1
#############################################
f_solver_func_millennial <- function(., t, y, parms) {
  #State parameters (i.e. calculated parameters)
  
  ##wcor FUNCTION##
  #Sw = (1 / (1+w1*exp(-w2*swc/.35)))
  
  ##tcor FUNCTION##
  #St = (t2 + (t3/pi)* atan(pi*t4*(Temp - t1))) / (t2 + (t3/pi)* atan(pi*t4*(Tref - t1)))
  
  #temperature dependence of CUE
  #ae = assimilation_efficiency; could rename this for consistency but assimilation efficiency is more correct
  #because maintenance respiration is modeled separately
  #function of temp
  #ae = (CUEref - (CUEt*(Temp-Taeref)))
  ae = (.super$pars$cue[[4]] - (.super$pars$millennial[['cuet']]*(.super$env$temp-.super$pars$millennial[['Taeref']])))
  
  #Maximum sorption capacity from Mayes et al. 2012
  #function of clay
  ##state par or calculate in sorption function##
  #Qmax = (BD * 10^(c1*log10(pclay) + c2))/1000 
  
  #Desorption function from Mayes et al. 2012
  #function of pH
  ##state par or calculate in sorption function##
  #Kdm = 10^(-.186*pH - .216)
  
  ######
  #Pools
  ######
  # C1 = POM
  # C2 = MB
  # C3 = MAOM
  # C4 = DOC
  # C5 = Aggregate
  
  #######
  #Fluxes
  #######
  
  # Fa = kb * A*Sw*St
  #Aggregate breakdown
  Fa = .$decomp.d5(t=t,C=y,i=5) * .$tcor(.) *.$wcor(.)
  #decomp = lin
  
  # Fpa = ( Vpa * P) / (Kpa + P)*(1 - A / Amax)*Sw*St
  #  # "The formation of aggregate C (A) from POM follows Michaelis–Menten dynamics,
  # where Vpa is the maximum rate of aggregate formation, Kpa is the half-
  # saturation constant of aggregate formation, and Amax is the maximum capacity
  # of C in soil aggregates.
  # NOT A FUNCTION OF MICROBIAL BIOMASS
  Fpa = .$aggform.a1(t=t,C=y,i=1) * .$tcor(.) *.$wcor(.)
  #aggform = ( Vpa * P) / (Kpa + P)*(1 - A / Amax) where...
  #Amax = 500 in the paper or is caculated based on a separate equation in the source code
  
  # Fpd = Vpd * (P/(Kpd + P)) * (B/(Kpe+B))*Sw*St
  #decomposition of POM governed by double Michaelis-Menten equation
  Fpd = .$decomp.d1(t=t,C=y,i=1) * .$tcor(.) *.$wcor(.)
  #decomp = Vpd * (P/(Kpd + P)) * (B/(Kpe+B))
  
  # Fd = kd * D*St*Sw
  #DOC loss via leaching #think we are missing a k for this...
  Fd = .$decomp.d4(t=t,C=y,i=4) * .$tcor(.) *.$wcor(.)
  #decomp = lin; this is technically leaching and not decomposition though...
  
  # Fdm = Sw*St*D * ((Kdm*Qmax*D)/(1+Kdm*D) - M) / Qmax
  #Sorption
  #This is a two-way equation reflecting the instantaneous net balance between sorption and desorption
  Fdm = .$sorp.s4(t=t,C=y,i=4) * .$tcor(.) *.$wcor(.)
  #sorp = D * ((Kdm*Qmax*D)/(1+Kdm*D) - M) / Qmax
  #Kdm and Qmax state parameters
  #Kdm may need to be calculated dynamically in function (based on pH which could be f(time))
  
  # Fdb = Vdm * D * B/(B+Kdb)*Sw*St
  #Uptake: Microbial uptake of DOC
  Fdb = .$docuptake(t=t,C=y,i=4)* .$tcor(.) *.$wcor(.)
  #uptake = Vdm * D * B/(B+Kdb) #reverse michaelis-menten
  
  # Fma = (Vma * M) / (Kma + M) * (1-A/Amax)*Sw*St
  #Aggregate formation from MAOM
  Fma = .$aggform.a3(t=t,C=y,i=3)* .$tcor(.) *.$wcor(.)
  #aggform = same as for Fpa
  
  # Fbm = kmm * B*Sw*St
  #turnover (and sorption) of microbial necromass
  #here it is more of a turnover process (i.e. independent of mineral properties)
  #but in the source code it depends on saturation level of MAOM (which is really weird because it implies
  #that mineral saturation limits microbial turnover)
  #Though I guess maintenance respiration would increase to compensate??? This would be worth investigating
  Fbm = .$sorp.s2(t=t,C=y,i=2, k_from_list = FALSE, k = .super$pars$millennial[['kmm']])* .$tcor(.) *.$wcor(.) 
  #sorp = lin
  
  # Fmr = km * B*Sw*St
  # Maintenance respiration
  Fmr = .$decomp.d2(t=t,C=y,i=2)* .$tcor(.) *.$wcor(.)
  #decomp = lin
  
  #check inputs
  #print(.$input(t)[[1]])
  #print(.$input(t)[[4]])
  
  #ODE system
  # dP <- pri*i + pa* Fa - Fpa  - Fpd
  dP <- .$input(t)[[1]] + .super$pars$millennial[['pa']] * Fa - Fpa - Fpd
  # dB <- ae * Fdb - Fbm - Fmr
  dB <- ae * Fdb - Fbm - Fmr
  # dM <- Fdm + Fbm - Fma + (1-pa)*Fa  
  dM <- Fdm + Fbm - Fma + (1-.super$pars$millennial[['pa']])*Fa
  # dD <- i * (1-pri) - Fd + Fpd - Fdm - Fdb
  dD <- .$input(t)[[4]] - Fd + Fpd - Fdm - Fdb
  # dA <- Fma + Fpa - Fa
  dA <- Fma + Fpa - Fa
  #print(.$aggform.a1(t=t,C=y,i=1))
  list(c(dP, dB, dM, dD, dA))
}


############## MEND2013 #################
# lsoda style function to solve MEND2013
#########################################
f_solver_func_mend2013 <- function(., t, y, parms) {
  # MEND <- function(t,state,params){
  #   with(as.list(c(state,params)),
  #        { 
  #Input to POC
  #          IP = I*(1-fid)
  IP = .$input(t)[[1]] #*(1-fid) done within input function
  #Input to DOC
  #          ID = I*fid
  ID = .$input(t)[[5]] #*fid done within input function
  #Decomp of POC
  #          F1 = Vp*EP*P/(Kp+P)
  F1 = .$decomp.d1(t=t,C=y,i=1, cat=6) #!!MM decomp
  #Decomp of MAOC
  #          F3 = Vm*EM*M/(Km+M)
  F3 = .$decomp.d2(t=t,C=y,i=2, cat=7) #!! MM decomp 
  #Sorption to Q pool
  #          F4 = Kads*(1-Q/Qmax)*D
  F4 = .$sorp.s5(t=t,C=y,i=5, sat_pool = 3) #!!add function to functions.R # sat function
  #Depsorption from Q pool
  #          F5 = Kdes*(Q/Qmax)
  F5 = .$desorp.ds3(t=t,C=y,i=3) #!!add function to functions.R
  #DOC uptake by microbes
  #          F6 = 1/Ec*(Vd+Mr)*D*B/(Kd+D)
  F6 = 1/.super$pars$cue[[5]]*.$docuptake(t=t,C=y,i=5) #!!add function to functions.R
  #Microbial growth respiration
  #          F9 = (1/Ec-1)*Vd*B*D/(Kd+D)
  F9 = (1/.super$pars$cue[[5]]-1)*.$growthresp(t=t,C=y,i=4) #!!add function to functions.R
  #Microbial Maintenance respiration
  #          F10 = (1/Ec-1)*Mr*B*D/(Kd+D)
  F10 = (1/.super$pars$cue[[5]]-1)*.$maintresp(t=t,C=y,i=4) #!!add function to functions.R
  #Microbial mortality
  #          F12 = (1-Pep-Pem)*Mr*B
  F12 = (1-.super$pars$mend[['Pep']]-.super$pars$mend[['Pem']])*.$decomp.d4(t=t,C=y,i=4) #!!linear decomp
  #Enzyme production
  #          F13EP = Pep*Mr*B
  F13EP = .super$pars$mend[['Pep']]*.$decomp.d4(t=t,C=y,i=4)
  #          F13EM = Pem*Mr*B
  F13EM = .super$pars$mend[['Pem']]*.$decomp.d4(t=t,C=y,i=4)
  
  #Enzyme turnover
  #          F14EP = Rep*EP
  F14EP = .$decomp.d6(t=t,C=y,i=6)
  #          F14EM = Rem*EM
  F14EM = .$decomp.d7(t=t,C=y,i=7)
  #          
  #          dP <- IP + (1-Gd) * F12 - F1
  dP <- IP + (1-.super$pars$mend[['Gd']]) * F12 - F1
  #          dM <- (1-Fd) * F1 - F3
  dM <- (1-.super$pars$mend[['Fd']]) * F1 - F3
  #          dQ <- F4 - F5
  dQ <- F4 - F5
  #          dB <- F6 - (F9 + F10) - F12 - (F13EP + F13EM)
  dB <- F6 - (F9 + F10) - F12 - (F13EP + F13EM)
  #          dD <- ID + Fd * F1 + Gd * F12 + F3 + (F14EP + F14EM) - F6 - (F4 - F5)
  dD <- ID + .super$pars$mend[['Fd']] * F1 + .super$pars$mend[['Gd']] * F12 + F3 + (F14EP + F14EM) - F6 - (F4 - F5)
  #          dEP <- F13EP - F14EP
  dEP <- F13EP - F14EP
  #          dEM <- F13EM - F14EM
  dEM <- F13EM - F14EM
  #          list(c(dP, dM, dQ, dB, dD, dEP, dEM))
  list(c(dP, dM, dQ, dB, dD, dEP, dEM))
  #print(F4)
  #        }
  #   )
  # }
}

#version of MEND where MAOM saturates (ideally would be able to do this with just a different process function in the previous solver function)
f_solver_func_mend2013_sat <- function(., t, y, parms) {
  IP = .$input(t)[[1]] #*(1-fid) done within input function
  ID = .$input(t)[[5]] #*fid done within input function
  F1 = .$decomp.d1(t=t,C=y,i=1, cat=6) #!!MM decomp
  F3 = .$decomp.d2(t=t,C=y,i=2, cat=7) #!! MM decomp 
  F4 = .$sorp.s5(t=t,C=y,i=5, sat_pool = 3) #!!add function to functions.R # sat function
  F5 = .$desorp.ds3(t=t,C=y,i=3) #!!add function to functions.R
  F6 = 1/.super$pars$cue[[5]]*.$docuptake(t=t,C=y,i=5) #!!add function to functions.R
  F9 = (1/.super$pars$cue[[5]]-1)*.$growthresp(t=t,C=y,i=4) #!!add function to functions.R
  F10 = (1/.super$pars$cue[[5]]-1)*.$maintresp(t=t,C=y,i=4) #!!add function to functions.R
  F12 = (1-.super$pars$mend[['Pep']]-.super$pars$mend[['Pem']])*.$decomp.d4(t=t,C=y,i=4) #!!linear decomp
  F13EP = .super$pars$mend[['Pep']]*.$decomp.d4(t=t,C=y,i=4)
  F13EM = .super$pars$mend[['Pem']]*.$decomp.d4(t=t,C=y,i=4)
  F14EP = .$decomp.d6(t=t,C=y,i=6)
  F14EM = .$decomp.d7(t=t,C=y,i=7)
  dP <- IP + (1-.super$pars$mend[['Gd']]) * F12 - F1 +(1-.super$pars$mend[['Fd']]) * F1 * (.super$state$cpools[[2]]/.super$pars$poolmax[[2]]) #returns un-sorbed to unprotected
  dM <- (1-.super$pars$mend[['Fd']]) * F1 * (1-.super$state$cpools[[2]]/.super$pars$poolmax[[2]]) - F3 #transfer into MAOM limited to poolmax par
  dQ <- F4 - F5
  dB <- F6 - (F9 + F10) - F12 - (F13EP + F13EM)
  dD <- ID + .super$pars$mend[['Fd']] * F1 + .super$pars$mend[['Gd']] * F12 + F3 + (F14EP + F14EM) - F6 - (F4 - F5)
  dEP <- F13EP - F14EP
  dEM <- F13EM - F14EM
  list(c(dP, dM, dQ, dB, dD, dEP, dEM))
}

############## MILLENNIALv2 #################
# lsoda style function to solve MILLENNIALv2
#############################################
f_solver_func_millennialV2 <- function(., t, y, parms) {
  #State parameters (i.e. calculated parameters)
  #Equation 4
  ##scalar_wd = (swc / porosity)^0.5
  ##this is now .$wcor(.)
  
  #Equation 15
  #scalar_wb = exp(lambda * -matpot) * (kamin + (1 - kamin) * ((porosity - swc) / porosity)^0.5) * scalar_wd
  ##this is now .$wcor2(.)
  
  #Equation 11
  #Q max; maximum MAOM capacity
  #not sure if this will work?
  .super$pars$poolmax[[3]] = .super$env$BD * .super$env$claysilt * .super$pars$millennialV2[['param_pc']]
  # param_qmax = .super$env$BD * .super$env$claysilt * .super$pars$millennialV2[['param_pc']]
  
  #Equation 10
  #Binding affinity parameter
  kaff_lm = exp(-.super$pars$millennialV2[['sorp_p1']] * .super$env$pH - .super$pars$millennialV2[['sorp_p2']]) * .super$pars$millennialV2[['kld']]
  
  #Equation 3
  #arrhenius modification of vmax for pom decay
  #vmax_pl = .super$pars$millennialV2[['alpha_pl']] * .$tcor.t1(i=1)
  
  #Equation 14
  #arrhenius modification of vmax for mic uptake
  #vmax_lb = .super$pars$millennialV2[['alpha_lb']] * .$tcor.t4(i=4)
  #          
  #          #Fluxes
  
  #define states so that if/then statements work here
  POM = .super$state$cpools[[1]]
  MIC = .super$state$cpools[[2]]
  MAOM = .super$state$cpools[[3]]
  LMWC = .super$state$cpools[[4]]
  AGG = .super$state$cpools[[5]]
  
  #Equation 6
  #agg breakdown
  # AGG -> MAOM + POM
  if(AGG>0){
    f_AG_break = .$decomp.d5(t = t, C=y, i=5) * .$wcor(.)
  }else{
    f_AG_break=0
  }
  
  #Equation 5
  #agg formation
  # POM -> AGG
  # no saturation of agg fraction in v2 I guess
  if(POM>0){
    f_PO_AG = .$aggform.a1(t = t, C=y, i=1) * .$wcor(.)
  }else{
    f_PO_AG=0
  }
  
  #Equation 2
  #RMM decay of POM
  # POM -> LMWC
  if(POM>0 && MIC>0){
    f_PO_LM = .$tcor.t1(i=1) * .$decomp.d1(t=t,C=y,i=1,cat=2) *.$wcor(.)
    #f_PO_LM = vmax_pl * .$wcor(.) * POM * MIC / (km1 + MIC)
  }else{
    f_PO_LM=0
  }
  
  #Equation 8
  #leaching
  # LMWC -> out of system leaching
  if(LMWC>0){
    f_LM_leach = .$decomp.d4(t = t, C=y, i=4) * .$wcor(.)
    #f_LM_leach = k4 * .$wcor(.) * LMWC
  }else{
    f_LM_leach=0
  }
  
  #Equation 9
  #sorption
  # LMWC -> MAOM #this is no longer a double-sided equation, seperate funcs for sorp and desorp
  if(LMWC>0 && MAOM>0){
    f_LM_MA = .$wcor(.) * .$sorp.s4(t = t, C=y, i=4, k_from_list = FALSE, k = kaff_lm, sat_pool = 3)
    #f_LM_MA = .$wcor(.) * kaff_lm * LMWC * (1 - MAOM / param_qmax)
  }else{
    f_LM_MA=0
  }
  
  
  #Equation 12
  #desorption
  # MAOM -> LMWC
  if(MAOM>0){
    f_MA_LM = .$desorp.ds3(t = t, C=y, i=3, k = .super$pars$millennialV2[['kld']], cat = 2)
    #f_MA_LM = kld * MAOM / param_qmax
  }else{
    f_MA_LM=0
  }
  
  #Equation 13
  #microbial uptake via michaelis menten equation
  # LMWC -> MIC
  if(LMWC>0 && MIC>0){
    f_LM_MB = .$wcor2(.) * .$tcor.t4(i=4) * .$docuptake(t=t,C=y,i=4,cat=2)
    #f_LM_MB = vmax_lb * .$wcor2(.) * MIC * LMWC / (km4 + LMWC)
  }else{
    f_LM_MB=0
  }
  
  #Equation 18
  #agg formation from maom
  # MAOM -> AGG
  if(MAOM>0){  
    f_MA_AG = .$aggform.a3(t = t, C=y, i=3) * .$wcor(.)
    #f_MA_AG = k3 * .$wcor(.) * MAOM
  }else{
    f_MA_AG=0
  }
  
  #Equation 16
  #microbial turnover
  # MIC -> MAOM/LMWC
  if(MIC>0){
    f_MB_turn = .$decomp.d2(t = t, C=y, i=2)
    #f_MB_turn = k2 * MIC^2.0
  }else{
    f_MB_turn=0
  }
  
  #Equation 22
  # microbial growth flux, but is not used in mass balance
  
  #Equation 21
  #mic respiration
  # MIC -> atmosphere
  if(MIC>0 && LMWC>0){ 
    f_MB_atm = f_LM_MB * (1 - (.super$pars$cue[[4]] - .super$pars$millennialV2[['cue_t']] * (.super$env$temp - .super$pars$millennialV2[['Taeref']]) ) )
  }else{
    f_MB_atm=0
  }
  #          
  #          #ODE system
  
  #Equation 1
  dPOM = .$input(t)[[1]] + f_AG_break * .super$pars$millennialV2[['pa']] - f_PO_AG - f_PO_LM
  #dPOM = forc_npp * pri + f_AG_break * pa - f_PO_AG - f_PO_LM
  
  #Equation 20
  dMIC = f_LM_MB - f_MB_turn - f_MB_atm
  #dMIC = f_LM_MB - f_MB_turn - f_MB_atm
  
  #Equation 19
  dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * .super$pars$millennialV2[['param_pb']] - f_MA_AG + f_AG_break * (1. - .super$pars$millennialV2[['pa']])
  #dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * param_pb - f_MA_AG + f_AG_break * (1. - pa)
  
  #Equation 7
  dLMWC = .$input(t)[[4]] - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - .super$pars$millennialV2[['param_pb']]) + f_MA_LM
  #dLMWC = forc_npp * (1. - pri) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - param_pb) + f_MA_LM
  
  #Equation 17
  dAGG = f_MA_AG + f_PO_AG - f_AG_break
  #dAGG = f_MA_AG + f_PO_AG - f_AG_break
  
  #print(.$sorp.s4(t = t, C=y, i=4, k_from_list = FALSE, k = kaff_lm, sat_pool = 3))
  #print(f_MA_LM)
  
  list(c(dPOM, dMIC, dMAOM, dLMWC, dAGG))
}

############## CENTURY #################
# lsoda style function to solve CENTURY
#######################################
# - parms is a dummy argument to work with lsoda
# Equations from Abramoff et al. 2021 Millennial v2 paper
f_solver_func_century <- function(., t, y, parms) {
  #Equation B1 
  #just use .$tcor(.) instead of t_scalar
  # t_scalar <- (t2 + (t3 / pi) * atan(pi * t4 * (forc_st - t1))) /
  #   (t2 + (t3 / pi) * atan(pi * t4 *(30.0 - t1)))
  
  #Equation B2
  #just use .$wcor(.) instead of w_scalar (porosity would need to be set to .39 to match Abramoff et al. 2021)
  # w_scalar <- 1.0 / (1.0 + w1 * exp(-w2 * forc_sw/0.39))
  
  #Equation B3
  f_TEX = .super$pars$century[['c1']] - .super$pars$century[['c2']]*.super$env$claysilt*.01
  
  #Equation B4
  # f_StrLitter = StrLitter * k_strlitter * t_scalar * w_scalar * exp(-3*LigFrac)
  f_StrLitter = .$decomp.d1(t = t, C=y, i=1) * .$tcor(.) * .$wcor(.) * exp(-3*.super$env$lignin)
  
  #Equation B5
  # f_MetLitter = MetLitter * k_metlitter * t_scalar * w_scalar  
  f_MetLitter = .$decomp.d2(t = t, C=y, i=2) * .$tcor(.) * .$wcor(.)
  
  #Equation B6
  # f_ACTIVE <- ACTIVE * k_active * t_scalar * w_scalar * f_TEX
  f_ACTIVE = .$decomp.d3(t = t, C=y, i=3) * .$tcor(.) * .$wcor(.) *f_TEX
  
  #Equation B7 
  # f_SLOW <- SLOW * k_slow * t_scalar * w_scalar
  f_SLOW = .$decomp.d4(t = t, C=y, i=4) * .$tcor(.) * .$wcor(.)
  
  #Equation B8
  # f_PASSIVE <- PASSIVE * k_passive * t_scalar * w_scalar
  f_PASSIVE = .$decomp.d5(t = t, C=y, i=5) * .$tcor(.) * .$wcor(.)
  
  #Equation B9
  dStrLitter = .super$pars$input_coefs[[1]] * .super$env$litter - f_StrLitter
  
  #Equation B10
  dMetLitter = .super$pars$input_coefs[[2]] * .super$env$litter - f_MetLitter
  
  #Equation B11
  dACTIVE <- (1-.super$env$lignin) * .super$pars$century[['strlitter_to_active']] * f_StrLitter + .super$pars$century[['metlitter_to_active']] * f_MetLitter  + f_SLOW * .super$pars$century[['slow_to_active']] + f_PASSIVE * .super$pars$century[['passive_to_active']] - f_ACTIVE
  
  #Equation B12
  dSLOW <-  .super$env$lignin * .super$pars$century[['strlitter_to_slow']] * f_StrLitter + f_ACTIVE * (1-f_TEX-.super$pars$century[['active_to_passive']]) - f_SLOW
  
  #Equation B13
  dPASSIVE <- f_ACTIVE * .super$pars$century[['active_to_passive']] + f_SLOW * .super$pars$century[['slow_to_passive']] - f_PASSIVE
  
  list(c(dStrLitter, dMetLitter, dACTIVE, dSLOW, dPASSIVE))
}

############## MEND2019 #################
# lsoda style function to solve MEND
#######################################
# - parms is a dummy argument to work with lsoda
# Equations from Wang et al. 2019 MEND soil moisture paper
f_solver_func_mend2019 <- function(., t, y, parms) {#MEND based on Wang et al. 2019
  
  # #ENV functions
  ## Convert soil water content to soil water potential according to vanGenuchten equation
  ##store this in soil_decomp_water_functions.R
  # fSWC2SWP <- function(SWC0,SWCres,SWCsat,alpha,n,SWPmin){
  #   const_cm2MPa = 98e-6
  #   rlim = 1.01
  #   m = 1-1/n
  #   if(SWC0<=SWCres*rlim) {
  #     SWC = SWCres*rlim
  #   } else {
  #     SWC = SWC0
  #   }
  #   if(SWC < SWCsat){
  #     eff_sat = (SWC - SWCres)/(SWCsat - SWCres)
  #     fSWC2SWP = (1/(eff_sat^(1/m)) - 1)^(1/n)/alpha
  #     fSWC2SWP = -1*fSWC2SWP*const_cm2MPa
  #   } else {
  #     fSWC2SWP = 0
  #   }
  #   fSWC2SWP
  # }
  SWP = water_unit_converter(.)
  # #pH 
  # #below parameters found in MOD_MEND.F90 (MEND-CN model)
  # #cellulases: pHopt = 5.3, pHsens = 1.7
  # #ligninases: pHopt = 4.2, pHsens = 1.4
  # #Mineral-associated org matter: pHopt = 4.8, pHsen = 1.6
  # #pH default value = 6.0 (mean pH of 763 soil samples)
  #should move this eventually...
  fpH <- function(enz, a = pH){
    if(enz == 'cel'){
      pHopt = 5.3
      pHsens = 1.7
    } else if(enz == 'lig'){
      pHopt = 4.2
      pHsens = 1.4
    } else if(enz == 'mom'){
      pHopt = 4.8
      pHsens = 1.6
    }
    exp(-((a - pHopt)/pHsens)^2)
  }
  # 
  # #Arrhenius temp sensitivity
  #adding to soil_decomp_temperature_functions.R as f_tcor_arrhenius_mend
  # fTemp <- function(Ea., Temp. = Temp, Tempref. = Tempref){
  #   # TKref = Tref + const_tmp_C2K
  #   TKref = Tempref. + 273.15
  #   # TK = T + const_tmp_C2K
  #   TK = Temp. + 273.15
  #   # fTArh0 = dexp(Ea*1.D3/const_R * (1.d0/TKref - 1.d0/TK))
  #   exp(Ea.*1000/8.314 * (1/TKref - 1/TK))
  #   # #exp(-Ea.*1000 / 8.314472 * (1 / (Temp.+273) - 1 / (Tempref.+ 273)))
  # }
  # 
  # #calculated par for fSWP
  # SWPD2A = tau * SWPA2D
  SWPD2A = .super$pars$tau * .super$pars$SWPA2D
  # 
  # #water scalars (4 different functions depending on process being modified)
  # fSWP <- function(proc, SWP. = SWP, SWP_FC. = SWP_FC, SWP_min. = SWP_min, b. = b, 
  #                  SWPA2D. = SWPA2D, SWPD2A. = SWPD2A, omega. = omega){
  #   #proc = the process that is being modified (lig, cel, dormancy, resuscitation)
  #   if(proc == 'lig'){ #modifies decomp by oxidative enzymes (specifically, ligninases)
  #     if(SWP. <= -10^2.5) 0
  #     else if(SWP. <= -10^-1.5) 0.625 - 0.25 * log10(-SWP.)
  #     else if(SWP. <= -10^-2.5) 1
  #     else if(SWP. <= -10^-4) (2.5 + 0.4 * log10(-SWP.))/1.5
  #     else 0.6
  #   } else if(proc == 'cel'){ #modifies decomp by hydrolytic enzymes (specifically, cellulases)
  #     if (SWP. <= SWP_min.) 0
  #     else if(SWP. <= SWP_FC.) 1 - (log(SWP./SWP_FC.)/log(SWP_min./SWP_FC.))^b.
  #     else 1
  #   } else if(proc == 'dormancy') abs(-SWP.)^omega. / (abs(-SWP.)^omega. + abs(-SWPA2D.)^omega.)
  #   else if(proc == 'resuscitation') abs(-SWPD2A.)^omega. / (abs(-SWP.)^omega. + abs(-SWPD2A.)^omega.)
  # }
  # 
  # #CUE temp sensitivity
  # Yg = Ygref - kYg * (Temp - Tempref)
  # 
  # #calculated parms 
  # kads = kdes * Kba
  # Vm = -alpha * Vg / (alpha - 1)
  # pEM = pEP * fpEM
  # 
  # #Fluxes
  # #Michaelis-Menten decay of POM1
  # #temp sensitivity of half-saturation constant
  # KP1m = KP1 * fTemp(Ea = Ea[['KP1']])
  # F1 = VdP1 * EP1 * P1 / (KP1m + P1) * fTemp(Ea = Ea[['VdP1']]) * fpH(enz = 'lig') * fSWP(proc = 'lig')
  # F1 = min(F1,P1)
  # 
  # #Michaelis-Menten decay of POM2
  # #temp sensitivity of half-saturation constant
  # KP2m = KP2 * fTemp(Ea = Ea[['KP2']])
  # F2 = VdP2 * EP2 * P2 / (KP2m + P2) * fTemp(Ea = Ea[['VdP2']]) * fpH(enz = 'cel') * fSWP(proc = 'cel')
  # F2 = min(F2,P2)
  # 
  # #Michaelis-Menten decay of MAOM
  # #temp sensitivity of half-saturation constant
  # KMm = KM * fTemp(Ea = Ea[['KM']]) 
  # F3 = VdM * EM * M / (KMm + M) * fTemp(Ea = Ea[['VdM']]) * fpH(enz = 'mom') * fSWP(proc = 'cel')
  # F3 = min(F3,M)
  # 
  # #########
  # ###adding flux early here so so that fraction of DOC can be reserved for microbial uptake as it is in main MEND Fortran code##
  # #DOC uptake by microbes (This is basically Michaelis-Menten where Vmax is (1/Yg*(Vg+Vm)))
  # #temp modification of Vg, Vm, and KD
  # Vgm = Vg * fTemp(Ea = Ea[['Vg']])
  # Vmm = Vm * fTemp(Ea = Ea[['Vm']])
  # KDm = KD * fTemp(Ea = Ea[['KD']])
  # F6 = (1 / Yg) * (Vgm + Vmm) * D * BA / (KDm + D)
  # F6 = min(F6,D)
  # # D = D-F6
  # 
  # # print(paste0(c('mr = ', gamma*Vmm)))
  # #########
  # 
  # #Saturating adsorption of DOC to Q pool
  # # F4 = kads * (1 - Q / Qmax) * (D) * fTemp(Ea = Ea[['kads']])
  # F4 = kads * (1 - Q / Qmax) * (D-F6) * fTemp(Ea = Ea[['kads']]) ##computes based on DOC after microbial uptake (i.e. D-F6 instead of D) as in main code
  # # print(paste(c("Kads = ", kads*fTemp(Ea = Ea[['kads']]))))
  # 
  # #Desorption of Q to DOC pool (increases as Q gets to saturation point)
  # F5 = kdes * (Q / Qmax) * fTemp(Ea = Ea[['kdes']])
  # # print(paste(c("kdes = ", kdes*fTemp(Ea = Ea[['kdes']]))))
  # 
  # #sorption-desorption if-then from FORTRAN script
  # #           if (sOUT % des > (sINP % adsorbent + sOUT % ads)) then
  # if(F5 > (Q+F4)){
  #   F5 = Q+F4
  # } else if(F4>(D+F5)) {
  #   F4 = D+F5
  # }
  # #          sOUT % des = sINP % adsorbent + sOUT % ads
  # #          elseif(sOUT % ads > (sINP % adsorbate + sOUT % des)) then
  # # sOUT % ads = sINP % adsorbate + sOUT % des
  # 
  # # #DOC uptake by microbes (This is basically Michaelis-Menten where Vmax is (1/Yg*(Vg+Vm)))
  # # #temp modification of Vg, Vm, and KD
  # # Vgm = Vg * fTemp(Ea = Ea[['Vg']])
  # # Vmm = Vm * fTemp(Ea = Ea[['Vm']])
  # # KDm = KD * fTemp(Ea = Ea[['KD']])
  # # F6 = (1 / Yg) * (Vgm + Vmm) * D * BA / (KDm + D)
  # # F6 = min(F6,D)
  # 
  # #Dormancy
  # F7 = (1 - D / (KDm + D)) * Vmm * BA * fSWP(proc = 'dormancy')
  # F7 = min(F7,BA)
  # 
  # #Reactivation
  # F8 = D /(KDm + D) * Vmm * BD * fSWP(proc = 'resuscitation')
  # F8 = min(F8,BD)
  # 
  # #Growth respiration for BA (active MBC)
  # F9 = (1 / Yg - 1) * Vgm * BA * D / (KDm + D) 
  # F9 = min(F9,BA)
  # 
  # #Maintenance respiration for BA (active MBC)
  # F10 = (1 / Yg - 1) * Vmm * BA * D / (KDm + D) 
  # F10 = min(F10,BA)
  # 
  # #Maintenance respiration for BD (dormant MBC)
  # F11 = Vm * BD * beta #Vm or Vmm??
  # F11 = min(F11,BD)
  # 
  # #Mortality for BA (active MBC)
  # F12 = gamma * Vmm * BA * fSWP(proc = 'dormancy')
  # F12 = min(F12,BA)
  # print(paste0(c("fSWP_dormancy = ",fSWP(proc = 'dormancy'))))
  # 
  # #Synthesis of Enzymes (lin)
  # F13EP1 = P1 / (P1 + P2 ) * pEP * Vmm * BA
  # F13EP2 = P2 / (P1 + P2 ) * pEP * Vmm * BA
  # F13EM = pEM * Vmm * BA
  # 
  # #Turnover of Enzymes
  # F14EP1 = rE * EP1
  # F14EP2 = rE * EP2
  # F14EM = rE * EM
  # 
  # #ODE system
  # dP1  <- IP1 + (1-gD) * F12 - F1
  # dP2  <- IP2 - F2
  # dM   <- (1-fD) * (F1 + F2) - F3
  # dQ   <- F4 - F5
  # dD   <- ID + fD * (F1 + F2) + gD * F12 + F3 + (F14EP1 + F14EP2 + F14EM) - F6 - (F4 - F5)
  # dBA  <- F6 - (F7 - F8) - (F9 + F10) - F12 - (F13EP1 + F13EP2 + F13EM)
  # dBD  <- (F7 - F8) - F11
  # dEP1 <- F13EP1 - F14EP1
  # dEP2 <- F13EP2 - F14EP2
  # dEM  <- F13EM - F14EM

  list(c(dP1, dP2, dM, dQ, dD, dBA, dBD, dEP1, dEP2, dEM))
}

### END ###