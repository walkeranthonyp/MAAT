################################
#
# MAAT soil_decomp generalisation functions 
# - based entirely on functions in SoilR package functions modifed to work with proto objects and MAAT    
# 
# Carlos Sierra, Markus Mueller (SoilR developers) 
# AWalker, Matt Craig, October 2019 
#
################################



# MODIFIED SoilR FUNCTIONS
################################

# input matrix, single column, rows = cpools_n
# - this is where inputs would be divided among pools
# - APW I don't think this needs to be part of the solver
f_input <- function(., t ) {
  #why is .$env$litter not .super$env$litter??
  .$env$litter * matrix(unlist(.super$pars$input_coefs)[1:.super$pars$n_pools], ncol=1 )
}

f_input_mimics <- function(., t ) {
  EST_LIT_in = .super$env$anpp / (365*24) # gC/m2/h (from gC/m2/y)
  EST_LIT    = EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h)
  EST_LIT/.super$env$depth 
 
}

## decomp matrix, single column, rows = n_pools
# - APW seems a little convoluted, why can't just use 1:n_pools for the id
f_DotO <- function(., C, t ) {
  # search fns proto object for functions named starting with 'decomp.' and use to make dnames 
  dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  # get integer id's of the decomp function for each pool
  id     <- sub('decomp.d', '', dnames )
  # call functions and create decomp matrix 
  m      <- matrix(ncol=1, nrow=.super$pars$n_pools )
  for(i in id) { 
    #print(i); print(.[[paste0('decomp.d',i)]]) 
    m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t,i=as.numeric(i))
  }
  m
}


## transfer matrix, square, n_pools extent
f_transfermatrix <- function(., C, t ) {
  # search fns proto object for functions named starting with 'transfer.' and use to make tnames 
  tnames <- grep('transfer\\.', names(.), value=T )
  # get integer id's of the transfer functions
  id     <- sub('transfer.t', '', tnames )
  # remove transfers that are not required, i.e. are 0 or from/to pools > n_pools
  id     <- id[!is.na(.super$fnames$transfer[paste0('t',id)])]
  # get integer id's of the from and to pools of the transfers 
  idm    <- apply(as.matrix(id,nrow=1), 1, function(i) as.numeric(unlist(strsplit(i,'_to_'))) )
  idm    <- apply(idm,1,function(v) v[v<=.super$pars$n_pools] )
  # call functions and create transfer matrix 
  m      <- -1 * diag(nrow=.super$pars$n_pools)
  #print(idm)
  #print(dim(idm)[1])
  #print(m)
  for(i in 1:dim(idm)[1]) {
    ss <- idm[i,]
    #print(i); print(ss); print(.[[paste0('transfer.t',ss[1],'_to_',ss[2])]]) 
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',ss[1],'_to_',ss[2])]](C=C,t=t,from=ss[1],to=ss[2])
  }
  m
}
  
  
# deSolve/lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t) 
  list(as.vector(YD))
}

# lsoda style function to solve CORPSE
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
  t3_to_5 = .super$pars$mimics[['fSOMp_r_p1']] * exp(.super$pars$mimics[['fSOMp_r_p2']]*.super$pars$clay) * 0.1 #0.1 manual calibration from Will's script #fSOMp_r
  t3_to_6 = .super$pars$mimics[['fSOMc_r_p1']] * exp(.super$pars$mimics[['fSOMc_r_p2']]*fmet)*.super$pars$mimics[['fSOMc_r_p3']]  #fSOMc_r
  t3_to_7 = 1 - (t3_to_5 + t3_to_6) #fSOMa_r
  t4_to_5 = .super$pars$mimics[['fSOMp_k_p1']] * exp(.super$pars$mimics[['fSOMp_k_p2']]*.super$pars$clay) * 0.1 #fSOMp_k
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
  dSOMc = i_SOMc + .$decomp.d3(t=t,C=y,i=3)*t3_to_6 + .$decomp.d4(t=t,C=y,i=4)*t4_to_5 - d6*.$tcor(.)
  dSOMa = .$decomp.d3(t=t,C=y,i=3)*t3_to_7 + .$decomp.d4(t=t,C=y,i=4)*t4_to_7 + .$desorp.ds5(t=t,C=y,i=5) + d6*.$tcor(.) - d7*.$tcor(.)
  list(c(dLITm, dLITs, dMICr, dMICk, dSOMp, dSOMc, dSOMa))
}

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
  #DOC loss via leaching
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

### END ###
