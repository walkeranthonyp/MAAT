################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# Matt Craig, AWalker October 2019 
# Carlos Sierra, Markus Mueller (SoilR developers) 
#
################################

source('soil_decomp_temperature_functions.R')
source('soil_decomp_water_functions.R')

### FUNCTIONS
################################
#scaling functions
#correct soil protection rates
f_scor_sulman <- function(.,C,t,i) (.super$env$clay/.super$pars$clayref)^.super$pars$qslope_mayes

#this function not incorporated into the 'water_functions' script because it is not normalized (e.g. 0-1)
#thus it is kind of integral to the corpse model under the current parameterization and not substitutable
f_wcor_sulman <- function(.,C,t,i){
  theta <- .super$env$vwc/.super$env$porosity
  theta^3 * (1-theta)^2.5
}

f_tcor_wieder <- function(.,C,t,i){
  exp(.super$env$temp * .super$pars$mimics[['V_slope']] + .super$pars$mimics[['V_int']]) * .super$pars$mimics[['aV']]
}

# f_tcor_abramoff <- function(.,C,t,i){
#   t1 = 15.4
#   t2 = 11.75
#   t3 = 29.7
#   t4 = 0.031
#   (t2 + (t3/pi)* atan(pi*t4*(.super$env$temp - t1))) / (t2 + (t3/pi)* atan(pi*t4*(.super$pars$reftemp - t1)))
# }

# f_wcor_abramoff <- function(.,C,t,i){
#   w1 = 30
#   w2 = 9
#   (1 / (1+w1*exp(-w2*.super$env$vwc/.35)))
#   #.35 is whc I think? this should be specified in env. Perhaps porosity as in sulman.
# }

# f_tcor_arrhenius <- function(.,C,t,i) {
#   # returns a scalar to adjust parameters from reference temp (Tr) to current temp (Ts) 
#   # Arrhenius equation
#   
#   # input parameters  
#   # Ea     -- rate of increase to optimum  (J mol-1)
#   # R      -- molar gas constant J mol-1 K-1
#   
#   # Tr     -- reference temperature (oC) 
#   # Trk    -- reference temperature (K) 
#   # Tsk    -- temperature to adjust parameter to (K) 
#   
#   #convert to Kelvin
#   Trk <- .super$pars$reftemp + 273.15
#   Tsk <- .super$env$temp + 273.15
#   
#   exp( .super$pars$ea[[i]]*(Tsk-Trk) / (.super$pars$R*Tsk*Trk) )
# }


# decay functions
###################

# linear decomp, Oleson 1963
#MEND O6,O7
f_decomp_lin        <- function(.,C,t,i, k_from_list = TRUE, k = NULL) {
  if(k_from_list == TRUE){
    C[i]*.super$pars$k[[i]]
  } else {
    C[i]*k
  }
}

#density-dependent decomp, used for density-dependent turnover as in Georgiou et al. 2017  
f_decomp_dd_georgiou         <- function(.,C,t,i) (C[i]^.super$pars$beta) * .super$pars$k[[i]]      

# non-linear Michaelis-Menten decomp, a function of microbial biomass
f_decomp_MM_microbe <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[2]*C[i]) / (.super$pars$km[[i]]+C[i])

#reverse Michaelis-Menten decomp
#C[4] is microbial biomass in CORPSE
f_decomp_rmm <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[4] + .super$pars$km[[i]])

# double Michaelis-Menten decomp
f_decomp_dmm <- function(.,C,t,i){
  .super$pars$vmax[[i]] * (C[1]/(.super$pars$km[[1]] + C[1])) * (C[2]/(.super$pars$rkm[[1]]+C[2]))
}

f_zero <- function(.,C,t,i) 0

#########
#MEND-SPECIFIC FUNCTIONS
#########
# non-linear Michaelis-Menten decomp of POM, a function of POM-specific enzymes
#MEND O1
f_decomp_MM_enzpom <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[6]*C[i]) / (.super$pars$km[[i]]+C[i])

#MEND losses of microbial biomass toward CO2 (Fr) and toward other pools (Fe)
#MEND O2
f_decomp_mbc_mend <- function(.,C,t,i) {
  ( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
    C[2]*.super$pars$mr #Fe
}

# non-linear Michaelis-Menten decomp of MAOM, a function of MAOM-specific enzymes
#MEND O3
f_decomp_MM_enzmaom <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[7]*C[i]) / (.super$pars$km[[i]]+C[7])

#density-dependent decomp, as used in MEND for desorption from the Q pool
#MEND O4
f_decomp_dd_mend <- function(.,C,t,i) .super$pars$k[[i]]*(C[i]/.super$pars$poolmax[[i]])

#MEND losses of dissolved organic matter toward mbc uptake (Fu) and adsorption to mineral surfaces (Fa)
#MEND O5
f_decomp_doc_mend <- function(.,C,t,i) {
  C[5]*((.super$pars$vmax[[i]]+.super$pars$mr)/.super$pars$cue[[i]])*(C[2]/(.super$pars$km[[i]]+C[5])) +  #Fu
  C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]])             #Fa
}

#########
#CORPSE-SPECIFIC FUNCTIONS
#########
f_decomp_rmm_sulman <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[4] + .super$pars$km[[i]]*(C[1]+C[2]+C[3]))

f_micturn_sulman <- function(.,C,t,i) {
  (C[i] - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[i]] 
}

#########
#alternative CORPSE hyps
#########
#michaelis-menten version of decomp function
#would need to change km parameterization
f_decomp_mm_sulman <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[4]*C[i]) / (C[i] + .super$pars$km[[i]]*(C[1]+C[2]+C[3]))

#CORPSE density dependence version of microbial turnover
f_micturn_sulman_dd <- function(.,C,t,i) {
  (C[i]^.super$pars$beta - .super$pars$minmic * (C[1]+C[2]+C[3]))/.super$pars$k[[i]] 
}

#adding a decomp function that saturates and returns excess to unprotected pool
f_decomp_lin_sat_corpse        <- function(.,C,t,i) {
  protected_max = 2 #arbitrary
    C[i]*.super$pars$k[[i]]* (1-(C[5] + C[6] + C[7])/protected_max)
}

# f_desorp_lin_sat_corpse <- function(.,C,t,i) {
#   C[i]*.super$pars$k[[i]]
# }

#########
#MILLENNIAL-SPECIFIC FUNCTIONS
#########
f_aggform_abramoff <- function(.,C,t,i, agg_pool = 5){
  if(i==1){
     # print(C[i])
     # print(.super$pars$millennial[['Vpa']])
     # print(.super$pars$millennial[['Kpa']])
     # print(.super$pars$poolmax[[5]])
    ( .super$pars$millennial[['Vpa']] * C[i]) / (.super$pars$millennial[['Kpa']] + C[i])*(1 - C[agg_pool] / .super$pars$poolmax[[5]])
  }
  else if(i==3){
    ( .super$pars$millennial[['Vma']] * C[i]) / (.super$pars$millennial[['Kma']] + C[i])*(1 - C[agg_pool] / .super$pars$poolmax[[5]])  
  }
}

f_docsorp_abramoff <- function(.,C,t,i){
  Qmax = (.super$env$BD * 10^(.super$pars$millennial[['c1']]*log10(.super$env$clay) + .super$pars$millennial[['c2']]))/1000
  Kdm = 10^(-.186*.super$env$pH - .216)
  C[i] * ((Kdm*Qmax*C[i])/(1+Kdm*C[i]) - C[3]) / Qmax
}

f_docuptake_abramoff <- function(.,C,t,i){
  #this is a reverse michaelis-menten equation so could potentially use a f_decomp_rmm function instead
  .super$pars$millennial[['Vdm']] * C[i] * C[2]/(C[2]+.super$pars$millennial[['Kdb']])
}


#########
#MIMICS-SPECIFIC FUNCTIONS
#########
#environmental controls are embedded,
#would be ideal to pull these out
#for now we could define all the parameters with if...then statements in the solver function
f_decomp_rmm_wieder <- function(.,C,t,i,cat_pool = 3){
  if(cat_pool == 3){
  if(i==7){
    pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
    Km = .super$pars$km[[7]] * pscalar
  } else if(i==6) {
    #km par same as structural litter (km2)
    Km = .super$pars$km[[2]] *.super$pars$mimics[['ko_r']]
  }  else {
    Km = .super$pars$km[[i]]
  }
  #correcting Km for temperature
  Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
  ###RMM equation
  C[i] * .super$pars$vmax[[i]] * C[cat_pool] / (Km_cor + C[cat_pool])
  ###
  } else {
      if(i==7){
        pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
        Km = .super$pars$km2[[7]] * pscalar
      } else if(i==6) {
        Km = .super$pars$km2[[2]] *.super$pars$mimics[['ko_k']]
      }  else {
        Km = .super$pars$km2[[i]]
      }
      Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
      ###RMM equation
      C[i] * .super$pars$vmax2[[i]] * C[cat_pool] / (Km_cor + C[cat_pool])
      ###
  }
}

# f_decomp_rmm_wieder <- function(.,C,t,i,cat = 'cat1', cat_pool = 3){
#   if(i==7){
#     pscalar = .super$pars$mimics[['pscalar_p1']] * exp(.super$pars$mimics[['pscalar_p2']]*sqrt(.super$env$clay))
#     Km = .super$pars$km[[7]][[cat]] * pscalar
#   } else if(i==6) {
#     #km par same as structural litter (km2)
#     Km = .super$pars$km[[2]][[cat]] *.super$pars$mimics[['ko']][[cat]]
#   }  else {
#     Km = .super$pars$km[[i]][[cat]]
#   }
#   #correcting Km for temperature
#   Km_cor = exp(.super$env$temp * .super$pars$mimics[['K_slope']] + .super$pars$mimics[['K_int']]) * .super$pars$mimics[['aK']] /Km
# 
#   C[i] * .super$pars$vmax[[i]][[cat]] * C[cat_pool] / (Km_cor + C[cat_pool])
# }

# transfer functions


#transfer functions
###################

# transfer all or nothing
f_transfer_all  <- function(.,C,t,from,to) 1
f_transfer_zero <- function(.,C,t,from,to) 0

# CUE or carbon transfer efficiency sets transfer from one pool to another
#MEND15
f_transfer_cue  <- function(.,C,t,from,to) .super$pars$cue[[from]]    

# CUE / transfer efficiency sets transfer subject to a maximum pool size 
# - can be used both for saturating MAOM pool and density dependent microbial growth efficiency
f_transfer_cue_sat <- function(.,C,t,from,to) .super$pars$cue[[from]] * (1-C[to]/.super$pars$poolmax[[to]])

#transfers remainder of cue function to another pool instead of CO2
#MEND13
f_transfer_cue_remainder <- function(.,C,t,from,to) (1-.super$pars$cue[[from]])

#transfer from mbc to pom in mend: (1-gd)(1-pe)(fe/f_decomp_mbc_mend)
f_transfer_mend21 <- function(.,C,t,from,to){
  (1-.super$pars$cue[[from]])*
  (1-(.super$pars$pep+.super$pars$pem))*
  ((C[2]*.super$pars$mr)/
     (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
     C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend25 <- function(.,C,t,from,to){
  (.super$pars$cue[[from]])*
    (1-(.super$pars$pep+.super$pars$pem))*
    ((C[2]*.super$pars$mr)/
       (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
       C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend26 <- function(.,C,t,from,to) {
  .super$pars$pep * 
    ((C[2]*.super$pars$mr)/
       (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
       C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend27 <- function(.,C,t,from,to) {
  .super$pars$pem * 
    ((C[2]*.super$pars$mr)/
       (( C[2]*((1/.super$pars$cue[[5]])-1)*((C[5]*(.super$pars$vmax[[5]]+.super$pars$mr))/(.super$pars$km[[5]]+C[5]))) + #Fr
       C[2]*.super$pars$mr)) #Fe)
}

f_transfer_mend52 <- function(.,C,t,from,to){
 # print(C)
  #print(from)
  #print(to)
  #print(t)
  #print(.super$pars$poolmax[[4]])
  #print(C[5])
  #print(C[4])
  (C[5]*((.super$pars$vmax[[from]]+.super$pars$mr)/.super$pars$cue[[from]])*(C[2]/(.super$pars$km[[from]]+C[5]))) / #Fu
  (C[5]*((.super$pars$vmax[[from]]+.super$pars$mr)/.super$pars$cue[[from]])*(C[2]/(.super$pars$km[[from]]+C[5])) +  #Fu
    C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]]))   #Fa
}

f_transfer_mend54 <- function(.,C,t,from,to){
  (C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]]))/ #Fu
    (C[5]*((.super$pars$vmax[[from]]+.super$pars$mr)/.super$pars$cue[[from]])*(C[2]/(.super$pars$km[[from]]+C[5])) +  #Fu
       C[5]*((.super$pars$Kads*(.super$pars$poolmax[[4]]-C[4]))/.super$pars$poolmax[[4]]))   #Fa
}


### END ###
