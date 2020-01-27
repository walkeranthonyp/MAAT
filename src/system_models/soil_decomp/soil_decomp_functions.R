################################
#
# MAAT soil_decomp process representation functions (PRFs)
# 
# Matt Craig, AWalker October 2019 
# Carlos Sierra, Markus Mueller (SoilR developers) 
#
################################



### FUNCTIONS
################################


# decay functions
###################

# linear decomp, Oleson 1963
#MEND O6,O7
f_decomp_lin        <- function(.,C,t,i) C[i]*.super$pars$k[[i]]

#density-dependent decomp, used for density-dependent turnover as in Georgiou et al. 2017  
f_decomp_dd_georgiou         <- function(.,C,t,i) (C[i]^.super$pars$beta) * .super$pars$k[[i]]      

# non-linear Michaelis-Menten decomp, a function of microbial biomass
f_decomp_MM_microbe <- function(.,C,t,i) (.super$pars$vmax[[i]]*C[2]*C[i]) / (.super$pars$km[[i]]+C[i])

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


# transfer functions
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
