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
f_input <- function(., t ) {
  
  .$env$litter * matrix(unlist(.super$pars$input_coefs)[1:.super$pars$n_pools], ncol=1 )
}


## decomp matrix, single column, rows = n_pools
f_DotO <- function(., C, t ) { 
  dnames <- grep('decomp\\.', names(.), value=T )[1:.super$pars$n_pools]
  id     <- sub('decomp.d', '', dnames )
  m      <- matrix(ncol=1, nrow=.super$pars$n_pools )
  # call functions and create decomp matrix 
  for(i in id) { m[as.numeric(i),] <- .[[paste0('decomp.d',i)]](C=C,t=t,i=as.numeric(i))
    #print(id); print(.[[paste0('decomp.d',i)]]) 
  }
  m
}


## transfer matrix, square, n_pools extent
f_transfermatrix <- function(., C, t ) {
  tnames <- grep('transfer\\.', names(.), value=T )
  id     <- sub('transfer.t', '', tnames )
  # remove transfers that are not required, i.e. are 0 or from/to pools > n_pools
  id     <- id[!is.na(.super$fnames$transfer[paste0('t',id)])]
  idm    <- apply(as.matrix(id,nrow=1), 1, function(i) as.numeric(unlist(strsplit(i,'_to_'))) )
  idm    <- apply(idm,1,function(v) v[v<=.super$pars$n_pools] )
  # call functions and create transfer matrix 
  m      <- -1 * diag(nrow=.super$pars$n_pools)
  #print(idm)
  #print(dim(idm)[1])
  #print(m)
  for(i in 1:dim(idm)[1]) {
    ss <- idm[i,]
    #print(i) 
    #print(ss)
    #print(.[[paste0('transfer.t',ss[1],'_to_',ss[2])]]) 
    m[matrix(rev(ss),nrow=1)] <- .[[paste0('transfer.t',ss[1],'_to_',ss[2])]](C=C,t=t,from=ss[1],to=ss[2])
  }
  m
}
  
  
# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func <- function(., t, y, parms) {
  YD = .$transfermatrix(y,t) %*% .$DotO(y,t) + .$input(t) 
  list(as.vector(YD))
}

# lsoda style function to solve
# - parms is a dummy argument to work with lsoda
f_solver_func_corpse <- function(., t, y, parms) {
  
  # #State parameters (i.e. calculated parameters)
  # Icr = root.litter*frroot + leaf.litter*frleaf      #partitions root inputs into simple and resistant pools
  # Ics = root.litter+leaf.litter - Icr                #partitions leaf inputs into simple and resistant pools
  # 
  # Q = (clay/5)^.4833 #The above equation simplifies to this where 5 is the reference clay and .4833 is the slope from Mayes et al. 2012
  # 
  # vmax = vmaxref/exp(-Ea/(8.314472*293.15)) * exp(-Ea/(8.314472*temp)) #Arrhenius modification of maximum decomp rate, this was not correct in Sulman supplement
  # litter = c(Cs, Cr, Cn) #trying to define litter as a vector to simplify code
  # fracmic = M/sum(litter) #microbial biomass as fraction of unprotected C
  # theta = vwc/porosity
  # protectedC = c(Ps, Pr, Pn)
  # 
  ############################################
  #f_decomp_rmm_sulman  # D = litter * vmaxref * M / (M + kC*sum(litter))
  # 
  ############################################
  #f_micturn_sulman #microbeTurnover = max(0,(M-minMicrobeC*sum(litter))/Tmic)    (make Tmic = 1/Tmic in pars then multiply)
  ############################################
  # 
  ############################################
  #f_decomp_lin #S = litter * protection_rate 
  #(* Q) this is f_scor_sulman
  ############################################
  # 
############################################
  #f_decomp_lin # U = protectedC/tProtected
  #make sure the tProtected parameter is (1/tProtected)
############################################
  
  
    # 
  # #ODE system
  # dCs <- Ics + U[1] - S[1] - D[1]
  dCs <- .$input[[1]](t) + .$desorp[[5]](.) - .$sorp[[1]](.)*.$scor(.) - .$decomp[[1]](.)*.$wcor(.)*.$tcor[[1]](.)
  # dCr <- Icr + U[2] - S[2] - D[2] 
  dCr <- .$input[[2]](t) + .$desorp[[6]](.) - .$sorp[[2]](.)*.$scor(.) - .$decomp[[2]](.)*.$wcor(.)*.$tcor[[2]](.)
  # dCn <- Om*et + U[3] - S[3] - D[3] 
  dCr <- .$decomp[[4]](.)*.super$pars$cue[[4]] + .$desorp[[7]](.) - .$sorp[[3]](.)*.$scor(.) - .$decomp[[3]](.)*.$wcor(.)*.$tcor[[3]](.)
  # dM  <- sum(D*eup) - Om
  dM <- .$decomp[[1]](.)*.super$pars$cue[[1]]+.$decomp[[2]](.)*.super$pars$cue[[2]]+.$decomp[[3]](.)*.super$pars$cue[[3]] - .$decomp[[4]](.)
  # dPs <- S[1] - U[1]
  dPs <- .$sorp[[1]]() - .$desorp[[5]]()
  # dPr <- S[2] - U[2]
  dPr <- .$sorp[[2]]() - .$desorp[[6]]()
  # dPn <- S[3] - U[3]
  dPn <- .$sorp[[3]]() - .$desorp[[7]]()
  # list(c(dCs, dCr, dCn, dM, dPs, dPr, dPn))
  list(c(dCs, dCr, dCn, dM, dPs, dPr, dPn))
}


### END ###
