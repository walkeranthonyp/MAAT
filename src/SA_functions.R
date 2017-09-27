################################
#
# Sensitivity Analysis functions for MAAT model 
# 
# AWalker August 2015
#
################################

#####################################################
calc_process_sensitivity <- function(df,res,col,colname,n,nA=1,nB=3,...) {
  
  # assume equal probability for each model
  mpA   <- rep(1/nA,nA)
  mpB   <- rep(1/nB,nB)
  
  # extract model output of interest (i.e. Delta)
  # f     <- as.vector(df[,which(names(df)==colname)])
  f     <- as.vector(df[[res]][,col])
  
  # total mean & variance in Delta
  # - variance should be normalised by the total sample size, not sample size - 1
  var_t <- var(f) 
  
  # reshape f vector into array dim = nA,n,nB,n - first dimension cycles first, like FORTRAN
  f     <- array(f,dim=c(n,nB,n,nA))
  f     <- aperm(f,4:1)
  
  # call function to calculate variance in Delta caused by variability in process x
  var_A <- func_process_sensitivity(1,f,nA,nB,mpA,mpB,n)
  var_B <- func_process_sensitivity(2,f,nA,nB,mpA,mpB,n)
  var_p <- c(var_A,var_B)
  
  # output a list of a scalar of total variance, and a matrix of variance caused by each process, and their proportion of total variance
  Tvar <- c(mean=mean(f),total_var=var_t,total_sd=var_t^0.5)
  
  list(Tvar=Tvar,par_var=var_p,sensitivity=var_p/var_t)
}



#####################################################
func_process_sensitivity <- function(ps,f,nA,nB,MA,MB,n) {
  # Process level variance decomposition after Ye et al (in review)  
  
  # if statements about which process is being quantified
  if(ps == 1) {
    nLoop1 <- nA
    nLoop2 <- nB  
    mLoop1 <- MA
    mLoop2 <- MB  
  } else {
    nLoop1 <- nB
    nLoop2 <- nA  
    mLoop1 <- MB
    mLoop2 <- MA  
  }
  
  # define the arrays
  E_tLoop2    <- array(0,dim=c(nLoop1,n,nLoop2))
  E_Loop2     <- matrix(0,nLoop1,n)
  E_tLoop1_2  <- matrix(0,nLoop1)
  E_tLoop1    <- matrix(0,nLoop1)
  
  # Loop over the models of the process in question
  for( i in 1:nLoop1 ) {  
    
    # Loop over the parameter realizations of the models of the process in question
    for( j in 1:n ) {      
      
      # Loop over the models of the other processes
      for( k in 1:nLoop2 ) {
        # construct subscripts matrix for output array
        smat            <- cbind(i,j,k,1:n)
        if(ps==2) smat  <- smat[,c(3,4,1,2)]  
        # Calculate the partial mean EA considering uncertain parameters            
        E_tLoop2[i,j,k] <- mean(f[smat])
      }
      
      # Calculate the model averaging of the other process
      E_Loop2[i,j] <- sum(mLoop2*E_tLoop2[i,j,])
    }
    
    # Calculate the partial mean of the process in question (the mean of each representation of the process in question) considering uncertain parameters    
    E_tLoop1_2[i] <- mean(E_Loop2[i,]^2)
    E_tLoop1[i]   <- mean(E_Loop2[i,])
  }
  
  # Model averaging of each representation of the process in question
  E_Loop1_2  <- sum(mLoop1*E_tLoop1_2)
  E_Loop1    <- sum(mLoop1*E_tLoop1)
  
  # Calculate the partial variance VB
  Var_Loop1 <- E_Loop1_2 - (E_Loop1)^2
  return(Var_Loop1)
}



#####################################################
func_sobol_sensitivity <- function(sn,k,ABout,ABiout,pnames=NULL) {
  # Calculates Sobol first order and total sensitivity indices using Saltelli method

  # expects: 
  # - sn     - the number of iterations in the basic mc
  # - k      - the number of parameters that have been varied
  # - ABout  - a vector of output from the 2sn sobol mc, i.e. from param matrices A and B
  # - ABiout - an sn * k matrix of output from the k ABi matrices
  
  # returns a list with objects:
  # - mean         - overall mean of model output
  # - total_var    - overall variance of model output
  # - total_sd     - overall standard deviation of model output
  # - pvar         - vector of length k, partial variance for each parameter
  # - pe           - vector of length k, partial mean for each parameter
  # - sensitivity1 - vector of length k, first order sensitivity (main effect)
  # - sensitivityT - vector of length k, total sensitivity (including all interactions)
  
  # initialise partial variance and mean vectors
  vv <- numeric(k)
  ev <- numeric(k)
  
  for(p in 1:k) {
    # Calculate the partial variance caused by single parameter
    # following the equation (b) in Table 2 Saltelli et al. (2010).
    vv[p] <- sum( ABout[(sn+1):(2*sn)]*(ABiout[1:sn,p] - ABout[1:sn]) ) / sn      

    # Calculate the partial mean for total effect sensitivity index
    # following the equation (f) in Table 2 Saltelli et al. (2010),
    ev[p] <- sum( (ABout[1:sn] - ABiout[1:sn,p])^2 ) / (2*sn)      
  }
  
  # calculate the total variance
  v_t <- var(ABout)
  
  #Calculate first-order sensitivity index
  si  <- vv/v_t
  
  #Calculate total-effect sensitivity index
  st  <- ev/v_t
  
  # name parameter SA vectors
  if(!is.null(pnames)) {
    names(si) <- pnames
    names(st) <- pnames
  }
   
  # a model and scenario name would be useful in this output list - would need to be added to the output from the wrapper 
  Tvar <- c(mean=mean(ABout),total_var=v_t,total_sd=v_t^0.5)

  # output list
  list(Tvar=Tvar,pvar=vv,pe=ev,sensitivity1=si,sensitivityT=st)
#   list(Tvar=Tvar,par_var=rbind(vv,ev),sensitivity=rbind(si,st))
#   list(mean=mean(ABout),total_var=v_t,total_sd=v_t^0.5,par_var=rbind(vv,ev),sensitivity=rbind(si,st))
}



########################################################
calc_parameter_sensitivity <- function(sn,outdata,delta) {
  # calcualtes Sobol sensitivity indices using output from a Saltelli algorithm
  # Expects the use of multiple models and scenarios

  # Expects outdata to be a list of two objects - AB: a 4D array & ABi: a 5D array  
  # AB, list element 1 - a 4D numeric array AB 
  # - dim 1 - models  (number - nmod)
  # - dim 2 - scenarios/environments  (number - nscen)
  # - dim 3 - parameter samples (number - 2*sn)
  # - dim 4 - output columns from the model (variable of interest subscript - delta)
  # ABi, list element 2 - a 5D numeric array ABi 
  # - dim 1 - models  (number - nmod)
  # - dim 2 - scenarios/environments  (number - nscen)
  # - dim 3 - parameter samples (number - sn)
  # - dim 4 - output columns from the model (variable of interest subscript - delta)
  # - dim 5 - parameters (number - k)
  
  # list processing function 
  # - creates a weighted average across list elements from a Sobol list 
  # - same dimensions as the list dimensions of ABi
  average_list <- function(inlist,prob) {
    total_var <- sum( sapply(inlist,function(l) l$Tvar[2]) * prob)
    olist <- 
      list(
        Tvar = c(
          mean        = sum( sapply(inlist,function(l) l$Tvar[1]) * prob),
          total_var   = total_var,
          total_sd    = total_var^0.5
        ),
        pvar          = apply( sapply(inlist,function(l) l$pvar) ,1,function(v) sum(v*prob)),
        pe            = apply( sapply(inlist,function(l) l$pe)   ,1,function(v) sum(v*prob)),
        sensitivity1  = apply( sapply(inlist,function(l) l$pvar) ,1,function(v) sum(v*prob))/total_var,
        sensitivityT  = apply( sapply(inlist,function(l) l$pe)   ,1,function(v) sum(v*prob))/total_var
      )
    
    # name parameter vectors
    if(!is.null(names(inlist[[1]]$sensitivity1))) {
      names(olist$sensitivity1) <- names(inlist[[1]]$sensitivity1)  
      names(olist$sensitivityT) <- names(inlist[[1]]$sensitivity1)  
    }
    
    olist
  }
  
  # calculate number of models and scenarios
  nmod    <- dim(outdata$AB)[1]
  nscen   <- dim(outdata$AB)[2]
  sn      <- dim(outdata$ABi)[3]
  k       <- dim(outdata$ABi)[5]
  
  # assume models are of equal probability
  pmod    <- rep(1/nmod,nmod)
  # assume scenarios are of equal probability
  pscen   <- rep(1/nscen,nscen)
  
  # calculate Sobol for each scenario and model combination
  n <- 1
  # model loop
  for(m in 1:nmod) {
    # scenario loop
    for(e in 1:nscen) {
      out1 <- list(func_sobol_sensitivity(sn, k, outdata$AB[m,e,,delta], outdata$ABi[m,e,,delta,], pnames=dimnames(outdata$ABi)[5] ))      
      out2 <- if(e==1) out1 else c(out2, out1 )
    }
    sensout <- if(m==1) list(out2) else c(sensout, list(out2) )
  }
  
  # calculate Sobol including scenario uncertainty
  # hold model constant  
  for(m in 1:nmod) {
    out1 <- average_list(sensout[[m]], pscen )
    outs <- if(m==1) list(out1) else c(outs,list(out1))
  }
  
  # calculate Sobol including model uncertainty
  # hold scenario constant
  for(e in 1:nscen) {
    out1 <- average_list(lapply(sensout,function(l) l[[e]]) , pmod)
    outm <- if(e==1) list(out1) else c(outm,list(out1))
  }
  
  # calculate Sobol including model & scenario uncertainty
  outsm <- average_list(outm,pscen)
  # could also read (I think):
  # outsm <- average_list(outs,pmod)
  
  # output list
  list(individual=sensout,incscenario=outs,incmodel=outm,incmodelscenario=outsm)
}


