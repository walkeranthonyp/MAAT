################################
#
# Sensitivity Analysis functions for MAAT model 
# 
# AWalker August 2015
#
################################



# Functions related to the process sensitivity index using Ye method (Dai et al 2017 WRR).
#####################################################


# Calculate the process sensitivity index using Ye method (Dai etal 2017 WRR).
calc_process_sensitivity <- function(delta) {
  
  # delta should be a 4D array comprised of only a single model output variable
  # - dim 1 (rows)      process(es) B parameter sample 
  # - dim 2 (columns)   process(es) B representation(s)
  # - dim 3 (slices)    process A parameter sample
  # - dim 4 (cube rows) process A representation
  
  # total variance in Delta
  # - variance should be normalised by the total sample size, not sample size - 1
  var_t <- var(delta) 
  
  # call function to calculate variance in Delta caused by variability in process x
  pvar  <- pvar_process_sensitivity(delta)

  # output a list of a scalar of total variance, and a matrix of variance caused by each process, and their proportion of total variance
  Tvar  <- c(mean=mean(delta), total_var=var_t, total_sd=var_t^0.5 )
  
  # output list
  list(Tvar=Tvar, par_var=pvar, sensitivity=pvar/var_t )
}


# Calculate partial variance for the process sensitivity index using Ye method (Dai etal 2017 WRR).
pvar_process_sensitivity <- function(delta) {

  # get dim extents, delta should be a 4D array
  # - dim 1 (rows)      process(es) B parameter sample 
  # - dim 2 (columns)   process(es) B representation(s)
  # - dim 3 (slices)    process A parameter sample
  # - dim 4 (cube rows) process A representation
  adim  <- dim(delta) 
  nA    <- adim[4]
  nB    <- adim[2]
  n     <- adim[1] # adim[1] and adim[3] should be the same

  # permute delta array so that dim order = nA,n,nB,n
  delta <- aperm(delta, 4:1 )
  
  # assume equal probability for each model
  mpA   <- rep(1/nA, nA )
  mpB   <- rep(1/nB, nB )

  # allocate the calculation arrays and matrices
  E_ThetaBMB   <- array(0,dim=c(nA,n,nB))
  E_MB         <- matrix(0,nA,n)
  E_ThetaAMA_2 <- matrix(0,nA)
  E_ThetaAMA   <- matrix(0,nA)
  
  # Loop over the models of the process in question
  for( i in 1:nA ) {  
    
    # Loop over the parameter realizations of the models of the process in question
    for( j in 1:n ) {      
      
      # Loop over the models of the other processes
      for( k in 1:nB ) {
        # construct subscripts matrix for output array
        smat              <- cbind(i,j,k,1:n)
        # Calculate the partial mean E considering uncertain parameters            
        E_ThetaBMB[i,j,k] <- mean(delta[smat])
      }
      
      # Calculate the model averaging of the other process
      E_MB[i,j] <- sum(mpB*E_ThetaBMB[i,j,])
    }
    
    # Calculate the partial mean of the process in question (the mean of each representation of the process in question) considering uncertain parameters    
    E_ThetaAMA_2[i] <- mean(E_MB[i,]^2)
    E_ThetaAMA[i]   <- mean(E_MB[i,])
  }
  
  # Model averaging of each representation of the process in question
  E_MA_2 <- sum(mpA*E_ThetaAMA_2)
  E_MA   <- sum(mpA*E_ThetaAMA)
  
  # Calculate the partial variance VB
  pvar   <- E_MA_2 - (E_MA)^2
  return(pvar)
}


# call resampling and process SA recalculation function
# - designed to be called from a vapply function
boot_psa <- function(nsamp, bootn, ... ) {
  vapply(rep(nsamp,bootn), redsamp_psa, 1, nsamp=nsamp, ... )
}

# resampling and process SA recalculation function
redsamp_psa <- function(nsamp, a1, ... ) {
  ss1   <- ceiling(runif(nsamp)*dim(a1)[1])
  ss2   <- ceiling(runif(nsamp)*dim(a1)[1])
  calc_process_sensitivity(a1[ss2,,ss1,,drop=F])$sensitivity
} 

# functions to convert Process SA list output into a dataframes
convert_to_df_1list_proc1 <- function(onelist, s=NA, proc=NA ) {
  df1 <- unlist(onelist)
  data.frame( scenario=s, variable=proc, t(df1) )
}

convert_to_df_3list_proc <- function(list1) {
  df1 <-
    do.call('rbind',
      lapply(1:length(list1), function(j, l1=list1[[j]] ) 
        do.call('rbind',
          lapply(1:length(l1), function(i, l=l1[[i]], j ) 
          convert_to_df_1list_proc1(l,s=i,proc=names(list1)[j]), j=j ))))
  names(df1)[3:7] <- c('mean','variance','sd','partial_variance','sensitivity1') 
  
  # calculate weighted mean sensitivities across scenarios
  ls    <- length(unique(df1$scenario))
  pscen <- rep(1/ls,ls) 
  df1p  <- subset(df1,scenario==1) 
  df1p$scenario <- -1
  for(p in unique(df1$variable)) {
    df1s <- subset(df1,variable==p)
    df1p[which(df1p$variable==p),3:7] <- apply(as.matrix(df1s[,3:7]), 2, function(v) sum(v*pscen) )
  }

  # calculate integrated sensitvity 
  df1p$sensitivity1 <- df1p$partial_variance / df1p$variance 

  rbind(df1p,df1)
}



# Functions related to the parameter sensitivity index using saltelli method (Saltelli et al 2010).
#####################################################


# wrapper to calculate and average Sobol sensitivity indices for multiple environment and model combinations
calc_parameter_sensitivity <- function(AB, ABi) {
  # Expects the use of multiple models and scenarios

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
    total_var <- sum( sapply(inlist, function(l) l$Tvar[2]) * prob)
    olist <- 
      list(
        Tvar = c(
          mean        = sum( sapply(inlist, function(l) l$Tvar[1]) * prob),
          total_var   = total_var,
          total_sd    = total_var^0.5
        ),
        pvar          = apply( sapply(inlist, function(l) l$pvar) , 1, function(v) sum(v*prob) ),
        pe            = apply( sapply(inlist, function(l) l$pe)   , 1, function(v) sum(v*prob) ),
        sensitivity1  = apply( sapply(inlist, function(l) l$pvar) , 1, function(v) sum(v*prob) )/total_var,
        sensitivityT  = apply( sapply(inlist, function(l) l$pe)   , 1, function(v) sum(v*prob) )/total_var
      )
    
    # name parameter vectors
    if(!is.null(names(inlist[[1]]$sensitivity1))) {
      names(olist$sensitivity1) <- names(inlist[[1]]$sensitivity1)  
      names(olist$sensitivityT) <- names(inlist[[1]]$sensitivity1)  
    }
    
    olist
  }
  
  # calculate number of models and scenarios
  nmod    <- dim(AB)[1]
  nscen   <- dim(AB)[2]

  # assume models are of equal probability
  pmod    <- rep(1/nmod,nmod)
  # assume scenarios are of equal probability
  pscen   <- rep(1/nscen,nscen)
  
  # calculate Sobol for each scenario and model combination
  # model loop
  for(m in 1:nmod) {
    # scenario loop
    for(e in 1:nscen) {
      # out1 <- list(func_sobol_sensitivity(AB[m,e,], ABi[m,e,,] ))      
      out1 <- list(func_sobol_sensitivity(as.numeric(AB[m,e,]), ABi[m,e,,] ))      
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



# Calculates Sobol first order and total sensitivity indices using Saltelli method for a single environment and model combination
func_sobol_sensitivity <- function(ABout,ABiout) {
  
  # expects: 
  # - ABout  - a vector of output from the 2sn sobol mc, i.e. from param matrices A and B
  # - ABiout - an sn * k matrix of output from the k ABi matrices
  
  # returns a list with objects:
  # - Tvar         - a vector of mean (overall mean of model output), total_var (overall variance of model output), total_sd (overall standard deviation of model output)
  # - pvar         - vector of length k, partial variance for each parameter
  # - pe           - vector of length k, partial mean for each parameter
  # - sensitivity1 - vector of length k, first order sensitivity (main effect)
  # - sensitivityT - vector of length k, total sensitivity (including all interactions)

  # calculate dimension extents
  sn <- dim(ABiout)[1]
  k  <- dim(ABiout)[2]
    
  # initialise partial variance and mean vectors
  vv <- numeric(k)
  ev <- numeric(k)
  
  for(p in 1:k) {
    # Calculate the partial variance for a single parameter
    # following the equation (c) in Table 2 Saltelli et al. (2010), from Jansen (1999)
    vv[p] <- sum( (ABout[(sn+1):(2*sn)] - ABiout[1:sn,p])^2 )
    
    # Calculate the partial mean for a single parameter
    # following the equation (f) in Table 2 Saltelli et al. (2010), from Jansen (1999)
    ev[p] <- sum( (ABout[1:sn] - ABiout[1:sn,p])^2 )
  }
  
  # calculate the total variance
  v_t <- var(ABout)
  
  # Calculate first-order sensitivity index (Si)
  # Complete the equation (c) in Table 2 Saltelli et al. (2010), from Jansen (1999)
  si  <- ( v_t - (vv/(2*sn)) ) / v_t
  
  # Calculate total-effect sensitivity index (STi)
  # Complete the equation (f) in Table 2 Saltelli et al. (2010), from Jansen (1999)
  st  <- (ev/(2*sn)) / v_t
  
  # name parameter SA vectors
  pnames <- dimnames(ABi)[[4]]      
  if(!is.null(pnames)) {
    names(si) <- pnames
    names(st) <- pnames
  }
  
  # a model and scenario name would be useful in this output list - would need to be added to the output from the wrapper 
  Tvar <- c(mean=mean(ABout),total_var=v_t,total_sd=v_t^0.5)
  
  # output list
  list(Tvar=Tvar,pvar=v_t - (vv/(2*sn)),pe=ev/(2*sn),sensitivity1=si,sensitivityT=st)
}

# call resampling and parameter SA recalculation function
# - designed to be called from a vapply function
boot_sa <- function(nsamp, bootn, k, ... ) {
  vapply(rep(nsamp,bootn), redsamp_sa, numeric(k), nsamp=nsamp, ... )
}

# resampling and parameter SA recalculation function
redsamp_sa <- function(nsamp, AB, ABi, model=1, env=1, ... ) {
  ss1   <- ceiling(runif(nsamp)*dim(ABi)[3])
  ss2   <- ss1 + dim(ABi)[3]
  calc_parameter_sensitivity(AB[,,c(ss1,ss2),drop=F], ABi[,,ss1,,drop=F] )$individual[[model]][[env]]$sensitivity1
} 

# functions to convert Sobol SA list output into a dataframes
convert_to_df <- function(list1, m=T, s=T ) {
  if(m&s) convert_to_df_1list(list1)
  else if(m|s) convert_to_df_2list(list1, m=m )
  else convert_to_df_3list(list1)
}

convert_to_df_1list <- function(onelist, m=-1, s=-1 ) {
  df1 <- 
    suppressWarnings(
      data.frame( model=m, scenario=s, variable=names(onelist$sensitivity1), mean=onelist$Tvar[1], variance=onelist$Tvar[2], 
                  vapply(onelist[2:5], function(v) v, numeric(length(onelist$pvar))) )
    )
  names(df1)[6:7] <- c('partial_variance','partial_mean')
  df1
}

convert_to_df_2list <- function(list1, m ) {
  do.call('rbind',
          lapply(1:length(list1), function(i,l=list1[[i]]) 
            if(m) convert_to_df_1list(l,s=i)
            else  convert_to_df_1list(l,m=i)
          ))
}

convert_to_df_3list <- function(list1) {
  do.call('rbind',
          lapply(1:length(list1), function(j, l1=list1[[j]] ) 
            do.call('rbind',
                    lapply(1:length(l1), function(i, l=l1[[i]], j ) 
                      convert_to_df_1list(l,m=j,s=i), j=j ))))
}



### END ###
