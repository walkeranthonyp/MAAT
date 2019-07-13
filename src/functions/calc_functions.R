################################
#
# Mathematical functions for MAAT model 
# 
# AWalker January 2017
#
################################

suppressPackageStartupMessages(library(randtoolbox))



#####################################################
res_calc <- function(df1,df2,rcols,type='abs') {
  dfo <- df1
  
  if(class(df1)=='list') {
    # if a list, expects rcols to be a single integer that indexes a vector/matrix element of the input 
    if(type=='abs')      dfo[[rcols]] <- df1[[rcols]] - df2[[rcols]]
    else if(type=='rel') dfo[[rcols]] <- df1[[rcols]] / df2[[rcols]]
    
  } else {
    # if a dataframe or matrix, expects rcols to be a vector of integer(s) that indexes columns of the input
    if(type=='abs')      dfo[,rcols] <- df1[,rcols] - df2[,rcols]
    else if(type=='rel') dfo[,rcols] <- df1[,rcols] / df2[,rcols]
    
  }# else stop(paste('no method for class:',class(df1))) 
  dfo
}



#####################################################
rsobol <- function(n,pars=c(1,0.1),norm=F) {
  # this function generates the Sobol sequence for either 
  # a normal distribution norm=T   or  a uniform distribution norm=F
  # the pars argument is a two element vector,
  # for normal distributions the first element is the mean, the second the coefficient of variation
  # for uniform distributions the first element is the lower value of the range, the second the upper value of the range
 
  return( 
    if(norm) {
      pars[1] * ((sobol(n,init=maat$wpars$sobol_init,normal=T) / (0.674489/pars[2]) ) + 1)
    } else {
      (pars[2]-pars[1]) * sobol(n,init=maat$wpars$sobol_init,normal=F) + pars[1]
    }
  )
  
  maat$wpars$sobol_init <- F 
}



#####################################################
runif1 <- function(n,pars=c(1,0.1),norm=F) {
  # this function generates a sequence for either 
  # a normal distribution norm=T   or  a uniform distribution norm=F
  # the pars argument is a two element vector,
  # for normal distributions the first element is the mean, the second the coefficient of variation as a proportion
  # for uniform distributions the first element is the lower value of the range, the second the upper value of the range
 
  return( 
    if(norm) {
      pars[1] * rnorm(n,1,pars[2]) 
    } else {
      runif(n,min(pars),max(pars))
    }
  )
  
  maat$wpars$sobol_init <- F 
}



#####################################################



### END ###
