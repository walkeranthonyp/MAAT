################################
#
# MAthematical functions for MAAT model 
# 
# AWalker January 2017
#
################################

library(randtoolbox)


#####################################################
res_calc <- function(df1,df2,rcols,type='abs') {
  dfo <- df1
  if(type=='abs')      dfo[,rcols] <- df1[,rcols] - df2[,rcols]
  else if(type=='rel') dfo[,rcols] <- df1[,rcols] / df2[,rcols]
  else { print(paste('no method for type:',type)); stop }
  dfo
}



#####################################################
rsobol <- function(n,pars=c(1,0.1),norm=F) {
  # this function generates the Sobol sequence for either 
  # a normal distribution norm=T   or  a uniform distribution norm=F
  # the pars argument is a two element vector,
  # for normal distributions the first element is the mean, the second the coefficient of variation
  # for uniform distributions the first element is the lower value of the range, the second the upper value of the range
  
  if(norm) {
    pars[1] * ((sobol(n,normal=T) / (0.674489/pars[2]) ) + 1)
  } else {
    (pars[2]-pars[1]) * sobol(n,normal=F) + pars[1]
  }
}

