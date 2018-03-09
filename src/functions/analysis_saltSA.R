################################
#
# MAAT - calculate Saltelli sensitivity index 
# 
# AWalker (walkerap@ornl.gov) 
# October 2017
#
################################

library(pryr)
source('SA_functions.R')



# process and calculate parameter sensitivity analysis 
###################################################################


# load output dataframe from Sobol analysis
setwd(wdd)
AB    <- readRDS(paste(runid,'_out_salt_AB.RDS',sep=''))
ABi   <- readRDS(paste(runid,'_out_salt_ABi.RDS',sep=''))

# identify model output of interest (delta)
delta <- which(dimnames(AB)[[4]]==delta_var)
if(dimnames(ABi)[[4]][delta]!=delta_var) stop('Sobol SA model output subscripts not the same in AB and ABi output arrays')

# need to code tests for the dimensions of lists and arrays are what they are expected to be
# dimensions of AB array from Saltelli ensemble
# - model cominations, environment combinations, parameter samples, output variable
dim(AB)
dimnames(AB)
object_size(AB)

# dimensions of ABi array from Saltelli ensemble
dim(ABi)
dimnames(ABi)
object_size(ABi)

# calculate response in AB and ABi arrays
if(!is.null(response)) {
  # this apply runs the function over the 2nd dimension (environment) of the array
  # the result is an array with the model and environment dimensions reversed compared with the original arrays, this is corrected by aperm
  ABi <- aperm(apply(ABi, c(1,3:5), res_function ), c(2,1,3:5) )
  gc()
  AB  <- aperm(apply(AB,  c(1,3:4), res_function ), c(2,1,3:4) )
}

# subset for model output variable, preserving dimensions
AB    <- AB[,,,delta,  drop=F]
ABi   <- ABi[,,,delta,,drop=F]

# drop model output variable dimension
dn_AB  <- dimnames(AB)
dn_ABi <- dimnames(ABi)
AB     <- array(AB,  dim(AB)[1:3]  )
ABi    <- array(ABi, dim(ABi)[c(1:3,5)] )
dimnames(AB)  <- dn_AB[1:3]
dimnames(ABi) <- dn_ABi[c(1:3,5)]

# Sobol for absolute values
sms   <- calc_parameter_sensitivity(AB, ABi )

# save sobol analysis
setwd(wdt)
saveRDS(sms, paste(runid_out,delta_var,'salt_list.RDS',sep='_') )

# remove large data structures
rm(list=ls(pattern='^AB'))
gc()



### END ###