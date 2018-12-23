################################
#
# MAAT - calculate Ye (Dai et al 2017) process sensitvity index 
# 
# AWalker (walkerap@ornl.gov) 
# March 2017
#
################################

library(pryr)
source('SA_functions.R')



# process and calculate process sensitivity analysis 
###################################################################

# allocate output list
salist        <- vector(mode="list", length=length(n_procs) )
names(salist) <- proc_names 

# process loop
for( p in 1:length(n_procs) ) {
  
  # set data directory
  setwd(wdd)
  
  # load data
  fname <- paste0('out_',runid,'_proc_',p)
  a1    <- readRDS(paste(fname,'.RDS',sep=''))
  delta <- which(dimnames(a1)[[1]]==delta_var)
  
  dim(a1)
  dimnames(a1)
  object_size(a1)

  # subset for delta, preserving dimensions
  a1    <- a1[delta,,,,,,drop=F]
  
  # drop model output variable dimension
  a1    <- array(a1, dim(a1)[2:6] )
  
  # calculate responses
  if(!is.null(response)) {
    a1  <- apply(a1, 2:5, res_function )
  }

  # calculate sensitivity under different environmental conditions 
  salist[[p]] <- apply(a1, 1, calc_process_sensitivity )

  # remove large datasets
  rm(a1); gc()
}

# save Ye SA analysis
setwd(wdt)
saveRDS(salist, paste0(paste(runid_out,delta_var,'psa_list',sep='_'),'.RDS'))

# create data frames of Process stats
salist_df <- convert_to_df_3list_proc(salist)
write.csv(salist_df, paste(runid_out,delta_var,'psa.csv',sep='_'), quote=F, row.names=F )



### END ###