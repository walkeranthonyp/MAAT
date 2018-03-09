################################
#
# MAAT - bootstrap to test convergence of Ye (Dai et al 2017) process sensitvity index
# 
# AWalker (walkerap@ornl.gov) 
# March 2017
#
################################

library(pryr)
library(lattice)
source('SA_functions.R')



# process and calculate bootstrapped process sensitivity analysis 
###################################################################

# allocate output list
salist        <- vector(mode="list", length=length(n_procs) )
names(salist) <- proc_names 

# process loop
for( p in 1:length(n_procs) ) {
  
  # set data directory
  setwd(wdd)
  
  # load data
  fname <- paste(runid,'_out_proc_',p,sep='')
  a1    <- readRDS(paste(fname,'.RDS',sep=''))
  delta <- which(dimnames(a1)[[1]]==delta_var)
  
  dim(a1)
  dimnames(a1)
  object_size(a1)

  # subset for model output variable and environment, preserving dimensions
  a1    <- a1[delta,e,,,,,drop=F]
  
  # drop model output variable and environment dimensions
  a1    <- array(a1, dim(a1)[3:6] )
  
  # bootstrap sensitivity calculation 
  psa     <- vapply(resamp_v_psa, boot_psa, rep(1,bootn), bootn=bootn, a1=a1 )
  
  # plots
  ppsi <- 
    xyplot(as.vector(t(psa))~rep(resamp_v_psa,bootn),
           ylab=paste0('Sk (',proc_names[p],')'), xlab='n',
           ylim=c(-0.02,1),
           panel=function(...){
             panel.abline(h=0:1, col='grey80' )
             panel.xyplot(...)
             panel.lines(apply(psa,2,mean)~resamp_v_psa, col='red', lwd=2 )
           })
  
  pvar   <- xyplot(apply(psa,2,var)^0.5~resamp_v_psa, ylab=paste0('sd in Sk (',proc_names[p],')'), xlab='n', ylim=c(0,0.05), abline=list(h=c(0.001,0.01),lty=2) )

  # output list
  salist[[p]]  <- list(pname=proc_names[p], data=psa, ppsi=ppsi, pvar=pvar )

  # remove large datasets
  rm(a1); gc()
}

# save Ye bootSA analysis
setwd(wdt)
saveRDS(salist, paste0(paste(runid_out,delta_var,'psaBoot',sep='_'),'.RDS'))

setwd(wdp)
pdf(paste0(paste(runid_out,delta_var,'psaBoot',sep='_'),'.pdf'),height=4)
for( p in 1:length(n_procs) ) {
  psa <- salist[[p]]$data
  print(salist[[p]]$ppsi)
  print(salist[[p]]$pvar)
}
dev.off()









