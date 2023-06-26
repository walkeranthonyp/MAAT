################################
#
# MAAT - analysis of Saltelli convergence 
# 
# AWalker (walkerap@ornl.gov) 
# October 2017
#
################################

library(pryr)
library(lattice)
source('SA_functions.R')



# process and calculate bootstrapped parameter sensitivity analysis 
###################################################################

# load output dataframe from Sobol analysis
setwd(wdd)
AB    <- readRDS(paste0('out_',runid,'_salt_AB.RDS',sep=''))
ABi   <- readRDS(paste0('out_',runid,'_salt_ABi.RDS',sep=''))

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

# bootstrap sensitivity calculation 
psa    <- vapply(resamp_v_sa, boot_sa, matrix(1,nrow=dim(ABi)[4],ncol=bootn), bootn=bootn, k=dim(ABi)[4], AB=AB, ABi=ABi, model=m, env=e )
olist  <- list(model=dn_AB[[1]][m], data=psa )

# save Saltelli Sobol bootSA analysis
setwd(wdt)
saveRDS(olist, paste0(paste(runid_out,delta_var,'saltBoot',sep='_'),'.RDS'))

# plots
setwd(wdp)
pdf(paste0(paste(runid_out,delta_var,'saltBoot',sep='_'),'.pdf'),height=4)
for(p in 1:dim(psa)[1] ) {
  p1 <- 
    xyplot(as.vector(t(psa[p,,]))~rep(resamp_v_sa,bootn),
           ylab=paste0('Si (',dn_ABi[[5]][p],')'), xlab='n',
           ylim=c(-0.02,1),
           panel=function(...){
             panel.abline(h=0:1, col='grey80' )
             panel.xyplot(...)
             panel.lines(apply(psa[p,,],2,mean)~resamp_v_sa, col='red', lwd=2 )
           })
  p2 <- 
    xyplot(apply(psa[p,,],2,var)^0.5~resamp_v_sa, ylab=paste0('sd in Si (',dn_ABi[[5]][p],')'), xlab='n', abline=list(h=c(0.001,0.01),lty=2), ylim=c(0,0.05) )
  print(p1)
  print(p2)
}
dev.off()

# remove large data structures
rm(list=ls(pattern='^AB'))
gc()



### END ###