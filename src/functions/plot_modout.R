################################
#
# MAAT - plot model output 
# 
# AWalker (walkerap@ornl.gov) 
# March 2017
#
################################

source('plotting_functions.R')



# load data
setwd(wdd)
fname <- paste(runid,'out','salt','AB',sep='_')
AB    <- readRDS(paste(fname,'.RDS',sep=''))
fname <- paste(runid,'out','salt','AB','dataf',sep='_')
l1    <- readRDS(paste(fname,'.RDS',sep=''))

# calculate response in AB and ABi arrays
if(!is.null(response)) {
  # this apply runs the function over the 2nd dimension (environment) of the array
  # the result is an array with the model and environment dimensions reversed compared with the original arrays, this is corrected by aperm
  AB  <- aperm(apply(AB,  c(1,3:4), res_function ), c(2,1,3:4) )
}

# subset array according to required environmental conditions
if(!is.null(esubs)) {
  AB   <- AB[,esubs,,]
  xlim <- 0:(length(esubs)+1)
}

# plot output variability
# ppp  <- plot_env_array(AB, yvar=delta_var, condvar1=condvar1, ylim=ylim, ylab=ylab, xlim=xlim, xlab=xlab, lout=c(3,2), gls=gls )
ppp  <- plot_env_array(AB, yvar=delta_var, condvar1=condvar1, ylim=ylim, ylab=ylab, xlim=xlim, lout=c(3,2), gls=gls )

setwd(wdp)
pdf(paste0(paste(runid_out,delta_var,sep='_'),'.pdf'),width=10,height=6);          print(ppp[[1]]); dev.off()
pdf(paste(runid_out,delta_var,'_bymodel.pdf',sep='_'),width=10,height=6);          print(ppp[[2]]); dev.off()
if(!is.null(condvar1)) pdf(paste0(paste(runid_out,delta_var,condvar1,sep='_'),'.pdf'),width=10,height=6); print(ppp[[3]]); dev.off()



# plot parameter sensitivities
setwd(wdp)
if(plot_pars) {
  
  # subscripts for which parameters to plot
  pp_subs <- if(is.null(pp_subs)) 1:dim(l1$pars)[2] else pp_subs

  pdf(paste(paste(runid,delta_var,'vs','pars',paste(pp_subs,collapse='-'),sep='_'),'.pdf',sep=''),width=10,height=10)
  for( pp in pp_subs ) {
    print(pp)
    
    if(is.null(pm_subs)) {
      
      cond1 <- rep(rep(l1$env[,1], each=dim(AB)[1]), dim(AB)[3] )
      cond2 <- rep(rep(l1$env[,2], each=dim(AB)[1]), dim(AB)[3] )
      ppp <- 
        plot_parDens(AB[,,,1], rep(l1$pars[,pp], each=prod(dim(AB)[1:2]) ),
                     cond1, cond2,
                     ylim=ylim, ylab=ylab, xlab=par_names[which(colnames(l1$pars)[pp]==names(par_names))],
                     gls=gls)
      print(ppp, newpage = T)
      
    } else {
      
      for( pm in pm_subs ) {
        cond1 <- rep(l1$env[,1], dim(AB)[3] )
        cond2 <- rep(l1$env[,2], dim(AB)[3] )
        
        ppp <- 
          plot_parDens(AB[pm,,,1], rep(l1$pars[,pp], each=dim(AB)[2] ),
                       cond1, cond2,
                       ylim=ylim, ylab=ylab, xlab=par_names[which(colnames(l1$pars)[pp]==names(par_names))],
                       gls=gls)
        
        print(ppp, newpage = T)
      }
    }
  }
  dev.off()

}

# remove large datasets
rm(AB,l1)
gc()



### END ###