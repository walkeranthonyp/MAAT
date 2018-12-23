################################
#
# MAAT - sensitvity analysis workflow 
# 
# AWalker (walkerap@ornl.gov) 
# March 2018
#
################################



# make working directories & file name strings
###################################################################

wdd <- paste(wd,'/results/',odate,sep='')
wdt <- paste0(wdd,'/tables/')
wdp <- paste0(wdd,'/plots/')
if(!file.exists(wdt)) dir.create(wdt)
if(!file.exists(wdp)) dir.create(wdp)

runid_out    <- paste(runid,output_label,sep='_')



# sensitivity analysis
###################################################################

# parameter SA
setwd(paste0(srcdir,'/functions'))
if(salt)      source('analysis_saltSA.R')
setwd(paste0(srcdir,'/functions'))
if(salt&boot) source('analysis_saltSAboot.R')

# process SA
setwd(paste0(srcdir,'/functions'))
if(proc)      source('analysis_procSA.R')
setwd(paste0(srcdir,'/functions'))
if(proc&boot) source('analysis_procSAboot.R')



# sensitivity analysis & model output plots
###################################################################

# plot SA metrics
setwd(paste0(srcdir,'/functions'))
if(plot_SA)   source('plot_SA.R')

# plot model output
setwd(paste0(srcdir,'/functions'))
if(plot_mout) source('plot_modout.R')

# run homemade plotting script 
setwd(wd)
if(plot_home) source(paste0('plot_homemade_',runid,'.R'))



### END ###