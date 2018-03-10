################################
#
# MAAT - analysis
# 
# AWalker (walkerap@ornl.gov) 
# March 2016
#
################################



# Currently this is a stand alone script to be run after the MAAT simulation
# this analysis script requires some user defined setup, e.g. in this case the making of co2 response data frames which require post-processing of model output 
# so this script is currently run line by line from within RStudio or similar
# but this script could be generalised to be called from run_MAAT.R
# this would allow the command line options of run_MAAT.R to be used by this script, and some of the lists in init_MAAT.R to inform the analysis 



# initialise analysis
###################################################################

rm(list=ls())

srcdir <- '##SRCDIR##' 
wd     <- '##PDIR##'



# set up run specific variables
###################################################################

runid      <- 'SA'
odate      <- '##OUTPUTDATE##'
  
# routine switches
salt       <- T    # post-process and calculate Saltelli parameter sensitivity index 
proc       <- T    # post-process and calculate Ye process sensitivity index
boot       <- T    # boot strap the above that are TRUE
plot_SA    <- T    # plot sensitvity indices 
plot_mout  <- F    # plot model output
plot_pars  <- F    # plot model out vs parameter values - only implemented if plot_mout is TRUE

# output of interest
# - 'h4' is hydraulic head at 6000 m in the model domain. 
# - the label four is because it is the fourth hyraulic head output from the model (i.e. h is not output across tyhe whole model domain)
# - the subset of the domain is requested as output in init_MAAT_SA.R with the variable out_hsub in the pars.static list
delta_var <- 'h4'

# response function for delta_var
response   <- NULL

# two vectors as long as the number of processes in process sensitivity analysis
# names of the processes
proc_names <- c('Recharge','Geology')
# number of process representations for each process in question 
n_procs    <- c(2,2)

# number of MC realizations for each process in process SA method 
n          <-  1000
# number of MC realizations for each parameter in parameter SA method 
sobol_n    <-  1000*n

# boot variables (only used if boot is TRUE)
# number of resamples per sample n
bootn        <- 20
# vector of sample n for process SA
resamp_v_psa <- c(1:10, seq(10,100,10), seq(200,n,100) )
# vector of sample n for parameter SA
resamp_v_sa  <- c(seq(1e1,1e2,1e1), seq(1e3,1e4,1e3), seq(1e5,sobol_n,1e5) )
# subscript for environmental condition to run the bootstrap on
e            <- 1
# subscript for model combination to run the bootstrap on (only relevent to parameter SA boot)
m            <- 1



# ensemble labels & plotting parameters
###################################################################

# vector of parameter names and their names (expressions) for plotting
par_names  <- c(
  gwater_rt.a  = expression(a),
  gwater_rt.b  = expression(b),
  gwater_rt.K  = expression(K),
  gwater_rt.K1 = expression(K[1]),
  gwater_rt.K2 = expression(K[2])
)

# names of the models (combinations of process representation)
mod_names <- c(
  f_rechrg_lin.f_trans_single   = expression(R[2]*G[1]),
  f_rechrg_power.f_trans_single = expression(R[1]*G[1]),
  f_rechrg_lin.f_trans_double   = expression(R[2]*G[2]),
  f_rechrg_power.f_trans_double = expression(R[1]*G[3])
)


# sensitivity index plotting parameters
##############################

# environment labels
evar1name <- 'Precipitation'
evar2name <- NULL
evar1     <- 1524
evar2     <- NULL

# type of sensitivty to plot (so far only relevant for parameters)
sens      <- 'sensitivity1'



# model output plotting parameters
##############################

# suscripts for environment subsetting
esubs <- NULL

# conditioning variable
condvar1 <- NULL

# axis names cex
acex  <- 1.5

# points cex
pcex  <- 0.8*acex

# axis labels
ylab    <- list(expression('Hydraulic Head [m]'), cex=acex )
xlab_e1 <- list(expression('Precipitation [mm]'), cex=acex )
xlab_e2 <- NULL

# axis extents
ylim    <- NULL
xlim    <- NULL

# y-axis gridlines
gls     <- NULL

# output id
output_label <- 'absolute'



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

# # plot model output
# setwd(paste0(srcdir,'/functions'))
# if(plot_mout) source('plot_modout.R')



### END ###
