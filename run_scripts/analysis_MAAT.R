################################
#
# MAAT Leaf Model - analysis
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
# however this script would require modification at the time of 



# initialise analysis
###################################################################

rm(list=ls())

srcdir <- "/mnt/disk2/models/maat/src" 
wd     <- "#PROJECTDIR#"
setwd(srcdir)
source('SA_functions.R')
source('general_functions.R')

runid  <- "#USRDEFINED#"

# this needs to be brought into the .RDS object and read from there, or from the init file
# a vector as long as the number of processes in process sensitivity analysis
# - each element is the number of process representations for the process in question 
n_procs <- c(2,3)

# number of MC realizations for each process It is assumed that each model
# - uses the same number of realizations, but this is not required for the method 
n       <-  100
sobol_n <-  1*n

# set working directory
setwd(wd)
setwd('./results/2017-01-18')



# Sobol analysis
###################################################################

# total number of parameters varied in the Sobol
pars  <- 7 
# output of interest (column index for the output matrix)
delta <- 1 

# load output dataframe from Sobol analysis
list1 <- readRDS(paste(runid,'_out_salt.RDS',sep=''))
# dimensions of AB array from Saltelli ensemble
# - model cominations, environment combinations, parameter samples, output variable
dim(list1$AB)

# dimensions of AB array from Saltelli ensemble
# - model cominations
length(list1$ABi)
# - environment combinations
length(list1$ABi[[1]])
# - parameter samples, parameter, output variable
dim(list1$ABi[[1]][[1]])

# Sobol for absolute values
sms   <- calc_parameter_sensitivity(sobol_n,pars,list1,delta)


# create response data for Sobol
# need to automate this process
# second dim of this array is the env variable

# response in AB matrix
ABtest       <- list1$AB[,1:2,,]
ABtest[,1,,] <- list1$AB[,2,,]-list1$AB[,1,,]
ABtest[,2,,] <- list1$AB[,3,,]-list1$AB[,2,,]

ABitest      <- list(
  lapply(list1$ABi,function(l) l[[2]]-l[[1]]),
  lapply(list1$ABi,function(l) l[[3]]-l[[2]])
)

# response in ABi matrix
# - I'm not sure this is correct, 
for(m in 1:prod(n_procs)) {
  l1 <- lapply(ABitest,function(l) l[[m]] )
  ABitest_t <- if(m==1) list(l1) else c(ABitest_t,list(l1))
}

dim(ABitest_t[[1]][[1]])
dim(ABitest_t)

length(ABitest_t)
length(ABitest_t[[3]])

list1_ares <- list(AB  = ABtest,
                   ABi = ABitest_t) 

# Sobol for response values
sms_ares <- calc_parameter_sensitivity(sobol_n,pars,list1_ares,delta)
# names(sms_ares)
# sms_ares[[2]]

# save sobol analysis
saveRDS(list(absolute=sms,response=sms_ares),paste(runid,'salt_list.RDS',sep='_'))

# remove large lists
rm(ls(pattern='^list1'))



# process sensitivity analysis
###################################################################

filter_data <- function(df1,...) {
  # subset out various environmental conditions - i.e. CO2 at 280, 400, 600
  df1a <- df1[seq(1,length(df1[,1]),3),]
  df1b <- df1[seq(2,length(df1[,1]),3),]
  df1c <- df1[seq(3,length(df1[,1]),3),]
  
  # calculate responses to CO2 increase
  df1d <- res_calc(df1b,df1a,res_cols)
  df1e <- res_calc(df1c,df1b,res_cols)
  df1f <- res_calc(df1b,df1a,res_cols,type='rel')
  df1g <- res_calc(df1c,df1b,res_cols,type='rel')
  
  list(df_280=df1a,df_400=df1b,df_600=df1c,df_400_280_ares=df1d,df_600_400_ares=df1e,df_400_280_res=df1f,df_600_400_res=df1g)
}

# process loop
p <- 1
for(p in 1:length(n_procs) ) {
  
  # load data
  df1   <- readRDS(paste(runid,'_out_proc_',p,'.RDS',sep=''))
  names(df1)
  res_cols <- 11:17
  for(c in res_cols) df1[,c] <- as.numeric(df1[,c])
  df1$leaf.ca_conc <- as.numeric(as.character(df1$leaf.ca_conc))
  
  # generate subsets based on environment
  dlist  <- filter_data(df1,res_cols=res_cols)
  
  # combined response dataframe - currently for plotting purposes only
  df_ares <- data.frame(rbind(dlist$df_400_280_ares,dlist$df_600_400_ares),res=rep(c('280-400','400-600'),each=length(dlist$df_400_280_ares[,1])))
  
  # calculate output
  salist <- lapply(dlist,calc_process_sensitivity,colname='A',n=length(n_procs),nA=n_procs[p],nB=prod(n_procs[-p]))
  
  # all processes
  olist  <- if(p==1) list(salist) else c(olist,list(salist))
}

# write output
saveRDS(olist,paste(runid,'psa_list.RDS',sep='_'))





# plot output
###################################################################

library(lattice)
ylab  <- expression('A ['*mu*mol*m^-2*s^-1*']')
ylabr <- expression('A response ['*mu*mol*m^-2*s^-1*']')
xlab  <- expression('Atmospheric CO'[2]*' concentration ['*mu*mol*mol^-1*']')
xlabr <- expression('CO'[2]*' increase ['*mu*mol*mol^-1*']')

p1 <- 
xyplot(A~leaf.ca_conc,df1,ylab=ylab,xlab=xlab,xlim=c(100,800),ylim=c(0,21),
       panel=function(...,box.width=100){
         panel.violin(...,horizontal=F,col = "lightblue",varwidth = FALSE, box.width = box.width)
         panel.bwplot(...,horizontal=F, col='black',coef=100,cex=0.8,pch='|',fill='gray',box.width = 0.1*box.width)
       },
       par.settings = list(box.rectangle=list(col='black'),box.umbrella=list(alpha=0,col='black',lty=1),
                           plot.symbol = list(pch='.', cex = 0.1))
)

p2 <- 
xyplot(A~leaf.ca_conc|leaf.etrans,df1,ylab=ylab,xlab=xlab,xlim=c(100,800),ylim=c(0,21),layout=c(3,1),
       panel=function(...,box.width=100){
         panel.violin(...,horizontal=F,col = "lightblue",varwidth = FALSE, box.width = box.width)
         panel.bwplot(...,horizontal=F, col='black',coef=100,cex=0.8,pch='|',fill='gray',box.width = 0.1*box.width)
       },
       par.settings = list(box.rectangle=list(col='black'),box.umbrella=list(alpha=0,col='black',lty=1),
                           plot.symbol = list(pch='.', cex = 0.1))
)


# response
p1r <- 
  bwplot(A~res,df_ares,ylab=ylabr,xlab=xlabr,ylim=c(0,7),
         horizontal=F,
         panel=function(...,box.width=1){
           panel.violin(...,col = "lightblue",varwidth = FALSE, box.width = box.width)
           panel.bwplot(..., col='black',coef=100,cex=0.8,pch='|',fill='gray',box.width = 0.1*box.width)
         },
         par.settings = list(box.rectangle=list(col='black'),box.umbrella=list(alpha=0,col='black',lty=1),
                             plot.symbol = list(pch='.', cex = 0.1))
  )

p2r <- 
  bwplot(A~res|leaf.etrans,df_ares,ylab=ylabr,xlab=xlabr,ylim=c(0,7),
         horizontal=F,layout=c(3,1),
         panel=function(...,box.width=1){
           panel.violin(...,col = "lightblue",varwidth = FALSE, box.width = box.width)
           panel.bwplot(..., col='black',coef=100,cex=0.8,pch='|',fill='gray',box.width = 0.1*box.width)
         },
         par.settings = list(box.rectangle=list(col='black'),box.umbrella=list(alpha=0,col='black',lty=1),
                             plot.symbol = list(pch='.', cex = 0.1))
  )


pdf(paste(runid,'plot_overall_res.pdf',sep='_'),width=12,height=6)
print(update(p1,main='absolute'),split=c(1,1,2,1),more=T)
print(update(p1r,main='response'),split=c(2,1,2,1),more=F)
dev.off()

pdf(paste(runid,'plot_etrans_res.pdf',sep='_'),width=12,height=12)
print(update(p2,main='absolute'),split=c(1,1,1,2),more=T)
print(update(p2r,main='response'),split=c(1,2,1,2),more=F)
dev.off()

# remove large datasets
rm(ls(pattern='^df'))


