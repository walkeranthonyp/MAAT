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
source('general_functions.R')

runid  <- 'fixedr_fixedVcmaxJe_sobol'

# number of representations/models in processes A & B
nA     <- 1                 
nB     <- 3                 

# number of MC realizations for each process It is assumed that each model
# - uses the same number of realizations, but this is not required for the method 
n       <-  1000
sobol_n <-  100*n

# set working directory
setwd(wd)
setwd('./results/2016-03-16')



# Sobol analysis
###################################################################

# total number of parameters varied in the Sobol
pars  <- 6 
# output of interest (column index for the output matrix)
delta <- 1 

# load output dataframe from Sobol analysis
list1 <- readRDS(paste(runid,'_out_sobol.RDS',sep=''))
# names(list1)
# dim(list1$AB)
# dim(list1$ABi[[1]][[1]])

# Sobol for absolute values
sms   <- calc_parameter_sensitivity(sobol_n,pars,list1,delta)
# sms[2]
# names(sms[[2]][[1]])

# create response data for Sobol
ABtest       <- list1$AB[,1:2,,]
ABtest[,1,,] <- list1$AB[,2,,]-list1$AB[,1,,]
ABtest[,2,,] <- list1$AB[,3,,]-list1$AB[,2,,]

ABitest      <- list(
  lapply(list1$ABi,function(l) l[[2]]-l[[1]]),
  lapply(list1$ABi,function(l) l[[3]]-l[[2]])
  )

for(m in 1:3) {
  l1 <- lapply(ABitest,function(l) l[[m]] )
  ABitest_t <- if(m==1) list(l1) else c(ABitest_t,list(l1))
}

length(ABitest_t)
length(ABitest_t[[3]])

list1_ares <- list(AB  = ABtest,
                   ABi = ABitest_t) 

# Sobol for response values
sms_ares <- calc_parameter_sensitivity(sobol_n,pars,list1_ares,delta)
# names(sms_ares)
# sms_ares[[2]]

# save sobol analysis
saveRDS(list(absolute=sms,response=sms_ares),paste(runid,'sobol_list.RDS',sep='_'))

# remove large lists
rm(ls(pattern='^list1'))



# process sensitivity analysis
###################################################################

# load data
df1   <- readRDS(paste(runid,'_out.RDS',sep=''))
names(df1)
res_cols <- 9:15
for(c in res_cols) df1[,c] <- as.numeric(df1[,c])
df1$leaf.ca_conc <- as.numeric(as.character(df1$leaf.ca_conc))

# subset out various environmental conditions - i.e. CO2 at 280, 400, 600
df_280 <- df1[seq(1,length(df1[,1]),3),]
df_400 <- df1[seq(2,length(df1[,1]),3),]
df_600 <- df1[seq(3,length(df1[,1]),3),]

# calculate responses to CO2 increase
df_400_280_ares <- res_calc(df_400,df_280,res_cols)
df_600_400_ares <- res_calc(df_600,df_400,res_cols)
df_400_280_res  <- res_calc(df_400,df_280,res_cols,type='rel')
df_600_400_res  <- res_calc(df_600,df_400,res_cols,type='rel')

# cobined response dataframe - currently for plotting purposes only
df_ares <- data.frame(rbind(df_400_280_ares,df_600_400_ares),res=rep(c('280-400','400-600'),each=length(df_400_280_ares[,1])))

# calculate variances and sensitivities
a280     <- calc_process_sensitivity(df_280,'A',n)
a400     <- calc_process_sensitivity(df_400,'A',n)
a600     <- calc_process_sensitivity(df_600,'A',n)

ar400280 <- calc_process_sensitivity(df_400_280_ares,'A',n)
ar600400 <- calc_process_sensitivity(df_600_400_ares,'A',n)

rr400280 <- calc_process_sensitivity(df_400_280_res, 'A',n)
rr600400 <- calc_process_sensitivity(df_600_400_res, 'A',n)

# set up output list
olist <- list(abs_280=a280,abs_400=a400,a600=a600,
              absres_400_280=ar400280,absres_600_400=ar600400,
              relres_400_280=rr400280,relres_600_400=rr600400)

# write output
saveRDS(olist,paste(runid,'psa_list.RDS',sep='_'))

# plot output
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


