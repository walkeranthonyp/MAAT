################################
#
# MAAT - write tables and plot sensitivity analysis (process and parameter) 
# 
# AWalker (walkerap@ornl.gov) 
# March 2017
#
################################

setwd(paste0(srcdir,'/functions'))
source('plotting_functions.R')
source('general_functions.R')



# environment labels
elabs <- 
  if(!is.null(evar1)&!is.null(evar2)) {
    paste(paste0(evar1name,'_',evar1),paste0(evar2name,'_',evar2))
  } else if(!is.null(evar1)) {
    paste0(evar1name,'_',evar1)
  } else {
    NULL
  }

# model labels
setwd(wdd)
fname <- paste('out',runid,'salt','AB','dataf',sep='_')
l1    <- readRDS(paste(fname,'.RDS',sep=''))
full_mod_names <- apply(l1$fnames,1,paste,collapse='.')
ss_names       <- if(!is.null(mod_names)) match(full_mod_names,names(mod_names)) else NULL
mlabs          <- if(!is.null(mod_names)) mod_names[ss_names] else full_mod_names
mlabs          <- if(length(mlabs)>12) paste0('M',1:length(mlabs)) else mlabs  



# plot process SA indices
############################################

# read process SA index list
setwd(wdt)
df1   <- read.csv(paste(runid_out,delta_var,'psa.csv',sep='_'))
df1t  <- sens_table(df1)

# add environmetal variables to sensitivity output matrix
df1te <- add_scenario_values(df1t)
  
# write latex table, & create pdf (if call=TRUE)
write_Latex(df1te, paste(runid_out,'procSAtables',sep='_'), write_table_Latex, call=T )

# integrated sensitivities
df1$averaging <- ifelse(df1$scenario==-1, 'Combined', 'Scenario' )

# plot
setwd(wdp)
pdf(paste(paste(runid_out,'PSAplots','radincs',sep='_'),'.pdf',sep=''), width=3.5, height=7)
# placeholder for if statement & function if under 4 variables
# plot_SA_4orless(yv=df1$sensitivity1, xv=df1$co2, cv1=df1$par, gv=df1$variable,ylim=ylim1, ylab=ylab1, xlab=xlab_e1 )
radar_plot(df1, var1='averaging', vnames=NULL, lnames=list(NULL, elabs), max_si=1.0 )
dev.off()



# parameter sensitivity plots 
############################################

# read data
setwd(wdt)
df1  <- read.csv(paste(runid_out,delta_var,'salt.csv',sep='_'))
df1t <- sens_table(df1)

# add environmetal variables to sensitivity output matrix
df1te <- add_scenario_values(df1t)

# add par names suitable for latex
stf <- NULL
if(!is.null(par_names_latex)) {
  ss_names <- match(names(df1te),names(par_names_latex))
  labs     <- par_names_latex[ss_names]
  names(df1te)[which(!is.na(labs))] <- labs[which(!is.na(labs))]
  stf <- function(x){x}
}

# write latex table, & create pdf (if call=TRUE)
# write_Latex(df1te, paste(runid_out,'saltSAtables',sep='_'), write_table_Latex, call=T, landscape=T, sanitize.text.function=stf )
# subset to individual tables
tlist <- list(
  subset(df1te, par=='int'|model=='int' ),
  subset(df1te, par!='int'&model!='int' )
)
write_Latex(tlist, paste(runid_out,'saltSAtables',sep='_'), write_tablelist_Latex, call=T, 
            landscape=T, sanitize.text.function=stf, floating=F, tabular.environment="longtable" )

# legend labels
lnames <- list(NULL, elabs, mlabs )

# integrated sensitivities
df1$averaging <- ifelse(df1$model==-1&df1$scenario==-1, 'Combined', ifelse(df1$scenario==-1, 'Model', ifelse(df1$model==-1, 'Scenario', NA )))
df1c <- subset(df1, scenario==-1|model==-1 )

# broken out by scenario within each model
df1s <- subset(df1, scenario!=-1&model!=-1 )

# plot
setwd(wdp)
pdf(paste(paste(runid_out,'SAplots','radincs',sep='_'),'.pdf',sep=''), width=10, height=7)
radar_plot(df1c, var1='averaging', lnames=lnames, max_si=1.0 )
radar_plot(df1s, var1='model', lnames=list(elabs), max_si=1.0 )
dev.off()

df1s <- subset(df1, scenario!=-1&model!=-1 )
radar_plot(df1s, var1='model', lnames=lnames, max_si=1.0 )





# plot output variability 
############################################

setwd(wdd)
fname <- paste('out',runid,'salt','AB',sep='_')
AB    <- readRDS(paste(fname,'.RDS',sep=''))

ppp  <- plot_env_array(AB, yvar=delta_var, condvar1='lim', ylim=ylim, ylab=ylab, xlim=xlim, lout=c(3,2), gls=gls )
ppp  <- plot_env_array(AB, yvar=delta_var, condvar1='lim', ylim=ylim, ylab=ylab, xlim=xlim, xlab=xlab, lout=c(3,2), gls=gls )
ppp  <- plot_env_array(AB, yvar=delta_var, condvar1='lim', ylim=ylim, ylab=ylab, xlim=xlim, lout=c(3,2), gls=gls )

setwd(wdp)
pdf(paste(runid,delta_var,'.pdf',sep='_'),width=10,height=6);                  print(ppp[[2]]); dev.off()
pdf(paste(runid,delta_var,'_bymodel.pdf',sep='_'),width=10,height=6);          print(ppp[[1]]); dev.off()
pdf(paste0(paste(runid,delta_var,condvar1,sep='_'),'.pdf'),width=10,height=6); print(ppp[[3]]); dev.off()



### END ###
