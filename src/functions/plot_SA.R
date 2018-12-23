################################
#
# MAAT - write tables and plot sensitivity analysis (process and parameter) 
# 
# AWalker (walkerap@ornl.gov) 
# March 2017
#
################################

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



# plot process SA indices
############################################

# read process SA index list
setwd(wdt)
# l1  <- readRDS(paste0(paste(runid_out,delta_var,'psa_list',sep='_'),'.RDS'))
df1   <- read.csv(paste(runid_out,delta_var,'psa.csv',sep='_'))
df1t  <- sens_table(df1)

# add environmetal variables to sensitivity output matrix
df1te <- add_scenario_values(df1t)
  
# write latex table, & create pdf (if call=TRUE)
write_Latex(df1te, paste(runid_out,'procSAtables',sep='_'), write_table_Latex, call=T )

# # plotting dataframe - make columnwise 
# dfpsi <- stack(as.data.frame(matpsi))
# dfpsi$par  <- dfpsi$ind
# dfpsi$ind  <- paste0('S',1:dim(matpsi)[1])
# dfpsi$type <- 'By Environment' 

# for less than 5 processes
# dfp     <- data.frame(si=as.numeric(df1), proc=rep(colnames(df1),each=length(l1[[1]])) ) 
#dfp     <- data.frame(si=as.numeric(df1), proc=rep(dimnames(df1)[[3]],each=length(l1[[1]])) ) 

# # legend labels
lnames <- list( elabs )

# plot
setwd(wdp)
pdf(paste(paste(runid_out,'PSAplots','radincs',sep='_'),'.pdf',sep=''), width=3.5, height=7)
# radar_plot(dfpsi, vnames=NULL, lnames=lnames, max_si=1.0 )
radar_plot(df1, var1=NULL, vnames=NULL, lnames=lnames, max_si=1.0 )
dev.off()



# parameter sensitivity plots 
############################################

# read info data
# - this is a bit of a fix so model names can be applied to data in l3, need to integrate this so that labels are passed through all processing stages
setwd(wdd)
fname <- paste('out',runid,'salt','AB','dataf',sep='_')
l1    <- readRDS(paste(fname,'.RDS',sep=''))

# read data
setwd(wdt)
# l3 <- readRDS(paste(runid_out,delta_var,'salt_list.RDS',sep='_'))
df1  <- read.csv(paste(runid_out,delta_var,'salt.csv',sep='_'))
df1t <- sens_table(df1)

# add environmetal variables to sensitivity output matrix
df1te <- add_scenario_values(df1t)

# write latex table, & create pdf (if call=TRUE)
write_Latex(df1te, paste(runid_out,'saltSAtables',sep='_'), write_table_Latex, call=T )

# # integrated over models and scenarios
# dfp1a      <- data.frame(values=l3$incmodelscenario[[sens]])
# dfp1a$ind  <- 'MS1'
# dfp1a$type <- 'Combined'
# dfp1a$par  <- row.names(dfp1a)
# dfms       <- dfp1a
# 
# # integrated over models and present against different environmental conditions
# dfp1  <- sapply(l3$incmodel, function(sl) sl[[sens]])
# colnames(dfp1) <- paste('S',1:dim(dfp1)[2], sep='')
# dfp1a <- stack(data.frame(dfp1))
# dfp1a$type <- 'By Environment' 
# dfp1a$par  <- row.names(dfp1)
# dfm        <- dfp1a
# 
# # integrated over scenarios and against each model combination
# dfp1  <- sapply(l3$incscenario, function(sl) sl[[sens]])
# colnames(dfp1) <- paste('M',1:dim(dfp1)[2], sep='')
# dfp1a <- stack(data.frame(dfp1))
# dfp1a$type <- 'By Model' 
# dfp1a$par  <- row.names(dfp1)
# dfs        <- dfp1a
# 
# # plotting dataframe
# dfsi    <- rbind(dfms,dfs,dfm)

# model labels
full_mod_names <- apply(l1$fnames,1,paste,collapse='.')
ss_names       <- if(!is.null(mod_names)) match(full_mod_names,names(mod_names)) else NULL
labs           <- if(!is.null(mod_names)) mod_names[ss_names] else full_mod_names
labs           <- if(length(labs)>12) paste0('M',1:length(labs)) else labs  

# legend labels
lnames <- list('weighted mean', labs, elabs )

# plot
setwd(wdp)
pdf(paste(paste(runid_out,'SAplots','radincs',sep='_'),'.pdf',sep=''), width=10, height=7)
radar_plot(dfsi, lnames=lnames, max_si=1.0 )
dev.off()

df1s <- subset(df1, scenario!=-1&model!=-1 )
radar_plot(df1s, var1='model', lnames=lnames, max_si=1.0 )




# broken out by scenario within each model combination
# model is on the outer list, scenario the inner list
# dfp1  <- lapply(l3, function(l) lapply(l$individual, function(sl) sapply(sl, function(ssl) ssl$sensitivity1) ))
# names(dfp1$absolute) <- paste('M',1:prod(n_procs),sep='')
# names(dfp1$response) <- paste('M',1:prod(n_procs),sep='')
# dfp1l <- lapply(dfp1, function(l) sapply(l, function(sl) sl) )
# dfp1a <- do.call('rbind',lapply(dfp1l,function(m) data.frame(stack(as.data.frame(m))) ))
# dfp1a$type <- output_label[1]
# dfp1a$type[which(grepl('res',row.names(dfp1a)))] <- output_label[2] 
# 
# # params
# dfp1a$par  <- row.names(dfp1$absolute[[1]])


# dfp1  <- lapply(l3$individual, function(sl) sapply(sl, function(ssl) ssl[[sens]]) )
# names(dfp1) <- paste('M',1:prod(n_procs),sep='')
# dfp1l <- vapply(dfp1, function(sl) sl, dfp1[[1]])
# 
# dfp1a <- stack(data.frame(dfp1))
# dfp1a$type <- 'incS' 
# dfp1a$par  <- row.names(dfp1[[1]])
# 
# 
# # env
# dfp1a$co2    <- NA
# dfp1a$PAR    <- NA
# dfp1a$co2[which(grepl('abs',row.names(dfp1a)))] <- rep(abs_co2,each=14) 
# dfp1a$PAR[which(grepl('abs',row.names(dfp1a)))] <- rep(abs_par,each=14) 
# dfp1a$co2[which(grepl('res',row.names(dfp1a)))] <- rep(res_co2,each=14) 
# dfp1a$PAR[which(grepl('res',row.names(dfp1a)))] <- rep(res_par,each=14) 
# 
# # plot
# pdf(paste(paste(runid_out,delta_var,'SAplots','radall',sep='_'),'.pdf',sep=''), width=1.5*10, height=2*6,pointsize=20)
# radar_plot(subset(dfp1a,type=='absolute'),var1='PAR',var2='co2')
# dev.off()



### END ###