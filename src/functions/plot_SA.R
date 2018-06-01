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
l1      <- readRDS(paste0(paste(runid_out,delta_var,'psa_list',sep='_'),'.RDS'))

# extract process SA index
ex_sens <- function(l) vapply(l, function(l) l$sensitivity, numeric(1) )
df1     <- vapply(l1, ex_sens, array(0, dim=c(1, length(l1[[1]]))) )
if(class(df1)=='numeric') df1 <- array(df1,c(1,length(df1)),dimnames=list(NULL,names(df1)))
matpsi <- apply(df1,3,function(m) m )

# extract process SA summary stats
ex_mean <- function(l) vapply(l, function(l) l$Tvar, numeric(3) )
df2     <- vapply(l1, ex_mean, array(0, dim=c(length(l1[[1]][[1]]$Tvar), length(l1[[1]]))) )
apply(df2,1:2,function(m) m )

# generate sensitivity output matrix
if(!is.null(evar1)&!is.null(evar2)) {
  # t1                <- cbind( evar1, evar2, t(round(df2[,,1],2)), round(df1,2) )
  t1                <- cbind( evar1, evar2, t(round(df2[,,1],2)), round( matpsi , 2 ) )
  colnames(t1)[1:2] <- c(evar1name, evar2name )
} else if(!is.null(evar1)) {
  # t1                <- cbind( evar1, t(round(df2[,,1],2)), round(df1,2) )
  t1                <- cbind( evar1, t(round(df2[,,1],2)), round( matpsi , 2 ) )
  colnames(t1)[1]   <- evar1name
} else {
  t1                <- cbind( t(round(df2[,,1],2)), round(df1,2) )
}

# write latex table, & create pdf (if call=TRUE)
write_Latex(t1, paste(runid_out,'procSAtables',sep='_'), write_table_Latex, call=T )

# plotting dataframe - make columnwise 
dfpsi <- stack(as.data.frame(matpsi))
dfpsi$par  <- dfpsi$ind
dfpsi$ind  <- paste0('S',1:dim(matpsi)[1])
dfpsi$type <- 'By Environment' 

# for less than 5 processes
dfp     <- data.frame(si=as.numeric(df1), proc=rep(colnames(df1),each=length(l1[[1]])) ) 

# # legend labels
lnames <- list( elabs )

# plot
setwd(wdp)
pdf(paste(paste(runid_out,'PSAplots','radincs',sep='_'),'.pdf',sep=''), width=3.5, height=7)
radar_plot(dfpsi, vnames=NULL, lnames=lnames, max_si=1.0 )
dev.off()



# parameter sensitivity plots 
############################################

# read info data
# - this is a bit of a fix so model names can be applied to data in l3, need to integrate this so that labels are passed through all processing stages
setwd(wdd)
fname <- paste(runid,'out','salt','AB','dataf',sep='_')
l1    <- readRDS(paste(fname,'.RDS',sep=''))

# read data
setwd(wdt)
l3 <- readRDS(paste(runid_out,delta_var,'salt_list.RDS',sep='_'))


# generate Latex tables of metrics 
# - do this from a called script, that way the script can be used without calling all the data processing routines again

# integrated over models and scenarios
dfp1a      <- data.frame(values=l3$incmodelscenario[[sens]])
dfp1a$ind  <- 'MS1'
dfp1a$type <- 'Combined'
dfp1a$par  <- row.names(dfp1a)
dfms       <- dfp1a

# integrated over models and present against different environmental conditions
dfp1  <- sapply(l3$incmodel, function(sl) sl[[sens]])
colnames(dfp1) <- paste('S',1:dim(dfp1)[2], sep='')
dfp1a <- stack(data.frame(dfp1))
dfp1a$type <- 'By Environment' 
dfp1a$par  <- row.names(dfp1)
dfm        <- dfp1a

# integrated over scenarios and against each model combination
dfp1  <- sapply(l3$incscenario, function(sl) sl[[sens]])
colnames(dfp1) <- paste('M',1:dim(dfp1)[2], sep='')
dfp1a <- stack(data.frame(dfp1))
dfp1a$type <- 'By Model' 
dfp1a$par  <- row.names(dfp1)
dfs        <- dfp1a

# plotting dataframe
dfsi    <- rbind(dfms,dfs,dfm)

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