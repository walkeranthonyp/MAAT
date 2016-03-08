################################
#
# General functions for MAAT model 
# 
# AWalker August 2015
#
################################

write_to_file <- function(df,ofile,app=F,type='csv') {
  if(type=='csv') write.table(format(df,width=12),paste(ofile,'.csv',sep=''),quote=F,row.names=F,col.names=!app,sep=',',append=app)  
  else if(type=='rds') saveRDS(df,paste(ofile,'RDS',sep='.'))
  else print(paste('No methods for output file format:',type))
} 










