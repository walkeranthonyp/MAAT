################################
#
# General functions for MAAT model 
# 
# AWalker August 2015
#
################################

library(xtable)


#####################################################
write_to_file <- function(df,ofile,app=F,type='csv') {
  if(type=='csv') write.table(format(df,width=12),paste(ofile,'.csv',sep=''),quote=F,row.names=F,col.names=!app,sep=',',append=app)  
  else if(type=='rds') saveRDS(df,paste(ofile,'RDS',sep='.'))
  else print(paste('No methods for output file format:',type))
} 



#####################################################
res_calc <- function(df1,df2,rcols,type='abs') {
  dfo <- df1
  if(type=='abs')      dfo[,rcols] <- df1[,rcols] - df2[,rcols]
  else if(type=='rel') dfo[,rcols] <- df1[,rcols] / df2[,rcols]
  else { print(paste('no method for type:',type)); stop }
  dfo
}


#####################################################
# write Latex doc
write_Latex <- function(obj,fname,func=write_list_Latex,call=F,...) {
  write("\\documentclass[10pt]{article}",fname)
  write("\\usepackage{booktabs}",fname,append=T)
  write("\\usepackage{placeins}",fname,append=T)
#   write("\\usepackage[section]{placeins}",fname,append=T)
  write("\\begin{document}",fname,append=T)

  func(obj,fname,...)
  
  write("\\end{document}",fname,append=T)  
  if(call) system(paste('pdflatex',fname)) 
}

write_table_Latex <- function(sub,fname,cnames=NULL,rnames=NULL,capname=NULL) {
  if(!is.null(cnames)) colnames(combtable)  <- cnames
  print.xtable.apw(sub,cap=capname,file=fname,append=T,include.rownames=!is.null(rnames))      
}

write_tablelist_Latex <- function(slist,fname,cnames=NULL,rnames=NULL,capname=NULL) {
  for(i in 1:length(slist)) {
    write_table_Latex(slist[[i]],fname,capname=capname[i])
  }
}

# write list recursively
write_list_Latex <- function(sublist,fname) {
  for(i in 1:length(sublist)) {
      if(typeof(sublist[[i]]) == "list") {
        write_list_Latex(sublist[[i]],fname)
      } else {
        if(class(sublist[[i]]) == "numeric") print.xtable.apw(t(as.matrix(sublist[[i]])),cap=names(sublist)[1],file=fname,append=T) else print.xtable.apw(sublist[[i]],cap=names(sublist)[1],file=fname,append=T) 
      } 

  }
  write("\\FloatBarrier",fname,append=T)
}

# write list recursively, combine lowest level list
write_list_Latex_comb <- function(sublist,fname,cnames=NULL,rnames=NULL,...) {
  if(typeof(sublist[[1]]) != "list") {
    write_list_Latex(sublist,fname)
    
  } else if(typeof(sublist[[1]][[1]]) == "list") {
    for(i in 1:length(sublist)) {
      write_list_Latex_comb(sublist[[i]],fname,cnames,rnames)
    }
    
  } else {
    for(li in 1:length(sublist[[1]])) {
      combtable <- do.call('rbind',lapply(sublist,function(l) l[[li]]))
      print(combtable)
      print(cnames)
      print(rnames)
      if(length(rnames)==dim(combtable)[1]) row.names(combtable) <- rnames
      if(length(cnames)==dim(combtable)[2]) colnames(combtable)  <- cnames
      print.xtable.apw(combtable,cap=names(sublist[[1]])[li],file=fname,append=T)      
    }
  } 
  write("\\FloatBarrier",fname,append=T)
}


print.xtable.apw <- function(x,cap=NULL,...){
  
  print(xtable(x,caption=cap,...),
        floating=T,
        caption.placement='top',
        hline.after=NULL,
        add.to.row=list(pos=list(-1,0, nrow(x)),
                        command=c(
                          '\\toprule\n',
                          '\\midrule\n',
                          '\\bottomrule\n')),
        table.placement='H',
        ...)
}



