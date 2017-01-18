################################
#
# General functions for MAAT model 
# 
# AWalker August 2015
#
################################

library(xtable)
library(XML)


#####################################################
write_to_file <- function(df,ofile,app=F,type='csv') {
  if(type=='csv') write.table(format(df,width=12),paste(ofile,'.csv',sep=''),quote=F,row.names=F,col.names=!app,sep=',',append=app)  
  else if(type=='rds') saveRDS(df,paste(ofile,'RDS',sep='.'))
  else print(paste('No methods for output file format:',type))
} 

readXML <- function(input.file=NULL) {
  settings.xml <- NULL
  
  # Parse input xml
  if (!is.null(input.file) && file.exists(input.file)) {
    settings.xml <- xmlParse(input.file)  
    # convert the xml to a list
    settings.list <- xmlToList(settings.xml)
    
  } else {
    print("***** WARNING: readXML called but with no input file defined *****")
  }
  
  # make sure something was loaded
  if (is.null(settings.xml)) {
    #log.error("Did not find any settings file to load.")
    stop(paste('readXML called, input file empty'))
  }
  
  # Remove comment or NULL fields
  settings.list <- settings.list[settings.list !="NULL" ]
  
  # Evaluate and return XML as a list
  return(evalXMLlist(settings.list))
  
}

# run over all elements of a list, returning a list of the same structure
# evaluate character strings to r objects
evalXMLlist <- function(sublist) {
  for(i in 1:length(sublist)) {
    sublist[[i]] <- 
      if(typeof(sublist[[i]]) == "list")  evalXMLlist(sublist[[i]]) else eval(parse(text=sublist[[i]]))
  }
  sublist
}

# save a list to an xml file
listtoXML <- function(fname,name,...) {
  # expects the argument - sublist - in ...
  
  root <- newXMLNode(name)
  
  rec <- function(node, sublist) {
    for(i in 1:length(sublist)) {
      child <- newXMLNode(names(sublist)[i], parent=node)
      
      if (typeof(sublist[[i]]) == "list") rec(child, sublist[[i]]) else xmlValue(child) <- sublist[[i]]
    }
  }  
  
  rec(root,...)
  saveXML(root,fname)
}

# overwrite values in mainlist with values in sublist
# - currently expects list to be of identical structure
# fuselists <- function(mainlist,sublist) {
#   for(i in 1:length(sublist)) {
#     mlsub <- which(names(mainlist)==names(sublist)[i])
#     if(length(mlsub)!=1) {
#       stop(paste("names mismatch when fusing initialisation lists:",mainlist,sublist))      
#     }
#     mainlist[[mlsub]] <- 
#       if(typeof(sublist[[i]]) == "list")  fuselists(mainlist[[mlsub]],sublist[[i]]) else sublist[[i]]
#   }
#   mainlist
# }

fuselists <- function(mainlist,sublist) {
  l <- length(sublist)
  n <- 0
  for(i in 1:l) {
    mlsub <- which(names(mainlist)==names(sublist)[i])
    if(length(mlsub)==1) {
      n <- n + 1
      mainlist[[mlsub]] <- 
        if(typeof(sublist[[i]]) == "list")  fuselists(mainlist[[mlsub]],sublist[[i]]) else sublist[[i]]
    }
    # print(c(n,l))
    # print(names(sublist)[i])
  }
  if(n!=l) stop(paste("\n names mismatch when fusing initialisation lists:",names(mainlist),names(sublist)[i]))      
  mainlist
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



