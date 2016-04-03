################################
#
# General functions for MAAT model 
# 
# AWalker August 2015
#
################################

library(XML)

write_to_file <- function(df,ofile,app=F,type='csv') {
  if(type=='csv') write.table(format(df,width=12),paste(ofile,'.csv',sep=''),quote=F,row.names=F,col.names=!app,sep=',',append=app)  
} 

readXML <- function(input.file=NULL) {
  settings.xml <- NULL
  
  ### Parse input settings file
  if (!is.null(input.file) && file.exists(input.file)) {
    settings.xml <- xmlParse(input.file)  
    # convert the xml to a list
    settings.list <- xmlToList(settings.xml)
    
  } else {
    print("***** WARNING: no settings file defined *****")
  }
  
  # make sure something was loaded
  if (is.null(settings.xml)) {
    #log.error("Did not find any settings file to load.")
    stop("Did not find any settings file to load :: STOP")
  }
  
  ### Remove comment or NULL fields
  settings.list <- settings.list[settings.list !="NULL" ]
  
  # Return settings file as a list
  return(settings.list)
  
}

# overwrite values in mainlist with values in sublist
fuselists <- function(mainlist,sublist) {
  for(i in 1:length(sublist)) {
    mlsub <- which(names(mainlist)==names(sublist)[i])
    if(length(mlsub)!=1) {
      stop("names mismatch when fusing initialisation lists :: STOP")      
    }
    mainlist[[mlsub]] <- 
      if(typeof(sublist[[i]]) == "list")  fuselists(mainlist[[mlsub]],sublist[[i]]) else sublist[[i]]
  }
  mainlist
}

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








