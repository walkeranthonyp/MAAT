################################
#
# Script to generate default and option XMLs from a newly created MAAT model object 
# call this from the command line inside the MAAT directory with RScript   
#
# AWalker November 2017
#
################################

mod_obj <- NULL

# parse command line arguments   
# - any one of the above objects can be specified as a command line argument using the syntax:
# - Rscript <nameofthisscript> "<object1> <- <value1>" 
# - e.g. Rscript create_default_XMLs.R "mod_obj <- 'template'" 

if(length(commandArgs(T))>=1) {
  for( ca in 1:length(commandArgs(T)) ) {
    eval(parse(text=commandArgs(T)[ca]))
  }
}

if(is.null(mod_obj)) {
  stop('mod_obj needs to be specified as a command line argument')
} 

# load function to creat XMLs from lists
setwd('./src/')
source('general_functions.R')

# open the newly created model object
setwd(mod_obj)
mo <- paste0(mod_obj,'_object')
source(paste0(mo,'.R'))

# create lists to convert to XML
l1 <- list(list(fnames = get(mo)[['fnames']],
                pars   = get(mo)[['pars']],
                env    = get(mo)[['env']]
              ))
names(l1) <- mod_obj

l2 <- list(list(env    = get(mo)[['env']]))
names(l2) <- mod_obj
l2[[1]][[1]][] <- 'column name of variable in metdata file'

# convert list to XMLs
listtoXML(paste(mod_obj,'default.xml',sep='_'), 'default',  sublist=l1 )
listtoXML(paste(mod_obj,'options.xml',sep='_'), 'options',  sublist=l1 )
setwd('init_files')
listtoXML(paste(mod_obj,'user_static.xml',sep='_'),  'static',               sublist=l1 )
listtoXML(paste(mod_obj,'user_dynamic.xml',sep='_'), 'dynamic',              sublist=l1 )
listtoXML(paste(mod_obj,'user_met.xml',sep='_'),     'met_data_translator',  sublist=l2 )

# open options XML to add labels and non-default options
setwd('..')
file.edit(paste(mod_obj,'options.xml',sep='_'))
