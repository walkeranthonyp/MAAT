################################
#
# Script to generate default and option XMLs from a newly created MAAT model object 
# call this from the command line inside the run_scripts directory with RScript   
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

# load function to create XMLs from lists
setwd('~/Documents/MAAT/src/functions/')
source('general_functions.R')

# open the newly created model object
setwd(paste('..','system_models',mod_obj,sep='/'))
mo <- paste0(mod_obj,'_object'); source(paste0(mo,'.R'))

# create lists to convert to XML
l1 <- list(fnames = list(get(mo)[['fnames']]),
           pars   = list(get(mo)[['pars']]),
           env    = list(get(mo)[['env']])
           )
names(l1$fnames) <- mod_obj
names(l1$pars)   <- mod_obj
names(l1$env)    <- mod_obj

# add double quote to character strings in fnames - this allows proper parsing when the XMLs are read by MAAT
#l1[[mod_obj]][['fnames']] <- lapply(l1[[mod_obj]][['fnames']],  function(c1) paste0("'",c1,"'")

# create met data list
l2 <- list(env = list(get(mo)[['env']]))
names(l2$env) <- mod_obj
l2[[1]][[1]][] <- 'column name of variable in metdata file'

# create eval data list
l3 <- list(state = list(as.list(unlist(get(mo)[['state']]))))
names(l3$state) <- mod_obj
l3[[1]][[1]][] <- 'column name of variable in evaldata file'

# read all options
mf <- paste0(mod_obj,'_system_functions.R'); source(mf)
mf <- paste0(mod_obj,'_functions.R');        source(mf)
l1opt     <- l1['fnames']
l1names   <- l1[['fnames']][[mod_obj]]
l1names[] <- names(l1names)
l1opt[['fnames']][[mod_obj]] <- 
  lapply(l1names, function(c1) { 
    print('', quote=F )
    print(paste('Searching for representations of process:',c1), quote=F ) 
    c2 <- ls(pos=1, pattern=paste0('f_',c1))
    print('Found:', quote=F ) 
    print(c2, quote=F ) 
    out <- do.call('paste', lapply(as.list(c2),  function(c3) paste0("'",c3,"'",',') ) )
    out <- substr(out, 1, nchar(out)-1)
    paste('c(',out,')')
  })

# replace all entires in l1 with NA for init static and dynamic
l1n <- rapply(l1, function(x) NA, how='replace' )  

# convert list to XMLs
print('', quote=F )
listtoXML(paste(mod_obj,'default.xml',sep='_'), 'default',  sublist=l1 )
listtoXML(paste(mod_obj,'options.xml',sep='_'), 'options',  sublist=l1opt )
setwd('init_files')
listtoXML(paste(mod_obj,'user_static.xml',sep='_'),  'static',               sublist=l1n )
listtoXML(paste(mod_obj,'user_dynamic.xml',sep='_'), 'dynamic',              sublist=l1n )
listtoXML(paste(mod_obj,'user_met.xml',sep='_'),     'met_data_translator',  sublist=l2 )
listtoXML(paste(mod_obj,'user_eval.xml',sep='_'),    'eval_data_translator', sublist=l3 )

# open options XML to add labels by hand 
setwd('..')
#file.edit(paste(mod_obj,'options.xml',sep='_'))



### END ###
