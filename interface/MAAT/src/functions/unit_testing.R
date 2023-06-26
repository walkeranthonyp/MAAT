###########################
#
# Unit testing of MAAT general functions
#
# AWalker
# Nov 2018
#
###########################



rm(list=ls())
library(lattice)
source('general_functions.R')

l1 <- readXML('test.xml')
l2 <- rapply(l1, function(x) if(is.na(x)) NULL else x, how='replace' )
class(unlist(l2))

ln <- list(NULL)
ln2 <- rapply(ln, function(x) if(is.na(x)) NULL else x, how='replace' )
unlist(ln2)
unlist(ln2$fnames$f)
