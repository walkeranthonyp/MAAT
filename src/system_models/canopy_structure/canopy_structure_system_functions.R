################################
#
# Canopy structure system functions 
# 
# AWalker October 2018
#
################################



# Leaf demography  
###############################
f_sys_leafdem <- function(.) {

  # set correct canopy output
  .super$canopy$cpars$output[] <- 'leaf_dem'    
  
  # initialise canopy_structure
  .super$state$lai[] <- .$lai()
  
  # pass LAI to canopy object
  .super$canopy$fnames$lai[] <- 'f_lai_constant' 
  .super$canopy$pars$lai[]   <- .super$state$lai 

  # calculate fractions of leaf cohorts in upper and lower canopy
  .super$state$upper_can_prop <- .$leafdem_upper()
  .super$state$lower_can_prop <- .$leafdem_lower()

  # fix vcmax as a constant at the leaf level 
  .super$canopy$leaf$fnames$vcmax  <- 'f_vcmax_constant'  

  # scale vcmax as decreasing exponential function through the canopy 
  .super$canopy$fnames$vcmax0      <- 'f_vcmax0_constant'  
  .super$canopy$fnames$scale_vcmax <- 'f_scale_vcmax_beerslaw'  

  # calculate vcmax for leaf cohorts
  .super$state$leafdem_traits$canopy.vcmax0 <- .$leafdem_vcmax0()

  # initialise layers
  .super$canopy$pars$layers <- .super$pars$layers     # this could be a function to dynamically specify the no. of layers 
  .super$init_vert(.=.super, l=.super$pars$layers)    # reallocating this memory is unnecessary in cases where layers is a fixed parameter. 
 
  # run canopy for each leaf cohort
  # create canopy parameter matrix
  lcohorts <- 3
  cohortmatrix <- vapply(.super$state$leafdem_traits[c('canopy.vcmax0','canopy.vcmax0')], function(v) v, numeric(lcohorts) )
  #print(cohortmatrix)
  #print(class(cohortmatrix))
  #print(colnames(cohortmatrix))

  # run canopy - need to output layer matrix 
  canopy_out <- vapply(1:lcohorts, .$run_canopy, matrix(0,.super$pars$layers,length(.super$canopy$state$vert$layer)), df=cohortmatrix )

  # assign data to canopy_structure object data structure
  for(vname in colnames(canopy_out)) {
    .super$state$vert$young[[vname]][]  <- canopy_out[,vname,1]
    .super$state$vert$mature[[vname]][] <- canopy_out[,vname,2]
    .super$state$vert$old[[vname]][]    <- canopy_out[,vname,3]
  }
  if(.super$cpars$verbose) {
    print('Leaf cohorts:', quote=F )
    print(cohortmatrix)
    print(canopy_out)
  }


  # combine leaf cohorts 
  linc  <- .super$state$lai / .super$pars$layers
  upper <- floor(.super$pars$lai_upper / linc) 
  for(vname in names(.super$state$vert$layer)) 
    .super$state$vert$layer[[vname]][1:upper] <- 
      .super$state$vert$young[[vname]][1:upper]  * .super$state$upper_can_prop[1] + 
      .super$state$vert$mature[[vname]][1:upper] * .super$state$upper_can_prop[2] +
      .super$state$vert$old[[vname]][1:upper]    * .super$state$upper_can_prop[3] 

  # is this the correct way to scale rs etc? --- no, best to convert to gs
  for(vname in names(.super$state$vert$layer)) 
    .super$state$vert$layer[[vname]][(upper+1):.super$pars$layers] <- 
       .super$state$vert$young[[vname]][(upper+1):.super$pars$layers]  * .super$state$lower_can_prop[1] + 
       .super$state$vert$mature[[vname]][(upper+1):.super$pars$layers] * .super$state$lower_can_prop[2] +
       .super$state$vert$old[[vname]][(upper+1):.super$pars$layers]    * .super$state$lower_can_prop[3] 

  # integrate canopy_structure layers
  # canopy_structure sum values
  for(vname in c('apar','A','gb','gs','gi','g','rd','Acg_lim','Ajg_lim','Apg_lim') ) {
    .super$state$integrated[[vname]][] <- sum(.super$state$vert$layer[[vname]]) * linc
  }
  # canopy_structure mean values 
  .super$state$integrated$cb[] <- sum(.super$state$vert$layer$cb) / .super$state$lai * linc
  .super$state$integrated$cc[] <- sum(.super$state$vert$layer$cc) / .super$state$lai * linc
  .super$state$integrated$ci[] <- sum(.super$state$vert$ci)       / .super$state$lai * linc
  
}



### END ###
