################################
#
# Canopy structure system functions 
# 
# AWalker October 2018
#
################################



# Leaf demography  
###############################
f_canstrctsys_leafdem <- function(.) {

  # set correct canopy output
  .$canopy$cpars$output <- 'leaf_dem'    
  
  # initialise canopy_structure
  .$state$lai <- get(.$fnames$lai)(.)
  
  # pass LAI to canopy object
  .$canopy$fnames$lai <- 'f_lai_constant' 
  .$canopy$pars$lai   <- .$state$lai 

  # calculate fractions of leaf cohorts in upper and lower canopy
  .$state$upper_can_prop <- get(.$fnames$leafdem_upper)(.)
  .$state$lower_can_prop <- get(.$fnames$leafdem_lower)(.)

  # fix vcmax as a constant at the leaf level 
  .$canopy$leaf$fnames$vcmax  <- 'f_vcmax_constant'  

  # scale vcmax as decreasing exponential function through the canopy 
  .$canopy$fnames$vcmax0      <- 'f_vcmax0_constant'  
  .$canopy$fnames$scale_vcmax <- 'f_scale_vcmax_beerslaw'  

  # calculate vcmax for leaf cohorts
  .$state$leafdem_traits$canopy.vcmax0 <- get(.$fnames$leafdem_vcmax0)(.)

  # initialise layers
  #linc           <- .$state$lai / .$pars$layers
  #ca_calc_points <- seq(linc, .$state$lai, linc ) 
  #layers         <- .$pars$layers # this could be a function to dynamically specify the no. of layers 
  .$canopy$pars$layers <- .$pars$layers     # this could be a function to dynamically specify the no. of layers 
  .$init_vert(l=.$pars$layers)              # reallocating this memory is unnecessary in cases where layers is a fixed parameter. 
 
  # run canopy for each leaf cohort
  # create canopy parameter matrix
  lcohorts <- 3
  cohortmatrix <- vapply(.$state$leafdem_traits[c('canopy.vcmax0','canopy.vcmax0')], function(v) v, numeric(lcohorts) )
  #print(cohortmatrix)
  #print(class(cohortmatrix))
  #print(colnames(cohortmatrix))

  # run canopy - need to output layer matrix 
  #canopy_out <- vapply(1:lcohorts, .$run_canopy, matrix(0,.$pars$layers,length(.$canopy$output()), df=cohortmatrix )
  canopy_out <- vapply(1:lcohorts, .$run_canopy, matrix(0,.$pars$layers,length(.$canopy$state$vert$layer)), df=cohortmatrix )
  #canopy_out <- vapply(1:lcohorts, .$run_canopy, .$canopy$output(), df=cohortmatrix )

  # assign data to canopy_structure object data structure
  for(vname in colnames(canopy_out)) {
    .$state$vert$young[[vname]][]  <- canopy_out[,vname,1]
    .$state$vert$mature[[vname]][] <- canopy_out[,vname,2]
    .$state$vert$old[[vname]][]    <- canopy_out[,vname,3]
  }
  if(.$cpars$verbose) {
    print('Leaf cohorts:', quote=F )
    print(cohortmatrix)
    print(canopy_out)
  }


  # combine leaf cohorts 
  linc  <- .$state$lai / .$pars$layers
  upper <- floor(.$pars$lai_upper / linc) 
  for(vname in names(.$state$vert$layer)) 
    .$state$vert$layer[[vname]][1:upper] <- 
      .$state$vert$young[[vname]][1:upper]  * .$state$upper_can_prop[1] + 
      .$state$vert$mature[[vname]][1:upper] * .$state$upper_can_prop[2] +
      .$state$vert$old[[vname]][1:upper]    * .$state$upper_can_prop[3] 

  # is this the correct way to scale rs etc? --- no, best to convert to gs
  for(vname in names(.$state$vert$layer)) 
    .$state$vert$layer[[vname]][(upper+1):.$pars$layers] <- 
       .$state$vert$young[[vname]][(upper+1):.$pars$layers]  * .$state$lower_can_prop[1] + 
       .$state$vert$mature[[vname]][(upper+1):.$pars$layers] * .$state$lower_can_prop[2] +
       .$state$vert$old[[vname]][(upper+1):.$pars$layers]    * .$state$lower_can_prop[3] 

  # integrate canopy_structure layers
  # canopy_structure sum values
  for(vname in c('apar','A','gb','gs','gi','g','rd','Acg_lim','Ajg_lim','Apg_lim') ) {
    .$state$integrated[[vname]][] <- sum(.$state$vert$layer[[vname]]) * linc
  }
  # canopy_structure mean values 
  .$state$integrated$cb[] <- sum(.$state$vert$layer$cb) / .$state$lai * linc
  .$state$integrated$cc[] <- sum(.$state$vert$layer$cc) / .$state$lai * linc
  .$state$integrated$ci[] <- sum(.$state$vert$ci)       / .$state$lai * linc
  
}



### END ###
