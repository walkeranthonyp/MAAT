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

  # initialise canopy_structure
  .$state$lai <- get(.$fnames$lai)(.)
  
  # pass LAI to canopy object
  .$canopy$fnames$lai <- 'f_lai_constant' 
  .$canopy$pars$lai   <- .$state$lai 

  # calculate fractions of leaf cohorts in upper and lower canopy
  .$state$upper_can_prop <- get(.$fnames$leafdem_upper)(.)
  .$state$lower_can_prop <- get(.$fnames$leafdem_lower)(.)

  # fix vcmax as a constant at the leaf level ... but can canopy handle passing vcmax as a constant? not yet, 
  .$canopy$leaf$fnames$vcmax <- 'f_vcmax_constant'  

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
  cohortmatrix <- vapply(.$state$leafdem_traits[c('canopy.vcmax0')], function(v) v, numeric(lcohorts) )
  # run canopy
  canopy_out <- vapply(1:lcohorts, .$run_canopy, .$canopy$output(), df=cohortmatrix )
  # assign data to canopy_structure object data structure
  for(vname in row.names(canopy_out)) {
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
  upper <- floor(.$pars$upper / linc) 
  for(vname in names(.$state$vert$layer)) .$state$vert$layer[[vname]][1:upper] <- .$state$vert$young[[vname]][1:upper]  * .$upper_can_prop[1] + 
                                                                                  .$state$vert$mature[[vname]][1:upper] * .$upper_can_prop[2] +
                                                                                  .$state$vert$old[[vname]][1:upper]    * .$upper_can_prop[3] 

  # is this the correct way to scale rs etc? --- no, best to convert to gs
  for(vname in names(.$state$vert$layer)) .$state$vert$layer[[vname]][(upper+1):layers] <- .$state$vert$young[[vname]][(upper+1):layers]  * .$upper_can_prop[1] + 
                                                                                           .$state$vert$mature[[vname]][(upper+1):layers] * .$upper_can_prop[2] +
                                                                                           .$state$vert$old[[vname]][(upper+1):layers]    * .$upper_can_prop[3] 

  # integrate canopy_structure layers
  # canopy_structure sum values
  .$state$integrated$A[]              <- sum(.$state$vert$layer$A) * linc
  .$state$integrated$respiration[]    <- sum(.$state$vert$layer$respiration) * linc
#  .$state$integrated$Acg_lim[]        <- sum(.$state$vert$layer$A * (.$state$lim==2)) * linc 
#  .$state$integrated$Ajg_lim[]        <- sum(.$state$vert$layer$A * (.$state$lim==3)) * linc
#  .$state$integrated$Apg_lim[]        <- sum(.$state$vert$layer$A * (.$state$lim==7)) * linc
#  .$state$integrated$layers_Acg_lim[] <- sum(.$state$vert$layer$lim==2)
#  .$state$integrated$layers_Ajg_lim[] <- sum(.$state$vert$layer$lim==3)
#  .$state$integrated$layers_Apg_lim[] <- sum(.$state$vert$layer$lim==7)
  .$state$integrated$rb[]             <- 1 / sum(1/.$state$vert$layer$rb * linc )
  .$state$integrated$rs[]             <- 1 / sum(1/.$state$vert$layer$rs * linc )
  .$state$integrated$ri[]             <- 1 / sum(1/.$state$vert$layer$ri * linc )
  # canopy_structure mean values -- where is cb?
  .$state$integrated$cc[]             <- sum(.$state$vert$layer$cc) / .$state$lai * linc
  .$state$integrated$ci[]             <- sum(.$state$vert$ci) / .$state$lai * linc
  
}



### END ###
