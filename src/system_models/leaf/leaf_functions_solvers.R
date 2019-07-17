###############################
#
# Leaf level physiology functions - photosynthesis solvers and associated functions
# 
# AWalker July 2019
# 13 C discrimination: Soumaya Belmechri, Juergen Knauer March 2019 
#
################################



### SOLVERS & RESIDUAL FUNCTIONS
################################

# Solver to find root of .$solver_func
#f_solver_brent <- function(., interval=c(.$pars$solver_min,.$pars$solver_max) ) {
f_R_Brent_solver <- function(.) {
  if(.super$cpars$verbose_loop) print(.super$env) 

  #.$puniroot(.$residual_func, interval=interval, extendInt='downX' )$root
  .super$solver_out <- .$puniroot(.$solver_func, interval=c(-0.002765326,50.1234), extendInt='downX' )
  .super$solver_out$root
}


#f_solver_brent_diag <- function(., interval=c(.$pars$solver_min,.$pars$solver_max) ) {
f_R_Brent_solver_diag <- function(.) {
  if(.super$cpars$verbose_loop) print(.super$env) 

  #.super$solver_out        <- .$puniroot(.$solver_func, interval=interval, extendInt='downX' )
  .super$solver_out        <- .$puniroot(.$solver_func, interval=c(-0.002765326,50.1234), extendInt='downX' )
  .super$state$iter[]      <- .super$solver_out$iter
  .super$state$estimprec[] <- .super$solver_out$estim.prec
  .super$solver_out$root
}

f_R_Brent_solver_diag_brackets <- function(., brackets ) {
  if(.super$cpars$verbose_loop) print(.super$env) 

  .super$solver_out        <- .$puniroot(.$solver_func, interval=brackets )
  #.super$state$iter[]      <- .super$solver_out$iter
  #.super$state$estimprec[] <- .super$solver_out$estim.prec
  .super$solver_out$root
}


# Residual function for solver to calculate assimilation
#f_residual_leaf_Ar <- function(., A ) {
f_A_r_leaf <- function(., A ) {
  # combines A, rs, ri, ci & cc eqs to a single f(A), 
  # combines all rate limiting processes
  #  -- for use with uniroot solver
  #  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation

  # calculate cc from ca, rb, rs, and ri
  # total resistance of a set of resistors in series is simply their sum 
  # assumes boundary layer and stomatal resistance terms are in h2o units
  # assumes mesophyll resistance is in co2 units
  .super$state$cc[] <- .$gas_diff( A , r=( 1.4*.super$state_pars$rb + 1.6*.$rs(A=A,c=.$gas_diff(A)) + .super$state_pars$ri ) )
  
  #print(c(A,.$assimilation()))
  
  # calculate residual of net A
  .$assimilation() - A
} 


# same as above function but with no stomatal resistance 
#f_residual_leaf_Ar_no_rs <- function(., A ) {
f_A_r_leaf_noRs <- function(., A ) {
  
  # calculate cc from ca, rb, and ri
  # assumes boundary layer resistance is in h2o units
  # assumes mesophyll resistance is in co2 units
  .super$state$cc <- .$gas_diff( A , r=( 1.4*.super$state_pars$rb + .super$state_pars$ri ) )
  
  # calculate residual of net A
  .$assimilation(.) - A
}


# semi-analytical solver
#f_solver_semiana_leaf <- function(.) {
f_A_r_leaf_semiana <- function(.) {
  # - finds the analytical solution assuming rb, ri, and g0 are zero to use as first guess (a1)
  # - or finds the analytical solution assuming rb and ri are zero to use as first guess (a1)
  # - calculate residual function value for a1 (fa1) 
  # - make a second guess (a2) from the sum of a1 + fa1   
  # - calculate residual function value for a2 (fa2) 
  # - make a third guess (a3) from the sum of a2 + fa2   
  # - calculate residual function value for a3 (fa3) 
  # - if selected, return third guess
  # - check fa1, fa2, and fa3 span the root
  # - if selected, fit a quadratic through the three sets of co-ordinates to find the root 
  # - if selected, use Brent to find the root using brackets 
 
  .super$state$aguess_flag <- 0 
 
  # find the analytical solution assuming rb and ri are zero to use as first guess (a1)
  .super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- .$analytical()

  # if f(e) < 0.1 then -999 is returned
  # this trap recalculates first guess using quadratic analytical solution
  if(.super$state$aguess[1] == -999) {
    .super$state$aguess_flag <- -999
    .super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- .$analytical_quad()
  } 

  # calculate residual of first guess
  .super$state$faguess[1] <- .$full_soln(.super$state$aguess[1]) 

  # if first guess is below zero return negative flag and leafsys routine will calculate A assuming gs = g0
  if(.super$state$aguess[1] < 0) { 
    .super$state$aguess_flag <- -99
    return(-99)
  }

  # if first residual is small return first guess
  else if(abs(.super$state$faguess[1]) < 1e-6) return(.super$state$aguess[1])
  else {

    # second guess, sum of residual and first guess 
    .super$state$aguess[2]  <- .super$state$faguess[1] + .super$state$aguess[1]  
    
    # if second guess or cc of second guess are below zero then return negative flag and leafsys routine will calculate A assuming gs = g0
    if(.super$state$aguess[2]<0 | .super$state$cc<0 ) {
      .super$state$aguess_flag <- -9
      return(-9)
  
    } else {
      # calculate residual of second guess
      .super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 
  
      # third guess
      .super$state$aguess[3]  <- .super$state$faguess[2] + .super$state$aguess[2]  
    
      # decide on final solver to use
      if(.super$fnames$semiana=='f_semiana_three') { 
    
        # return third guess
        .$semiana()
        
      } else { 
        
        .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
        
        # check that fa2 and fa1 span 0, if not reguess for fa3 until it and fa1 span zero and then take the mean a1 and a3 to find a2
        if( min(.super$state$faguess)*max(.super$state$faguess) >= 0 ) {
  
          # find guesses that span the root 
          .$iterate_guesses()
  
          if( prod(.super$state$faguess[1:2]) >= 0 ) {
            # reporting for development
            print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
            print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))
            print(unlist(.super$fnames)); print(unlist(.super$pars)); print(unlist(.super$state_pars)); print(unlist(.super$env)) 
            print(.$residual(centrala=.super$state$aguess[1]))
            stop(paste('Solver error'))
          }
        } 
      
        # call semi-ana solver final method 
        .$semiana()
       
      } 
    } 
  }
}


# return the third guess  
f_semiana_three <- function(.) .super$state$aguess[3]
  

# find the root using the Brent solver with brackets that span the root of the residual function
f_semiana_brent <- function(.) {

  brackets <- c( .super$state$aguess[which.min(.super$state$faguess)], .super$state$aguess[which.max(.super$state$faguess)] )
  # there could be a refinment step in here to select the narrowest bracket
  #.$solver_brackets <- f_R_Brent_solver_diag_brackets ###
  #f_R_Brent_solver_diag_brackets(., brackets )
  .$solver_brackets(brackets)
}  


# fit a quadratic through the three sets of co-ordinates a la Lomas (one step of Muller's method) 
f_semiana_quad <- function(.) {
  # Muller, David E., "A Method for Solving Algebraic Equations Using an Automatic Computer," Mathematical Tables and Other Aids to Computation, 10 (1956)    
  #bx <- ((a1+a0)*(fa2-fa1)/(a2-a1)-(a2+a1)*(fa1-fa0)/(a1-a0))/(a0-a2)
  #ax <- (fa2-fa1-bx*(a2-a1))/(a2**2-a1**2)
  #cx <- fa1-bx*a1-ax*a1**2
      
  # calculate residual of third guess
  .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
   
  bx <- ((.super$state$aguess[2]+.super$state$aguess[1])*(.super$state$faguess[3]-.super$state$faguess[2]) / 
  (.super$state$aguess[3]-.super$state$aguess[2])-(.super$state$aguess[3]+.super$state$aguess[2])*(.super$state$faguess[2]-.super$state$faguess[1])/(.super$state$aguess[2]-.super$state$aguess[1])) / 
  (.super$state$aguess[1]-.super$state$aguess[3])
  ax <- (.super$state$faguess[3]-.super$state$faguess[2]-bx*(.super$state$aguess[3]-.super$state$aguess[2]))/(.super$state$aguess[3]**2-.super$state$aguess[2]**2)
  cx <- .super$state$faguess[2]-bx*.super$state$aguess[2]-ax*.super$state$aguess[2]**2
       
  # find the root of the quadratic that lies between guesses 
  .super$state$assim[] <- quad_sol(ax,bx,cx,'both')
  ss <- which(.super$state$assim>min(.super$state$aguess[1:3])&.super$state$assim<max(.super$state$aguess[1:3]))
       
  # catch potential errors
  if(length(ss)==0)      stop('no solution,',.super$state$assim,'; within initial 3 guesses,',.super$state$aguess[1:3],'fguesses,',.super$state$faguess[1:3] ) 
  else if(length(ss)==2) stop('both solutions,',.super$state$assim,'; within initial 3 guesses,',.super$state$aguess[1:3],'fguesses,',.super$state$faguess[1:3] ) 
  else return(.super$state$assim[ss])
} 

  
#   
f_iterate_guesses <- function(.) {

  print('')
  print('fa1 and fa2 do not span zero')
  print(.super$state_pars); print(.super$env) 
  
  print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
  print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))
 
#  # recalculate third guess and residual of third guess
#  .super$state$aguess[3]  <- .super$state$faguess[2] + .super$state$aguess[2]  
#  .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
#  print('')
#  print('recalculate a3 as mean of a1 and a2 ')
#  print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
#  print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))

  # home in on guesses that span the root
  if(.super$state$faguess[1] > 0 ) { 
    while(.super$state$faguess[3]>0) {
      .super$state$aguess[3]  <- .super$state$aguess[3] + 0.01  
      .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
    }
  } else {
    while(.super$state$faguess[3] < 0) {
      .super$state$aguess[3]  <- .super$state$aguess[3] - 0.01  
      .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
    }
  }

  .super$state$aguess[2]  <- mean(.super$state$aguess[c(1,3)])
  .super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 

  # reporting for development
  #print(f_residual(.,centrala=.super$state$aguess[1]))
  print(.$residual(centrala=.super$state$aguess[1]))
 
  # over-cautious safety - with while condition above should be no way that this condition is true
  if(prod(.super$state$faguess[c(1,3)]) > 0) stop('failed')
} 


# calculate residual function around a central guess of A 
f_residual <- function(., centrala, range=2, inc=0.1 ) {
  a     <- seq(centrala-range, centrala+range, inc )
  resid <- numeric(length(a))

  for( i in 1:length(a) ) {
    #resid[i] <- .$full_soln(A=a[i])
    resid[i] <- .$solver_func(A=a[i])
  }

  cbind(a,resid)
}


# write resdidual for numerical solution
f_write_residual <- function(., fname='resid.csv', ... ) { 
  resid <- .$residual(...)
  write.csv(data.frame(resid), fname, append=T, row.names=F, quote=F )
}


## Semi-analytical solution 
#f_A_r_leaf_semiana_quad <- function(.) {
#  # - finds the analytical solution assuming rb and ri are zero to use as first guess (a1)
#  # - make a second guess (a2) by calculating Cc using a1 and actual values of rb and ri   
#  # - make a third guess (a3) by taking the mean of a1 and a2  
#  # - calculate residual function value for these three guesses (fa1, fa2, fa3) 
#  # - check fa1 and fa2 span the root
#  # - if not, sequentially add 0.01 to the third guess (a3) until fa1 and fa3 span the root
#  # - fit a quadratic through the three sets of co-ordinates to find the root 
# 
#  # find the analytical solution assuming rb and ri are zero to use as first guess (a1)
#  #.super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- f_A_r_leaf_analytical_quad(.)
#  #.super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- f_A_r_leaf_analytical(.)
#  .super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- .$analytical()
#  #.super$state$faguess[1] <- f_A_r_leaf(., .super$state$aguess[1] ) 
#  .super$state$faguess[1] <- .$full_soln(.super$state$aguess[1]) 
#
#  if(.super$state$aguess[1] == -999) {
#    #afinal <- f_A_r0_leaf_analytical_quad(.)      
#    afinal <- .$analytical_quad()   
#    ### can call residual function in here   
#    a <- seq(afinal-2,afinal+2,0.05 )
#    resid <- numeric(length(a))
#    for( i in 1:length(a) ) {
#    resid[i] <- f_A_r_leaf(., A=a[i] ) ###
#    }
#    print(c(''))
#    fe <- .$rs_fe()
#    print(c('f(e) <= 1.6','g0 only calulation:', afinal, ', f(e)', fe ))
#    print(cbind(a,resid))
#    write.csv(data.frame(afinal, fe, a, resid), 'fe_break_resid.csv', append=T, row.names=F, quote=F )
#    # if f(e) < 1.6 then -999 is returned, triggering a calcultaion where gs = g0, 
#    # but this doesn't always give the correct solution, A is often > 0 so gs > g0
#    # could use full soln here, or alternative method to get first guess, quadratic solve?
#    return(-999) 
#  } else if(.super$state$aguess[1] < 0) return(-99)
#  else if(abs(.super$state$faguess[1]) < 1e-6) return(.super$state$aguess[1])
#  else {
# 
#    # second guess, also analytical
#    #.super$state$cc[] <- .$gas_diff( .super$state$aguess[1] , 
#    #  r=( 1.4*.super$state_pars$rb + 1.6*.$rs(A=.super$state$aguess[1],c=.$gas_diff(.super$state$aguess[1])) + .super$state_pars$ri ) )
#    #.super$state$aguess[2]  <- .$assimilation() 
#    .super$state$aguess[2]  <- .super$state$faguess[1] + .super$state$aguess[1]  
#    if(.super$state$cc<0&F) {
#      print('')
#      print(c(.super$state$cc, .super$state$aguess[1], .super$state$aguess[2] ))
#    }
#    #.super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 
#  
#
#  # if both a1 and a2 are below zero then return -1 and leafsys routine will calculate A assuming gs = g0
#  #if(.super$state$cc < 0) return(-1.01)
#  if(.super$state$aguess[1]<0 & (.super$state$aguess[2]<0|.super$state$cc<0) ) return(-1.01)
#  else {
#
#    # calculate residual of second guess
#    .super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 
#
#    # third guess
#    #.super$state$aguess[3]  <- mean(.super$state$aguess[1:2])  
#    .super$state$aguess[3]  <- .super$state$faguess[2] + .super$state$aguess[2]  
#    .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
#
#    # check that fa2 and fa1 span 0, if not reguess for fa3 until it and fa1 span zero and then take the mean a1 and a3 to find a2
#    # under conditions tested so far this is only necessary for Ball-Berry gs
#    if(prod(.super$state$faguess[1:2]) > 0) {
#      print('')
#      print('fa1 and fa2 do not span zero quad')
#      print(.super$state_pars); print(.super$env) 
#      print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
#      print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))
#     
#
#      if(.super$state$faguess[1] > 0 ) { 
#        while(.super$state$faguess[3]>0) {
#          .super$state$aguess[3]  <- .super$state$aguess[3] + 0.01  
#          .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
#        }
#      } else {
#        while(.super$state$faguess[3] < 0) {
#          .super$state$aguess[3]  <- .super$state$aguess[3] - 0.01  
#          .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
#        }
#      }
#
#      .super$state$aguess[2]  <- mean(.super$state$aguess[c(1,3)])
#      .super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 
#
#      # reporting for development 
#      ### can call residual function in here   
#      a  <- seq(.super$state$aguess[1]-2,.super$state$aguess[1]+2,0.1 )
#      resid <- numeric(length(a))
#      for( i in 1:length(a) ) {
#        resid[i] <- f_A_r_leaf(., A=a[i] )
#      }
#      print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
#      print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))
#      print(cbind(a,resid))
#
#      # over-cautious safety - with while condition above should be no way that this condition is true
#      if(prod(.super$state$faguess[c(1,3)]) > 0) stop('failed')
#    }
#
#    # check f(guesses) span the root, if not print diagnostics and stop
#    if( min(.super$state$faguess[1:3]) * max(.super$state$faguess[1:3]) >= 0 ) {
#    print(unlist(.super$fnames)); print(unlist(.super$pars)); print(unlist(.super$state_pars)); print(unlist(.super$env)) 
#    print(c('faguess 1-3,',round(.super$state$faguess[1:3],4),'aguess 1-3,',round(.super$state$aguess[1:3],4)))
#    ### can call residual function in here   
#    sinput <- seq(.super$state$aguess[1]-2,.super$state$aguess[1]+2,0.1 )
#    out    <- numeric(length(sinput))
#    for( i in 1:length(sinput) ) {
#      out[i] <- f_A_r_leaf(., A=sinput[i] )
#    }
#      print(cbind(sinput,out))
#      stop(paste('Solver error'))
#    } 
#   
#    # fit a quadratic through the three sets of co-ordinates a la Lomas (one step of Muller's method) 
#    # Muller, David E., "A Method for Solving Algebraic Equations Using an Automatic Computer," Mathematical Tables and Other Aids to Computation, 10 (1956)    
#    #bx <- ((a1+a0)*(fa2-fa1)/(a2-a1)-(a2+a1)*(fa1-fa0)/(a1-a0))/(a0-a2)
#    #ax <- (fa2-fa1-bx*(a2-a1))/(a2**2-a1**2)
#    #cx <- fa1-bx*a1-ax*a1**2
#
#    bx <- ((.super$state$aguess[2]+.super$state$aguess[1])*(.super$state$faguess[3]-.super$state$faguess[2]) / 
#      (.super$state$aguess[3]-.super$state$aguess[2])-(.super$state$aguess[3]+.super$state$aguess[2])*(.super$state$faguess[2]-.super$state$faguess[1])/(.super$state$aguess[2]-.super$state$aguess[1])) / 
#      (.super$state$aguess[1]-.super$state$aguess[3])
#    ax <- (.super$state$faguess[3]-.super$state$faguess[2]-bx*(.super$state$aguess[3]-.super$state$aguess[2]))/(.super$state$aguess[3]**2-.super$state$aguess[2]**2)
#    cx <- .super$state$faguess[2]-bx*.super$state$aguess[2]-ax*.super$state$aguess[2]**2
# 
#    # find the root of the quadratic that lies between guesses 
#    .super$state$assim[] <- quad_sol(ax,bx,cx,'both')
#    ss <- which(.super$state$assim>min(.super$state$aguess[1:3])&.super$state$assim<max(.super$state$aguess[1:3]))
# 
#    # reporting for development 
#    .super$state$fA_ana_final[1] <- .$solver_func( .super$state$assim[1] ) 
#    .super$state$fA_ana_final[2] <- .$solver_func( .super$state$assim[2] ) 
#    .super$state$aguess[4]  <- .super$state$assim[ss]  
#    .super$state$faguess[4] <- .$solver_func( .super$state$assim[ss]) 
#    
#    # catch potential errors
#    if(length(ss)==0)      stop('no solution,',.super$state$assim,'; within initial 3 guesses,',.super$state$aguess[1:3],'fguesses,',.super$state$faguess[1:3] ) 
#    else if(length(ss)==2) stop('both solutions,',.super$state$assim,'; within initial 3 guesses,',.super$state$aguess[1:3],'fguesses,',.super$state$faguess[1:3] ) 
#    else return(.super$state$assim[ss])
#  }}
#}
#
#
## Clean semi-analytical solution 
#f_A_r_leaf_semiana_brent <- function(.) {
#  # - finds the analytical solution assuming rb and ri are zero to use as first guess (a1)
#  # - make a second guess (a2) by calculating Cc using a1 and actual values of rb and ri   
#  # - calculate residual function value for these two guesses (fa1, fa2) 
#  # - check these span the root
#  # - if not, make a third guess (a3) by taking the mean of a1 and a2 and calculate residual (fa3) 
#  # - check fa1 and fa3 span the root
#  # - if not, sequentially add 0.01 to the third guess (a3) until fa1 and fa3 do span the root
#  # - pass the two guesses (either a1 and a2 OR a1 and a3) as brackets to the Brent solver 
# 
#  # find the analytical solution assuming rb and ri are zero to use as first guess (a1)
#  #.super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- f_A_r_leaf_analytical_quad(.)
#  ###.super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- f_A_r_leaf_analytical(.)
#  .super$state$aguess[1]  <- .super$state$A_ana_rbzero[] <- .$analytical()
#  ###.super$state$faguess[1] <- f_A_r_leaf(., .super$state$aguess[1] ) 
#  .super$state$faguess[1] <- .$full_soln(.super$state$aguess[1]) 
#  #.super$state$faguess[1] <- .$solver_func( .super$state$aguess[1] ) 
# 
#  if(.super$state$aguess[1] == -999 ) return(-1.01) 
#  else if(.super$state$aguess[1] < 0) return(-1.01)
#  else if(abs(.super$state$faguess[1]) < 1e-6) return(.super$state$aguess[1])
#  else {
# 
#    # second guess, also analytical
#    #.super$state$cc[] <- .$gas_diff( .super$state$aguess[1] , 
#    #  r=( 1.4*.super$state_pars$rb + 1.6*.$rs(A=.super$state$aguess[1],c=.$gas_diff(.super$state$aguess[1])) + .super$state_pars$ri ) )
#    #.super$state$aguess[2]  <- .$assimilation() 
#    #.super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 
#    .super$state$aguess[2]  <- .super$state$faguess[1] + .super$state$aguess[1]  
#
#  # if both a1 and a2 are below zero then return -1 and leafsys routine will calculate A assuming gs = g0
#  #if(.super$state$aguess[1]<0 & .super$state$aguess[2]<0) return(-1.01)
#  #else {
#  if(.super$state$aguess[1]<0 & (.super$state$aguess[2]<0|.super$state$cc<0) ) return(-1.01)
#  #if(.super$state$cc < 0) return(-1.01)
#  else {
#
#    # calculate residual of second guess
#    .super$state$faguess[2] <- .$solver_func( .super$state$aguess[2] ) 
# 
#    # check that fa2 and fa1 span 0, if not reguess for fa3 until it and fa1 span zero and then take the mean a1 and a3 to find a2
#    if(prod(.super$state$faguess[1:2]) > 0) {
#      print('fa1 and fa2 do not span zero brent')
#      
#      # third guess
#      .super$state$aguess[3]  <- mean(.super$state$aguess[1:2])  
#      #.super$state$aguess[3]  <- .super$state$faguess[2] + .super$state$aguess[2]  
#      .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
#      print(.super$state_pars); print(.super$env) 
#      print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
#      print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))
#      
#      while(.super$state$faguess[3]>0) {
#        .super$state$aguess[3]  <- .super$state$aguess[3] + 0.01  
#        .super$state$faguess[3] <- .$solver_func( .super$state$aguess[3] ) 
#      }
#
#      .super$state$aguess[2]  <- mean(.super$state$aguess[c(1,3)])
#      .super$state$faguess[2] <- .$solver_func(.super$state$aguess[2]) 
#      
#      # reporting for development
#      ### can use residual function in here 
#      a  <- seq(.super$state$aguess[1]-2,.super$state$aguess[1]+2,0.1 )
#      resid <- numeric(length(a))
#      for( i in 1:length(a) ) {
#        resid[i] <- f_A_r_leaf(., A=a[i] )
#      }
#      print(c('faguess 1-3:',round(.super$state$faguess[1:3],4)))
#      print(c( 'aguess 1-3:',round(.super$state$aguess[1:3],4)))
#      print(cbind(a,resid))
#  
#      # over-cautious safety - with while condition above should be no way that this condition is true
#      if(prod(.super$state$faguess[c(1,3)]) > 0) stop('failed')
#     
#      # return root brackets  
#      brackets <- .super$state$aguess[c(1,3)]
#    } else {
#      brackets <- .super$state$aguess[c(1,2)]
#    }
#  
#    # Now brackets are available that span the root of the residual function, find the root using the Brent solver
#    ###f_R_Brent_solver_diag_brackets(., brackets )
#    .$solver_brackets(brackets)
#  }}
#}



### ANALYTICAL SOLUTIONS
################################

# solves A analytically by assuming g0 = 0 in the stomatal resistance function, and rb and ri also = zero 
#f_solver_analytical_leaf_simple <- function(.) {
f_A_r_leaf_analytical <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes
  
  # calculate cb, ci & cc
  .super$state$cb <- .super$state$ca
  fe              <- .$rs_fe() 
  if(fe<=1e-1) {
    print(c('Ana f(e)',fe)) 
    return(-999)
  } else {
    .super$state$ci      <- .super$state$ca * (1 - (1.6 / fe) )
    .super$state$cc      <- .super$state$ci 
 
    # calculate net A
    Anet <- .$assimilation()
  
    # calculate rs
    .super$state_pars$rs <- .super$state$ca / (fe * Anet * .super$env$atm_press*1e-6) 
  
    # return net A
    return(Anet)
  }
}


# solves A analytically by assuming rb and ri are zero 
#f_solver_analytical_leaf_quad <- function(.) {
f_A_r_leaf_analytical_quad <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes

  # set cb, rb & ri
  .super$state$cb[]      <- .super$state$ca
  #.super$state_pars$rb[] <- 0
  #.super$state_pars$ri[] <- 0

  # correctly assign g0 if rs functions assume g0 = 0
  #g0_hold <- .super$pars$g0
  g0  <- if(.super$fnames$rs=='f_rs_cox1998'|.super$fnames$rs=='f_rs_constantCiCa') 0 else .super$pars$g0 
  gsd <- .$rs_fe() / .super$state$ca
  p   <- .super$env$atm_press*1e-6

  # calculate coefficients of quadratic to solve A
  .$assim_quad_soln <- function(., V, K, g0, gsd, p ) {
    a   <- p*( 1.6 - gsd*(.super$state$ca + K) )
    b   <- p*gsd*( .super$state$ca*(V - .super$state$rd) - .super$state$rd*K - V*.super$state_pars$gstar ) - g0*(.super$state$ca + K) + 1.6*p*(.super$state$rd - V)
    c   <- g0*( V*(.super$state$ca - .super$state_pars$gstar) - .super$state$rd*(K + .super$state$ca) )
 
    # return A 
    quad_sol(a,b,c,'upper')
  }

  # calculate Anet for each limiting rate, use Agross variables as place holders
  .super$state$Acg[] <- .$assim_quad_soln(V=.super$state_pars$vcmaxlt, K=.super$state_pars$Km,                                   g0=g0, gsd=gsd, p=p )
  .super$state$Ajg[] <- .$assim_quad_soln(V=(.super$state$J/4),        K=(2*.super$state_pars$gstar),                            g0=g0, gsd=gsd, p=p )
  .super$state$Apg[] <- .$assim_quad_soln(V=(3*.super$state_pars$tpu), K=(-(1+3*.super$pars$Apg_alpha)*.super$state_pars$gstar), g0=g0, gsd=gsd, p=p )

  # determine rate limiting cycle 
  # - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  # - check which is correct here
  Amin        <- .$Alim() 
  #Amin        <- .$Alim() - .super$state$rd
  if(Amin==0) Amin <- 1e-6 
  
  # determine cc/ci based on Amin
  # calculate rs
  #.super$pars$g0[]     <- g0_hold 
  .super$state_pars$rs[] <- .$rs(A=Amin)
  .super$state$cc[] <-.super$state$ci[] <- .$gas_diff(A=Amin, r=1.6*.super$state_pars$rs )
    
  # recalculate Ag for each limiting process
  # necessary if Alim is Collatz smoothing as it reduces A, decoupling A from cc calculated in the quadratic solution 
  .super$state$Acg[] <- .$Acg() * .super$state$cc
  .super$state$Ajg[] <- .$Ajg() * .super$state$cc
  .super$state$Apg[] <- .$Apg() * .super$state$cc

  # return net A
  Amin
}


# solves A analytically by assuming rs is equal to 1/g0 
#f_solver_analytical_leaf_quad_r0 <- function(.) {
f_A_r0_leaf_analytical_quad <- function(.) {
  # combines A, rs, ci, cc eqs to a single f(), 
  # combines all rate limiting processes

  r0 <- .$rs_r0() 
  r  <- 1.4*.super$state_pars$rb + 1.6*r0 + .super$state_pars$ri
  p  <- .super$env$atm_press*1e-6

  # calculate coefficients of quadratic to solve A
  .$assim_quad_soln <- function(.,V,K,r,p) {
    a   <- -p*r
    b   <- .super$state$ca + K - .super$state$rd*p*r + V*p*r
    c   <- .super$state$ca*(.super$state$rd-V) + .super$state$rd*K + V*.super$state_pars$gstar 
    
    # return A 
    quad_sol(a,b,c,'lower')
  }

  # calculate Anet for each limiting rate, use Agross variables as place holders
  .super$state$Acg[] <- .$assim_quad_soln(V=.super$state_pars$vcmaxlt, K=.super$state_pars$Km,                                   r=r, p=p ) + .super$state$rd
  .super$state$Ajg[] <- .$assim_quad_soln(V=(.super$state$J/4),        K=(2*.super$state_pars$gstar),                            r=r, p=p ) + .super$state$rd
  .super$state$Apg[] <- .$assim_quad_soln(V=(3*.super$state_pars$tpu), K=(-(1+3*.super$pars$Apg_alpha)*.super$state_pars$gstar), r=r, p=p ) + .super$state$rd

  # determine rate limiting cycle - this is done based on carboxylation, not net assimilation (Gu etal 2010).
  Amin        <- .$Alim() - .super$state$rd 
  
  # determine cc/ci based on Amin
  .super$state_pars$rs[] <- r0 
  .super$state$cb[]      <- .$gas_diff(A=Amin)
  .super$state$cc[]      <-.super$state$ci <- .$gas_diff(A=Amin, r=1.6*.super$state_pars$rs, c=.super$state$cb )
    
  # recalculate Ag for each limiting process
  # necessary if Alim is Collatz smoothing as it reduces A, decoupling A from cc calculated in the quadratic solution 
  .super$state$Acg[] <- .$Acg() * .super$state$cc
  .super$state$Ajg[] <- .$Ajg() * .super$state$cc
  .super$state$Apg[] <- .$Apg() * .super$state$cc

  # return net A
  Amin
}


# Calculate assimilation assuming zero resistance to CO2 diffusion from the atmosphere to the site of carboxylation
#f_solver_analytical_leaf_no_r <- function(.,...) {
f_A_r_leaf_noR <- function(.,...) {
  # This function can be used to calculate the stomatal limitation to photosynthesis when rb and ri are assumed zero 
  # (when ri and rb are non-zero, use below function 'f_A_r_leaf_noRs' )
  
  # combines all rate limiting processes
  
  # assume cc = ca
  .super$state$cc <- .super$state$ca
  
  # calculate net A
  .$assimilation()
}



### END ###
