################################
#
# MAAT soil_decomp system representation functions (SRFs) 
# 
# AWalker, Matt Craig, October 2019 
#
################################



################################
# soil_decomp system function 

# a system of n pools, composed of a single dynamic solver function
f_sys_npools <- function(.) {
  
  # ------------------------------------------------------------
  # Save current state for mass-balance check
  # ------------------------------------------------------------
  
  .super$state_pars$previous_state <- .super$state$cpools
  
  
  # ------------------------------------------------------------
  # Calculate state-dependent parameters
  # ------------------------------------------------------------
  
  .$calc_state_pars()
  
  
  # ------------------------------------------------------------
  # Run solver
  #
  # Request output every 0.001 time units so that the trajectory
  # can be used to integrate the time-varying loss fluxes.
  # ------------------------------------------------------------
  
  solver_times <- seq(0, 1, by = 0.001)
  
  .super$state_pars$solver_out <-
    tryCatch(
      .$solver(
        .super$state$cpools,
        solver_times,
        .$solver_func
      ),
      error = function(c) {
        print(c)
        
        matrix(
          ncol = .super$pars$n_pools + 1
        )
      }
    )
  
  
  # ------------------------------------------------------------
  # Store solver output locally
  # ------------------------------------------------------------
  
  sol <- .super$state_pars$solver_out
  
  
  # ------------------------------------------------------------
  # Assign final solved values to state
  # ------------------------------------------------------------
  
  .super$state$cpools[,1] <-
    sol[
      nrow(sol),
      2:(.super$pars$n_pools + 1)
    ]
  
  
  # ------------------------------------------------------------
  # Calculate instantaneous external losses along the
  # entire solver trajectory
  # ------------------------------------------------------------
  
  loss_total <- numeric(nrow(sol))
  loss_leaching <- numeric(nrow(sol))
  
  
  for(j in seq_len(nrow(sol))) {
    
    # Current solver time
    t <- as.numeric(sol[j, 1])
    
    # Current pool state
    C <- as.numeric(sol[j, -1])
    
    
    # Calculate outfluxes for this state and time
    .$outfluxes(C, t)
    
    
    # Calculate transfer matrix
    T <- .$transfermatrix(C, t)
    
    
    # Calculate summed decomposition/outflux from each pool
    O <- as.numeric(.$DotO(C, t))
    
    
    # Fraction of each pool's outflux that leaves the system
    pool_loss_frac <- abs(colSums(T))
    
    pool_loss_frac[pool_loss_frac < 1e-15] <- 0
    
    
    # Total instantaneous external loss
    loss_total[j] <-
      sum(pool_loss_frac * O)
    
    
    # ----------------------------------------------------------
    # Leaching is outflux 8 in this model
    # ----------------------------------------------------------
    
    if(!is.na(.super$pars$leaching_outflux)) {
      
      leach_name <-
        paste0(
          "of",
          .super$pars$leaching_outflux
        )
      
      loss_leaching[j] <-
        .super$state$outflux[[leach_name]]
      
    } else {
      
      loss_leaching[j] <- 0
      
    }
    
  }
  
  
  # ------------------------------------------------------------
  # Integrate total external loss over the solver interval
  #
  # Trapezoidal integration:
  #
  # integral = sum(
  #   dt * (loss[j] + loss[j+1]) / 2
  # )
  # ------------------------------------------------------------
  
  integrated_loss <-
    sum(
      diff(sol[,1]) *
        (
          head(loss_total, -1) +
            tail(loss_total, -1)
        ) / 2
    )
  
  
  # ------------------------------------------------------------
  # Integrate leaching separately
  # ------------------------------------------------------------
  
  integrated_leaching <-
    sum(
      diff(sol[,1]) *
        (
          head(loss_leaching, -1) +
            tail(loss_leaching, -1)
        ) / 2
    )
  
  
  # ------------------------------------------------------------
  # Respiration = total external loss - leaching
  # ------------------------------------------------------------
  
  integrated_respiration <-
    integrated_loss - integrated_leaching
  
  
  # ------------------------------------------------------------
  # Save integrated quantities for mass balance
  # ------------------------------------------------------------
  
  .super$state_pars$integrated_loss <-
    integrated_loss
  
  .super$state_pars$integrated_leaching <-
    integrated_leaching
  
  .super$state_pars$integrated_respiration <-
    integrated_respiration
  
  
  # ------------------------------------------------------------
  # Calculate final-state instantaneous loss fluxes
  #
  # This uses the original f_calc_loss_fluxes().
  # These are useful as the model's final flux outputs, but
  # are NOT used for the annual mass-balance calculation.
  # ------------------------------------------------------------
  
  .$calc_loss_fluxes()
  
  
  # ------------------------------------------------------------
  # Mass balance check
  # ------------------------------------------------------------
  
  .$mass_balance()
  
}


# steady-state solver for a system of n pools
f_steadystate_npools <- function(.) {

  # calculate state parameters
  .$calc_state_pars()

  # run solver
  .super$state_pars$solver_steadystate_out <- 
    tryCatch(
      .$solver_steadystate(.super$state$cpools, 0, func=.$solver_func ),
      error = function(c) {
        print(c)
        list(y=rep(NA,.super$pars$n_pools))
      })

  # assign solved values to state
  .super$state$cpools[,1] <- .super$state_pars$solver_steadystate_out$y
 
  # calculate respiration and other loss fluxes (e.g. leaching)
  .$calc_loss_fluxes()

  # mass balance check 
  # - at steady state inputs = outputs
  .$mass_balance(steadystate=T)
}



# NULL solver to prevent steady-state initialisation, i.e. to initialise with values specified in input
f_steadystate_null   <- function(.) NULL 



### END ###
