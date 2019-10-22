################################
#
# General functions for MAAT model
#
# AWalker August 2015
#
################################

library(deSolve)


#####################################################
# uniroot function editing to work as a proto method

puniroot <-
  function (., f, interval, ..., lower = min(interval), upper = max(interval),
            f.lower = f(lower, ...), f.upper = f(upper, ...),
            extendInt = c("no", "yes", "downX", "upX"), check.conv = FALSE, tol = .Machine$double.eps^0.25,
            maxiter = 1000, trace = 0, NAproceed = T ) {
    if (!missing(interval) && length(interval) != 2L)
      stop("'interval' must be a vector of length 2")
    if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
      stop("lower < upper  is not fulfilled")
    if (is.na(f.lower))
      #if(NAproceed) {print("f.lower = f(lower) is NA"); print(paste("f.upper = ",f.upper)); f.lower <- f.upper } else
      if(NAproceed) {
        print("f.lower = f(lower) is NA")
        return(list(root = NA, f.root = NA, iter = 0, init.it = 0, estim.prec = NA))
      } else {
        stop("f.lower = f(lower) is NA")
      }
    if (is.na(f.upper))
      #if(NAproceed) {print("f.upper = f(upper) is NA"); print(paste("f.lower = ",f.lower)); f.upper <- f(20.3453245) } else
      if(NAproceed) {
        print("f.upper = f(upper) is NA")
        return(list(root = NA, f.root = NA, iter = 0, init.it = 0, estim.prec = NA))
      } else {
        stop("f.upper = f(upper) is NA")
      }
    Sig <- switch(match.arg(extendInt), yes = NULL, downX = -1, no = 0, upX = 1, stop("invalid 'extendInt'; please report"))
    truncate <- function(x) pmax.int(pmin(x, .Machine$double.xmax), -.Machine$double.xmax)
    f.low. <- truncate(f.lower)
    f.upp. <- truncate(f.upper)
    doX <- (is.null(Sig) && f.low. * f.upp. > 0 || is.numeric(Sig) && (Sig * f.low. > 0 || Sig * f.upp. < 0))
    if (doX) {
      if (trace)
        cat(sprintf("search in [%g,%g]%s", lower, upper,
                    if (trace >= 2)
                      "\\n"
                    else " ... "))
      Delta <- function(u) 0.01 * pmax(1e-04, abs(u))
      it <- 0L
      if (is.null(Sig)) {
        delta <- Delta(c(lower, upper))
        while (isTRUE(f.lower * f.upper > 0) && any(iF <- is.finite(c(lower,upper)))) {
          if ((it <- it + 1L) > maxiter) stop(gettextf("no sign change found in %d iterations", it - 1), domain = NA)
          if (iF[1]) {
            ol <- lower
            of <- f.lower
            if (is.na(f.lower <- f(lower <- lower - delta[1], ...))) {
              lower <- ol
              f.lower <- of
              delta[1] <- delta[1]/4
            }
          }
          if (iF[2]) {
            ol <- upper
            of <- f.upper
            if (is.na(f.upper <- f(upper <- upper + delta[2], ...))) {
              upper <- ol
              f.upper <- of
              delta[2] <- delta[2]/4
            }
          }
          if (trace >= 2)
            cat(sprintf(" .. modified lower,upper: (%15g,%15g)\\n", lower, upper))
          delta <- 2 * delta
        }
      }
      else {
        delta <- Delta(lower)
        while (isTRUE(Sig * f.lower > 0)) {
          if ((it <- it + 1L) > maxiter) stop(gettextf("no sign change found in %d iterations", it - 1), domain = NA)
          f.lower <- f(lower <- lower - delta, ...)
          if (trace >= 2) cat(sprintf(" .. modified lower: %g\\n", lower))
          delta <- 2 * delta
        }
        delta <- Delta(upper)
        while (isTRUE(Sig * f.upper < 0)) {
          if ((it <- it + 1L) > maxiter) stop(gettextf("no sign change found in %d iterations", it - 1), domain = NA)
          f.upper <- f(upper <- upper + delta, ...)
          if (trace >= 2) cat(sprintf(" .. modified upper: %g\\n", upper))
          delta <- 2 * delta
        }
      }
      if (trace && trace < 2)
        cat(sprintf("extended to [%g, %g] in %d steps\\n", lower, upper, it))
    }
    if (!isTRUE(as.vector(sign(f.lower) * sign(f.upper) <= 0)))
      stop(if (doX)
        "did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0"
        else "f() values at end points not of opposite sign")
    if (check.conv) {
      val <- tryCatch(.External2(stats:::C_zeroin2,
                                 function(arg) f(arg, ...), lower, upper, f.lower, f.upper, tol, as.integer(maxiter)),
                      warning = function(w) w)
      if (inherits(val, "warning"))
        stop("convergence problem in zero finding: ", conditionMessage(val))
    }
    else {
      val <- .External2(stats:::C_zeroin2, function(arg) f(arg, ...), lower, upper, f.lower, f.upper, tol, as.integer(maxiter))
    }
    iter <- as.integer(val[2L])
    if (iter < 0) {
      (if (check.conv)
        stop
       else warning)(sprintf(ngettext(maxiter, "_NOT_ converged in %d iteration",
                                      "_NOT_ converged in %d iterations"), maxiter), domain = NA)
      iter <- maxiter
    }
    if (doX)
      iter <- iter + it
    else it <- NA_integer_
    list(root = val[1L], f.root = f(val[1L], ...), iter = iter,
         init.it = it, estim.prec = val[3L])
}



##############################
# deSolve lsoda function editing to work as a proto method

plsoda <- 
  # function (., y, times, func, parms, rtol = 1e-06, atol = 1e-06, 
  function (., y, times, func, parms=NULL, rtol = 1e-06, atol = 1e-06, 
            jacfunc = NULL, jactype = "fullint", rootfunc = NULL, verbose = FALSE, 
            nroot = 0, tcrit = NULL, hmin = 0, hmax = NULL, hini = 0, 
            # ynames = TRUE, maxordn = 12, maxords = 5, bandup = NULL, 
            ynames = FALSE, maxordn = 12, maxords = 5, bandup = NULL, 
            banddown = NULL, maxsteps = 5000, dllname = NULL, initfunc = dllname, 
            initpar = parms, rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, 
            forcings = NULL, initforc = NULL, fcontrol = NULL, events = NULL, 
            lags = NULL, ...) 
  {
    # if (!is.null(rootfunc)) 
    #   return(lsodar(y, times, func, parms, rtol, atol, jacfunc, 
    #                 jactype, rootfunc, verbose, nroot, tcrit, hmin, hmax, 
    #                 hini, ynames, maxordn, maxords, bandup, banddown, 
    #                 maxsteps, dllname, initfunc, initpar, rpar, ipar, 
    #                 nout, outnames, forcings, initforc, fcontrol, events, 
    #                 lags, ...))
    # if (is.list(func)) {
    #   if (!is.null(jacfunc) & "jacfunc" %in% names(func)) 
    #     stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
    #   if (!is.null(initfunc) & "initfunc" %in% names(func)) 
    #     stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
    #   if (!is.null(dllname) & "dllname" %in% names(func)) 
    #     stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
    #   if (!is.null(initforc) & "initforc" %in% names(func)) 
    #     stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
    #   if (!is.null(events$func) & "eventfunc" %in% names(func)) 
    #     stop("If 'func' is a list that contains eventfunc, argument 'events$func' should be NULL")
    #   if ("eventfunc" %in% names(func)) {
    #     if (!is.null(events)) 
    #       events$func <- func$eventfunc
    #     else events <- list(func = func$eventfunc)
    #   }
    #   if (!is.null(func$jacfunc)) 
    #     jacfunc <- func$jacfunc
    #   if (!is.null(func$initfunc)) 
    #     initfunc <- func$initfunc
    #   if (!is.null(func$dllname)) 
    #     dllname <- func$dllname
    #   if (!is.null(func$initforc)) 
    #     initforc <- func$initforc
    #   func <- func$func
    # }
    print(hmax)
    hmax <- deSolve:::checkInput(y, times, func, rtol, atol, jacfunc, tcrit, 
                       hmin, hmax, hini, dllname)
    n <- length(y)
    # if (!is.numeric(maxordn)) 
    #   stop("`maxordn' must be numeric")
    # if (maxordn < 1 || maxordn > 12) 
    #   stop("`maxord' must be >1 and <=12")
    # if (!is.numeric(maxords)) 
    #   stop("`maxords' must be numeric")
    # if (maxords < 1 || maxords > 5) 
    #   stop("`maxords' must be >1 and <=5")
    if (jactype == "fullint") 
      jt <- 2
    # else if (jactype == "fullusr") 
    #   jt <- 1
    # else if (jactype == "bandusr") 
    #   jt <- 4
    # else if (jactype == "bandint") 
    #   jt <- 5
    # else stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr' or 'bandint'")
    # if (jt %in% c(4, 5) && is.null(bandup)) 
    #   stop("'bandup' must be specified if banded Jacobian")
    # if (jt %in% c(4, 5) && is.null(banddown)) 
    #   stop("'banddown' must be specified if banded Jacobian")
    if (is.null(banddown)) 
      banddown <- 1
    if (is.null(bandup)) 
      bandup <- 1
    # if (jt %in% c(1, 4) && is.null(jacfunc)) 
    #   stop("'jacfunc' NOT specified; either specify 'jacfunc' or change 'jactype'")
    Ynames    <- attr(y, "names")
    JacFunc   <- NULL
    flist     <- list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    ModelInit <- NULL
    Eventfunc <- NULL
    events    <- deSolve:::checkevents(events, times, Ynames, dllname)
    #print(events)
    if (!is.null(events$newTimes)) 
      times <- events$newTimes
    if (jt == 4 && banddown > 0) 
      erow <- matrix(data = 0, ncol = n, nrow = banddown)
    else erow <- NULL
    # if (is.character(func) | class(func) == "CFunc") {
    #   DLL <- checkDLL(func, jacfunc, dllname, initfunc, verbose, 
    #                   nout, outnames)
    #   ModelInit <- DLL$ModelInit
    #   Func <- DLL$Func
    #   JacFunc <- DLL$JacFunc
    #   Nglobal <- DLL$Nglobal
    #   Nmtot <- DLL$Nmtot
    #   if (!is.null(forcings)) 
    #     flist <- checkforcings(forcings, times, dllname, 
    #                            initforc, verbose, fcontrol)
    #   if (is.null(ipar)) 
    #     ipar <- 0
    #   if (is.null(rpar)) 
    #     rpar <- 0
    #   Eventfunc <- events$func
    #   if (is.function(Eventfunc)) 
    #     rho <- environment(Eventfunc)
    #   else rho <- NULL
    # }
    # else {
      if (is.null(initfunc)) 
        initpar <- NULL
      rho <- environment(func)
      # if (ynames) {
      #   Func <- function(time, state) {
      #     attr(state, "names") <- Ynames
      #     unlist(func(time, state, parms, ...))
      #   }
      #   Func2 <- function(time, state) {
      #     attr(state, "names") <- Ynames
      #     func(time, state, parms, ...)
      #   }
      #   JacFunc <- function(time, state) {
      #     attr(state, "names") <- Ynames
      #     rbind(jacfunc(time, state, parms, ...), erow)
      #   }
      #   if (!is.null(events$Type)) 
      #     if (events$Type == 2) 
      #       Eventfunc <- function(time, state) {
      #         attr(state, "names") <- Ynames
      #         events$func(time, state, parms, ...)
      #       }
      # }
      # else {
        Func    <- function(time, state) unlist(func(time, state, parms, ...))
        Func2   <- function(time, state) func(time, state, parms, ...)
        # Func2   <- function(time, state) func(time, state, ...)
        # Func2   <- function(., time, state) .$func(time, state, parms, ...)
        JacFunc <- function(time, state) rbind(jacfunc(time, state, parms, ...), erow)
        
        if (!is.null(events$Type)) 
          if (events$Type == 2) 
            Eventfunc <- function(time, state) events$func(time, state, parms, ...)
      # }
      FF      <- deSolve:::checkFunc(Func2, times, y, rho) # this works
      # FF      <- deSolve:::checkFunc(.$Func2, times, y, rho)
      # FF      <- deSolve:::checkFunc(func, times, y, rho) # this works
      # FF      <- deSolve:::checkFunc(.$lsexamp, times, y, rho) # this works
      # FF      <- .$pcheckFunc(times, y, rho) # this works
      Nglobal <- FF$Nglobal
      Nmtot   <- FF$Nmtot
      #print(FF)
      # print(events)
      # if (!is.null(events$Type)) 
      #   if (events$Type == 2) 
      #     checkEventFunc(Eventfunc, times, y, rho)
      # if (jt %in% c(1, 4)) {
      #   tmp <- eval(JacFunc(times[1], y), rho)
      #   if (!is.matrix(tmp)) 
      #     stop("Jacobian function, 'jacfunc' must return a matrix\n")
      #   dd <- dim(tmp)
      #   if ((jt == 4 && any(dd != c(bandup + banddown + banddown + 
      #                               1, n))) || (jt == 1 && any(dd != c(n, n)))) 
      #     stop("Jacobian dimension not ok")
      # }
    # }
    if (jt %in% c(1, 2)) 
      lmat <- n^2 + 2
    # else if (jt %in% c(4, 5)) 
    #   lmat <- (2 * banddown + bandup + 1) * n + 2
    lrn = 20 + n * (maxordn + 1) + 3 * n
    lrs = 20 + n * (maxords + 1) + 3 * n + lmat
    lrw = max(lrn, lrs)
    liw = 20 + n
    iwork <- vector("integer", 20)
    rwork <- vector("double", 20)
    rwork[] <- 0
    iwork[] <- 0
    iwork[1] <- banddown
    iwork[2] <- bandup
    iwork[6] <- maxsteps
    if (maxordn != 12) 
      iwork[8] <- maxordn
    if (maxords != 5) 
      iwork[9] <- maxords
    if (verbose) 
      iwork[5] = 1
    if (!is.null(tcrit)) 
      rwork[1] <- tcrit
    rwork[5] <- hini
    rwork[6] <- hmax
    rwork[7] <- hmin
    if (!is.null(times)) 
      itask <- ifelse(is.null(tcrit), 1, 4)
    else itask <- ifelse(is.null(tcrit), 2, 5)
    if (is.null(times)) 
      times <- c(0, 1e+08)
    if (verbose) 
      printtask(itask, func, jacfunc)
    storage.mode(y) <- storage.mode(times) <- "double"
    IN <- 1
    lags <- deSolve:::checklags(lags, dllname)
    on.exit(.C("unlock_solver"))
    out <- .Call("call_lsoda", y, times, Func, initpar, rtol, 
                 atol, rho, tcrit, JacFunc, ModelInit, Eventfunc, as.integer(verbose), 
                 as.integer(itask), as.double(rwork), as.integer(iwork), 
                 as.integer(jt), as.integer(Nglobal), as.integer(lrw), 
                 as.integer(liw), as.integer(IN), NULL, 0L, as.double(rpar), 
                 as.integer(ipar), 0L, flist, events, lags, PACKAGE = "deSolve")
    out <- deSolve:::saveOut(out, y, n, Nglobal, Nmtot, func, Func2, 
                             iin = c(1,12:21), iout = c(1:3, 14, 5:9, 15:16), nr = 5)
    attr(out, "type") <- "lsoda"
    if (verbose) diagnostics(out)
    out
}



### END ###
