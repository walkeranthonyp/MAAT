################################
#
# General functions for MAAT model
#
# AWalker August 2015
#
################################



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



### END ###
