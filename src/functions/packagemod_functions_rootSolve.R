################################
#
# General functions for MAAT model
#
# AWalker July 2021 
#
################################

library(rootSolve)


##############################
# rootSolve stode function editing to work as a proto method

pstode <- 
  function(., y, time=0, func, parms=NULL, rtol=1e-06, atol=1e-08, 
           ctol=1e-08, jacfunc=NULL, jactype="fullint", verbose=FALSE, 
           #bandup=1, banddown=1, positive=FALSE, maxiter=100, 
           bandup=1, banddown=1, positive=TRUE, maxiter=100, 
           ynames=TRUE, dllname=NULL, initfunc=dllname, initpar=parms, 
           rpar=NULL, ipar=NULL, nout=0, outnames=NULL, forcings=NULL, 
           initforc=NULL, fcontrol=NULL, times=time, ... ) {
  
    if (is.list(func)) {
        if (!is.null(jacfunc) & "jacfunc" %in% names(func)) 
            stop("If 'func' is a list that contains jacfunc, argument 'jacfunc' should be NULL")
        if (!is.null(initfunc) & "initfunc" %in% names(func)) 
            stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be NULL")
        if (!is.null(dllname) & "dllname" %in% names(func)) 
            stop("If 'func' is a list that contains dllname, argument 'dllname' should be NULL")
        if (!is.null(initforc) & "initforc" %in% names(func)) 
            stop("If 'func' is a list that contains initforc, argument 'initforc' should be NULL")
        if (!is.null(func$jacfunc)) 
            jacfunc <- func$jacfunc
        if (!is.null(func$initfunc)) 
            initfunc <- func$initfunc
        if (!is.null(func$dllname)) 
            dllname <- func$dllname
        if (!is.null(func$initforc)) 
            initforc <- func$initforc
        func <- func$func
    }
    if (!is.numeric(y)) 
        stop("`y' must be numeric")
    n <- length(y)
    if (!is.null(times) && !is.numeric(times)) 
        stop("`times' must be NULL or numeric")
    if (length(times) > 1) 
        warning("`times' should be one number - taking first value")
    time <- times[1]
    if (!rootSolve:::CheckFunc(func)) 
        stop("`func' must be a function or character vector or a compiled function")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname))) 
        stop("You need to specify the name of the dll or shared library where 'func' can be found (without extension)")
    if (!is.numeric(maxiter)) 
        stop("`maxiter' must be numeric")
    if (as.integer(maxiter) < 1) 
        stop("'maxiter' must be >=1")
    if (!is.numeric(rtol)) 
        stop("`rtol' must be numeric")
    if (!is.numeric(atol)) 
        stop("`atol' must be numeric")
    if (!is.numeric(ctol)) 
        stop("`ctol' must be numeric")
    if (!is.null(jacfunc) & !rootSolve:::CheckFunc(jacfunc)) 
        stop("`jacfunc' must be a function or character vector or a compiled function")
    if (length(atol) > 1 && length(atol) != n) 
        stop("`atol' must either be a scalar, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n) 
        stop("`rtol' must either be a scalar, or as long as `y'")
    if (length(ctol) > 1) 
        stop("`ctol' must be a scalar")
    itol <- 1
    if (length(atol) == n && length(rtol) == n) 
        itol <- 4
    else if (length(atol) == n && length(rtol) != n) 
        itol <- 2
    else if (length(atol) != n && length(rtol) == n) 
        itol <- 3
    if (jactype == "fullint") 
        imp <- 22
    else if (jactype == "fullusr") 
        imp <- 21
    else if (jactype == "bandusr") 
        imp <- 24
    else if (jactype == "bandint") 
        imp <- 25
    else if (jactype %in% c("1D", "1Dint")) 
        imp <- 0
    else stop("'jactype' must be one of 'fullint', 'fullusr', 'bandusr', or 'bandint'")
    if (imp == 0) {
        nspec <- bandup
        ndim <- banddown
        banddown <- nspec
    }
    else {
        nspec <- 0
        ndim <- 0
    }
    if (imp %in% c(21, 24) && is.null(jacfunc)) 
        stop("stode: cannot estimate steady-state: *jacfunc* NOT specified; either specify *jacfunc* or change *jactype*")
    if (imp %in% c(24, 25)) 
        nabd <- 1 + 2 * bandup + banddown
    else nabd <- n
    if (verbose) {
        print("Steady-state settings")
        if (is.character(func)) 
            print(paste("model function a DLL: ", func))
        if (is.character(jacfunc)) 
            print(paste("jacobian specified as a DLL: ", jacfunc))
        print("jacobian method")
        df <- c("method flag, mf", "jsv", "meth", "miter", "itask")
        if (imp == 22) 
            txt <- "full jacobian, calculated internally"
        else if (imp == 21) 
            txt <- "full jacobian, specified by user function"
        else if (imp == 24) 
            txt <- "banded jacobian, specified by user function"
        else if (imp == 0) 
            txt <- "banded jacobian, 1-D, specified internally"
        else txt <- "banded jacobian, calculated internally"
        print(data.frame("implicit method", value = imp, message = txt))
    }
    Ynames    <- attr(y, "names")
    JacFunc   <- NULL
    ModelInit <- NULL
    ModelForc <- NULL
    Forc      <- NULL
    if(rootSolve:::is.compiled(func)) {
        if(!is.null(initfunc)) {
            if(inherits(initfunc, "CFunc")) 
                ModelInit <- body(initfunc)[[2]]
            else if(is.loaded(initfunc, PACKAGE = dllname, type = "") || 
                is.loaded(initfunc, PACKAGE = dllname, type = "Fortran")) 
                ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
        }
        if (!is.null(initforc)) {
            if (inherits(initforc, "CFunc")) 
                ModelForc <- body(initforc)[[2]]
            else if (is.loaded(initforc, PACKAGE = dllname, type = "") || 
                is.loaded(initforc, PACKAGE = dllname, type = "Fortran")) 
                ModelForc <- getNativeSymbolInfo(initforc, PACKAGE = dllname)$address
            if (is.list(forcings)) {
                Forc <- NULL
                for (i in 1:length(forcings)) if (!is.null(fcontrol)) 
                  Forc <- c(Forc, do.call(approx, list(x = forcings[[i]], 
                    y = NULL, xout = times, fcontrol))$y)
                else Forc <- c(Forc, do.call(approx, list(x = forcings[[i]], 
                  y = NULL, xout = times))$y)
            }
            else Forc <- forcings
        }
    }
    if(rootSolve:::is.compiled(func)) {
        funcname <- func
        if (inherits(func, "CFunc")) 
            Func <- body(func)[[2]]
        else if (is.loaded(funcname, PACKAGE = dllname)) {
            Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
        }
        else stop(paste("cannot calculate steady-state: dyn function not loaded: ", 
            funcname))
        if (!is.null(jacfunc)) {
            if (!rootSolve:::is.compiled(jacfunc)) 
                stop("If 'func' is dynloaded, so must 'jacfunc' be")
            jacfuncname <- jacfunc
            if (inherits(jacfunc, "CFunc")) 
                JacFunc <- body(jacfunc)[[2]]
            else if (is.loaded(jacfuncname, PACKAGE = dllname)) {
                JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
            }
            else stop(paste("cannot calculate steady-state: jacobian function not loaded ", 
                jacfunc))
        }
        Nglobal <- nout
        rho <- NULL
        if (is.null(outnames)) {
            Nmtot <- NULL
        }
        else if (length(outnames) == nout) {
            Nmtot <- outnames
        }
        else if (length(outnames) > nout) 
            Nmtot <- outnames[1:nout]
        else Nmtot <- c(outnames, (length(outnames) + 1):nout)
    }
    else {
        rho <- environment(func)
        if (ynames) {
            Func <- function(time, state) {
                attr(state, "names") <- Ynames
                func(time, state, parms, ...)[1]
            }
            Func2 <- function(time, state) {
                attr(state, "names") <- Ynames
                func(time, state, parms, ...)
            }
            JacFunc <- function(time, state) {
                attr(state, "names") <- Ynames
                jacfunc(time, state, parms, ...)
            }
        }
        else {
            Func <- function(time, state) func(time, state, parms, 
                ...)[1]
            Func2 <- function(time, state) func(time, state, 
                parms, ...)
            JacFunc <- function(time, state) jacfunc(time, state, 
                parms, ...)
        }
        tmp <- eval(Func2(time, y), rho)
        if (!is.list(tmp)) 
            stop("Model function must return a list\n")
        if (length(tmp[[1]]) != length(y)) 
            stop(paste("The number of derivatives returned by 'func() (", 
                length(tmp[[1]]), "must equal the length of the initial conditions vector (", 
                length(y), ")", sep = ""))
        if (any(is.na(tmp[[1]]))) 
            stop("Model function must return a list of values, of which first element has length =length of y\n ")
        Nglobal <- if (length(tmp) > 1) 
            length(unlist(tmp[-1]))
        else 0
        Nmtot <- attr(unlist(tmp[-1]), "names")
        if (imp %in% c(21, 24)) {
            tmp <- eval(JacFunc(time, y), rho)
            if (!is.matrix(tmp)) 
                stop("Jacobian function must return a matrix\n")
            dd <- dim(tmp)
            if ((imp == 24 && any(dd != c(bandup + banddown + 
                1, n))) || (imp == 21 && any(dd != c(n, n)))) 
                stop("Jacobian dimension not ok")
        }
    }
    storage.mode(y) <- "double"
    storage.mode(rtol) <- storage.mode(atol) <- storage.mode(ctol) <- "double"
    if (is.null(ipar)) 
        ipar <- 0
    if (is.null(rpar)) 
        rpar <- 0
    Pos <- FALSE
    if (is.logical(positive)) {
        Pos <- positive
    }
    else {
        if (!is.vector(positive)) 
            stop("'positive' should either be TRUE/FALSE or\n           a VECTOR with indices to the state variables that have to be positive")
        if (max(positive) > n) 
            stop("the elements of 'positive' should be < the number of state variables")
        if (min(positive) < 1) 
            stop("the elements of 'positive' should be >0")
    }
    if (is.null(initfunc)) 
        initpar <- NULL
    out <- .Call("call_dsteady", y, as.double(time), Func, as.double(initpar), 
        as.double(Forc), ctol, atol, rtol, as.integer(itol), 
        rho, JacFunc, ModelInit, ModelForc, as.integer(verbose), 
        as.integer(imp), as.integer(bandup), as.integer(banddown), 
        as.integer(maxiter), as.integer(Pos), as.integer(positive), 
        as.integer(Nglobal), as.integer(nabd), as.integer(nspec), 
        as.integer(ndim), as.double(rpar), as.integer(ipar), 
        PACKAGE = "rootSolve")
    precis <- attr(out, "precis")
    steady <- attr(out, "steady")
    attributes(out) <- NULL
    if (Nglobal > 0) {
        if (!is.character(func) & !inherits(func, "CFunc")) {
            y <- out
            if (ynames) 
                attr(y, "names") <- Ynames
            out2 <- Func2(time, y)[-1]
            out <- c(list(y = y), out2)
        }
        else {
            var = out[(n + 1):(n + Nglobal)]
            cnames <- Nmtot
            unames <- unique(Nmtot)
            var <- lapply(unames, FUN = function(x) var[which(cnames == 
                x)])
            names(var) <- unames
            y <- out[1:n]
            out <- c(y = 1, var)
            out[[1]] <- y
        }
    }
    else {
        if (ynames) 
            attr(out, "names") <- Ynames
        out <- list(y = out)
    }
    attr(out, "precis") <- precis
    attr(out, "steady") <- (steady == 1)
    if (!steady) 
        warning("steady-state not reached")
    if (verbose) {
        print("precision at each steady state step")
        print(precis)
    }
    return(out)
}



### END ###
