################################
#
# General functions for MAAT model
#
# AWalker August 2015
#
################################

library(xtable)
library(XML)


#####################################################
transpose_list <- function(l) {
  # takes a list of lists
  # each element of the outer list must be of equal length
  # each corresponding element of the inner lists should be of the same type
  # returns a list of vectors of equal type

  lt <- lapply(1:length(l[[1]]), function(inner) unlist(lapply(l, function(e) e[[inner]]))  )
  names(lt) <- names(l[[1]])
  lt
}


#####################################################
readXML <- function(input.file=NULL) {
  settings.xml <- NULL

  # Parse input xml
  if (!is.null(input.file) && file.exists(input.file)) {
    settings.xml <- xmlParse(input.file)
    # convert the xml to a list
    settings.list <- xmlToList(settings.xml)

  } else {
    print("***** WARNING: readXML called but with no input file defined *****")
  }

  # make sure something was loaded
  if (is.null(settings.xml)) {
    #log.error("Did not find any settings file to load.")
    stop(paste('readXML called, input file empty'))
  }

  # Remove comment or NULL fields
  settings.list <- settings.list[settings.list !="NULL" ]

  # Evaluate and return XML as a list
  return(evalXMLlist(settings.list))

}


#####################################################
# run over all elements of a list, returning a list of the same structure
# evaluate character strings to r objects

evalXMLlist <- function(sublist) {

  # subset sublist to elements that are non-null
  sub <- which(sapply(sublist,is.null))
  if(length(sub)>0) sublist <- sublist[-sub]

  # recursively evaluate all elements of sublist
  for(i in 1:length(sublist)) {
    sublist[[i]] <-
      if(typeof(sublist[[i]]) == "list")  evalXMLlist(sublist[[i]]) else eval(parse(text=sublist[[i]]))
  }
  sublist
}


#####################################################
# save a list to an xml file
listtoXML <- function(fname,name,...) {
  # expects the argument - sublist - in ...

  rec <- function(node, sublist) {
    for(i in 1:length(sublist)) {
      if(!is.null(unlist(sublist[[i]]))) {
        child <- newXMLNode(names(sublist)[i], parent=node)

        if (typeof(sublist[[i]]) == "list") {
          if(is.null(names(sublist[[i]]))) {
            #warning(paste(fname,names(sublist)[i],'list or sublist defined with no named elements, likely an empty list'))
            sublist[[i]] <- NA
            names(sublist[[i]]) <- 'nada'
          }
          rec(child, sublist[[i]])
        }
        else {
          xmlValue(child) <- if(class(sublist[[i]])=='character') paste0("'",sublist[[i]],"'") else sublist[[i]]
        }
      }
    }
  }

  root <- newXMLNode(name)
  rec(root,...)
  saveXML(root,fname)
}

# overwrite values in mainlist with values in sublist
# - currently expects list to be of identical structure
# fuselists <- function(mainlist,sublist) {
#   for(i in 1:length(sublist)) {
#     mlsub <- which(names(mainlist)==names(sublist)[i])
#     if(length(mlsub)!=1) {
#       stop(paste("names mismatch when fusing initialisation lists:",mainlist,sublist))
#     }
#     mainlist[[mlsub]] <-
#       if(typeof(sublist[[i]]) == "list")  fuselists(mainlist[[mlsub]],sublist[[i]]) else sublist[[i]]
#   }
#   mainlist
# }


#####################################################
# This function recursively fuses two lists
# by overwriting any named entries in mainlist with the same names as those in the sublist
# the overwriting occurs within the equivalent hierarchy of each list
# i.e. a named element in the mainlist will not be overwritten by the same named element in the sublist if the elements are in different levels of the lists
fuselists <- function(mainlist,sublist) {
  # length of sublist
  l <- length(sublist)
  n <- 0

  # loop over sublist
  for(i in 1:l) {
    mlsub <- which(names(mainlist)==names(sublist)[i])
    if(length(mlsub)==1) {
      n <- n + 1
      mainlist[[mlsub]] <-
        if(typeof(sublist[[i]]) == "list") fuselists(mainlist[[mlsub]],sublist[[i]])
          else if(is.null(sublist[[i]]))  mainlist[[mlsub]]
          else                            sublist[[i]]
    } else {
      error_i <- if(exists('error_i')) c(error_i,i) else i
    }
    # print(c(n,l))
    # print(names(sublist)[i])
  }

  # if mlsub != 1 then stop
  # i.e. if a name is found in sublist that does not occur in mainlist
  #      or if a names is found multiple times in the mainlist
  if(n!=l) {
    print(c(l,n))
    # stop(paste("\n names mismatch when fusing initialisation lists:",names(mainlist)[mlsub],names(sublist)[i]))
    stop(paste('\n names mismatch when fusing initialisation lists, sublist element:', names(sublist)[error_i],
               '; not found in mainlist \n This is likely a mis-spelling of a variable name in the input file'))
  }
  mainlist
}


#####################################################
# write Latex doc
write_Latex <- function(obj, fname, func=write_list_Latex, call=F, landscape=F, ... ) {
  fnamec <- paste0(fname,'.tex')

  # if(landscape) write("\\documentclass[10pt,landscape]{article}",fnamec)
  # else          write("\\documentclass[10pt]{article}",fnamec)
  write("\\documentclass[10pt]{article}",fnamec)
  if(landscape) write("\\usepackage[margin=0.75in,landscape]{geometry}",fnamec,append=T)
  write("\\usepackage{booktabs}",fnamec,append=T)
  write("\\usepackage{placeins}",fnamec,append=T)
  write("\\usepackage{longtable}",fnamec,append=T)
  #   write("\\usepackage[section]{placeins}",fname,append=T)
  write("\\begin{document}",fnamec,append=T)

  func(obj, fnamec, ... )

  write("\\end{document}",fnamec,append=T)
  if(call) {
    system(paste('pdflatex',fnamec))
    system(paste('rm', paste0(fname,'.log'), paste0(fname,'.aux') ))
  }
}


#####################################################
write_table_Latex <- function(sub, fname, cnames=NULL, rnames=NULL, capname=NULL, ... ) {
  if(!is.null(cnames)) colnames(combtable)  <- cnames
  print.xtable.apw(sub, cap=capname, file=fname, append=T, include.rownames=!is.null(rnames), ... )
}


#####################################################
write_tablelist_Latex <- function(slist, fname, cnames=NULL, rnames=NULL, capname=NULL, ... ) {
  for(i in 1:length(slist)) {
    write_table_Latex(slist[[i]], fname, capname=capname[i], ... )
  }
}


#####################################################
# write list recursively
write_list_Latex <- function(sublist, fname, ... ) {
  for(i in 1:length(sublist)) {
      if(typeof(sublist[[i]]) == "list") {
        write_list_Latex(sublist[[i]], fname, ... )
      } else {
        if(class(sublist[[i]]) == "numeric") print.xtable.apw(t(as.matrix(sublist[[i]])), cap=names(sublist)[1], file=fname, append=T, ... )
        else                                 print.xtable.apw(sublist[[i]], cap=names(sublist)[1], file=fname, append=T, ... )
      }

  }
  write("\\FloatBarrier",fname,append=T)
}


#####################################################
# write list recursively, combine lowest level list
write_list_Latex_comb <- function(sublist,fname,cnames=NULL,rnames=NULL,...) {
  if(typeof(sublist[[1]]) != "list") {
    write_list_Latex(sublist,fname)

  } else if(typeof(sublist[[1]][[1]]) == "list") {
    for(i in 1:length(sublist)) {
      write_list_Latex_comb(sublist[[i]],fname,cnames,rnames)
    }

  } else {
    for(li in 1:length(sublist[[1]])) {
      combtable <- do.call('rbind',lapply(sublist,function(l) l[[li]]))
      print(combtable)
      print(cnames)
      print(rnames)
      if(length(rnames)==dim(combtable)[1]) row.names(combtable) <- rnames
      if(length(cnames)==dim(combtable)[2]) colnames(combtable)  <- cnames
      print.xtable.apw(combtable,cap=names(sublist[[1]])[li],file=fname,append=T)
    }
  }
  write("\\FloatBarrier",fname,append=T)
}


#####################################################
print.xtable.apw <- function(x, cap=NULL, ... ){

  print(xtable(x,caption=cap,...),
        # floating=T,
        caption.placement='top',
        hline.after=NULL,
        add.to.row=list(pos=list(-1,0, nrow(x)),
                        command=c(
                          '\\toprule\n',
                          '\\midrule\n',
                          '\\bottomrule\n')),
        table.placement='h',
        ... )
}


#####################################################
# covert stacked data frame to table for latex table output
sens_table <- function(df1) {
  s1    <- unstack(df1, sensitivity1~variable)
  tm    <- apply(as.matrix(unstack(df1, mean~variable)),1,mean)
  tv    <- apply(as.matrix(unstack(df1, variance~variable)),1,mean)
  scen  <- unstack(df1, scenario~variable)[,1]
  df2   <- data.frame(scenario=scen, mean=tm, variance=tv, s1 )
  if(any(names(df1)=='model')) {
    mod  <- unstack(df1, model~variable)[,1]
    df2  <- data.frame(model=mod, df2 )
  }
  df2
}


# add environmetal variables to sensitivity output matrix
add_scenario_values <- function(df1) {
  dfe <- df1['scenario']
  if(!is.null(evar1)&!is.null(evar2)) {
    dfe <- cbind(dfe,dfe)
    names(dfe) <- c(evar1name, evar2name )
    dfe1 <- dfe
    for(r in 1:length(dfe[,1])) dfe1[r,] <- if(dfe[r,1]>0) cbind(evar1,evar2)[dfe[r,1],] else rep('int',2)
  } else if(!is.null(evar1)) {
    names(dfe) <- evar1name
    dfe1 <- dfe
    for(r in 1:length(dfe[,1])) dfe1[r,] <- if(dfe[r,]>0) evar1[dfe[r,]] else 'int'
  }
  df1 <- df1[-which(names(df1)=='scenario')]
  df1 <- data.frame(dfe1, round(df1,2) )
  if(any(names(df1)=='model')) df1$model[df1$model<0] <- 'int'
  df1
}


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
