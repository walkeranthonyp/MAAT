################################
#
# MAAT Leaf Model - useful ploting functions
# 
# AWalker (walkerap@ornl.gov) 
# September 2017
#
################################

library(viridis)
library(lattice)
library(latticeExtra)
library(stringr)
# library(hexbin)
library(fmsb)
# library(gamlss)



# radar plots
############################################

# radar plot
ppSA_radar <- function(df1, alphap=100, max_si=0.6,
                       title='Test Radar SA', printleg=T, vnames=par_names, lnames=NULL,
                       ... ) {
  
  # can pass arguments like:
  # - centerzero=T; which puts 0 at the center of the plot
  # - vlcex       ; size of axis labels
  
  # add max and min to data
  n   <- dim(df1)
  df1 <- rbind(max=rep(max_si,n[2]) , min=rep(0,n[2]) , df1)
  
  # convert to data frame
  df1 <- as.data.frame(df1)
  
  # print(df1)
  
  # sort colour scheme
  colors_border <- viridis(n[1])
  col_poly_rgb  <- col2rgb(colors_border,alpha=T)
  colors_poly   <- rgb(t(col_poly_rgb[1:3,]),alpha=alphap, maxColorValue=255)
  
  # get parameter names
  if(!is.null(vnames)) {
    ss_names <- match(names(df1),names(vnames))
    labs     <- vnames[ss_names]
  } else labs <- names(df1)

  # plot radar
  radarchart( df1, axistype=0, maxmin=T,
              title=title,
              # custom polygon
              pcol=colors_border,
              pfcol=colors_poly,
              plwd=3, plty=1,
              # custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.4,0.1), cglwd=0.5,
              # custom labels
              vlcex=0.9,
              vlabels=labs
  )
  
  # plot legend
  if(is.null(lnames)) lnames <- rownames(df1[-c(1,2),])
  if(printleg) legend(x=-1.3, y=1.4, legend = lnames, bty = "n", 
                      pch=20 , col=colors_border , text.col = "grey", cex=1, pt.cex=1)
}


# organise radar plots
radar_plot <- function(df1, var1='type', var2=NULL, lnames=NULL,
                       index='scenario', values='sensitivity1', par='variable',
                       add_space=0, ... ) {
  
  # below variables break df1 into subsets to plot on a single radar plots
  # determine if a first subsetting variable has been given
  if(!is.null(var1)) {
    var1_c <- which(names(df1)==var1)
    var1_l <- unique(df1[,var1_c])
    print(var1_l)
  } else {
    var1_c <- NULL
    var1_l <- -999
  }
  
  # determine if a second subsetting variable has been given
  if(!is.null(var2)) {
    var2_c <- which(names(df1)==var2)
  } else {
    var2_c <- NULL
    var2_l <- -999
  }

  # rename variable
  vss <- which(names(df1)==values)
  names(df1)[vss] <- 'values' 
  pss <- which(names(df1)==par)
  names(df1)[pss] <- 'par' 
  
  # labels
  if(is.null(lnames)) lnames <- rownames(df1[-c(1,2),])

  # define plot layout
  lo <- if(var2_l[1]==-999) c(ceiling(length(var1_l)/3),3) else c(length(var1_l),length(var2_l))
  # add space for legend & additional plots
  lo[1] <- lo[1] + 1 + add_space
  # define plot layouts, plot in columns first, then rows
  print(lo)
  lmat <- 
    if(add_space > 0) t(matrix(1:prod(lo),ncol=lo[1]))[c(3:(2+add_space),1,2),]
    else              t(matrix(1:prod(lo),ncol=lo[1]))
  layout(lmat, widths=rep(2,prod(lo)), heights=c(rep(2,prod(lo)-lo[2]),rep(1,lo[2])), TRUE )
  par(mar=c(0,1,0,1), oma=c(0,0,0,0) )  
  
  # subsetting loops
  for(v1 in var1_l) {
  # for(v1 in 1) {
    if(v1!=-999) df1_sub <- subset(df1,df1[,var1_c]==v1) 
    else         df1_sub <- df1
    
    # reget var2 labels in case there are different numbers of uniques values of var2 within each value of var1
    if(!is.null(var2)) var2_l <- unique(df1_sub[,var2_c])
    for(v2 in var2_l) {

      if(v2!=-999) df1_sub2 <- subset(df1_sub,df1_sub[,var2_c]==v2) 
      else         df1_sub2 <- df1_sub
      
      iss <- 
        if(df1_sub2[1,var1] == 'Combined' )      which(names(df1_sub2)=='scenario') 
        else if(df1_sub2[1,var1] == 'Scenario' ) which(names(df1_sub2)=='scenario') 
        else if(df1_sub2[1,var1] == 'Model' )    which(names(df1_sub2)=='model') 
        else                                     which(names(df1_sub2)==index)
      names(df1_sub2)[iss] <- 'ind' 
      
      # print(df1_sub2)
      # df1_sub2$ind <- as.character(df1_sub2$ind)
      df1_sub_tab  <- xtabs(values~ind+par, df1_sub2 )
      # print(df1_sub_tab)
      
      # radar plot
      title <- if(lo[2] == 1) NULL else if(is.null(var2)) v1 else paste(v1,v2) 
      ppSA_radar(df1_sub_tab, title=paste('\n',title), printleg=F, ... )
    }      
  }
  
  # plot legend
  for( l in lnames ) {
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE, ann=FALSE )

    if(!is.null(l)) {
      n    <- length(l)
      ncol <- if(n>12) ceiling(n/12) else 1
      colors_border <- viridis(n[1])
      legend(x=1, y=9, legend = l, bty = "n", ncol=ncol, 
             pch=20 , col=colors_border , text.col = "grey", cex=1, pt.cex=1)
    }
  }
}



# plot sensitivity with 4 or less variables
############################################

plot_SA_4orless <- function(yv, xv, cv1=NULL, cv2=NULL, gv=NULL,
                            gls = seq(0,25,5), scales_x = NULL,
                            ... ) {  
  
  form <- 
    if(is.null(cv1))      (yv ~ xv) 
    else if(is.null(cv2)) (yv ~ xv | cv1)
    else                  (yv ~ xv | cv1 * cv2)
  
  if(is.numeric(xv)) {

    xyplot(form, groups=gv, ...,
           type='b',lwd=2,
           layout=c(3,1), scales=list(tck=c(-0.5,0), alternating=F),
           par.settings = simpleTheme(col=col,pch=19),
           strip=strip.custom(bg='grey90',fg='grey90',var.name=expression(I),strip.levels=c(T,T)),
           panel=function(...){
             panel.abline(h=gls,col='gray90')
             panel.xyplot(...)
           },
           auto.key=list(x=0.05,y=0.85,points=F,lines=T)
    )

  } else {

    bwplot(form, groups=gv, ...,
           layout=c(3,1), scales=list(tck=c(-0.5,0), alternating=F),
           par.settings = simpleTheme(col=col,pch=19,lwd=2),
           strip=strip.custom(bg='grey90',fg='grey90',var.name=expression(I),strip.levels=c(T,T)),
           panel=function(x,y,...){
             panel.abline(h=gls,col='gray90')
             panel.xyplot(x,y,...,type='b')
           },
           auto.key=list(x=0.05,y=0.85,points=F,lines=T)
    )
  }
}


# data variance plots - violin plots
############################################

# plot violin plots of variables
plot_SA_modvar <- function(yv, xv, cv1=NULL, cv2=NULL, pv=panel.violin,
                           lout = c(3,1), bxw = c(100,0.1), gls = seq(0,25,5),
                           scales_x = NULL,
                           xlim = NULL, ylim= NULL, xlab=NULL, ylab=NULL, ... ) {
  
  form <- 
    if(is.null(cv1))      (yv ~ xv) 
    else if(is.null(cv2)) (yv ~ xv | cv1)
    else                  (yv ~ xv | cv1 * cv2)
  
  scales_arg <- 
    if(!is.null(scales_x)) list(tck=c(-0.5,0),alternating=F,x=scales_x) 
    else                   list(tck=c(-0.5,0),alternating=F) 
  
  if(is.numeric(xv)) {
    
    xyplot(form ,
           xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
           scales=scales_arg, layout=lout, 
           strip=strip.custom(fg='grey90', bg='grey90', par.strip.text=list(cex=0.8), var.name=expression(I), strip.levels=c(T,T) ),
           panel=function(..., box.width=bxw ){
             panel.abline(h=gls, col='gray90' )
             pv(..., horizontal=F, col = "lightblue", varwidth = FALSE, box.width = box.width[1] )
             panel.bwplot(..., horizontal=F, col = 'black', coef=100, cex=0.8, pch='|', fill='gray', box.width = prod(box.width) )
           },
           par.settings = list(box.rectangle=list(col='black'), box.umbrella=list(alpha=0,col='black',lty=1),
                               plot.symbol = list(pch='.', cex = 0.1))
    )
    
  } else {
    
    bxw[1] <- bxw[1] / 100
    bwplot(form ,
           # xlim=xlim, 
           ylim=ylim, xlab=xlab, ylab=ylab,
           scales=scales_arg, layout=lout, 
           strip=strip.custom(fg='grey90', bg='grey90', par.strip.text=list(cex=0.8), var.name=expression(I), strip.levels=c(T,T) ),
           panel=function(..., box.width=bxw ){
             panel.abline(h=gls, col='gray90' )
             # pv(..., horizontal=F, col = "lightblue", varwidth = FALSE, box.width = box.width[1] )
             # panel.bwplot(..., horizontal=F, col = 'black', coef=100, cex=0.8, pch='|', fill='gray', box.width = prod(box.width) )
             pv(..., col = "lightblue", varwidth = FALSE, box.width = box.width[1] )
             panel.bwplot(..., col = 'black', coef=100, cex=0.8, pch='|', fill='gray', box.width = prod(box.width) )
           },
           par.settings = list(box.rectangle=list(col='black'), box.umbrella=list(alpha=0,col='black',lty=1),
                               plot.symbol = list(pch='.', cex = 0.1))
    )
    
  }
}


# violin plots of model output at different environmental conditions
# - broken out by process represenation and (if requested) condvar1  
plot_env_array <- function(a1, yvar='A', condvar1='lim', ..., esubs=NULL, lout=c(1,1), bxw=c(1,0.1), gls=seq(0,25,5) ) {
  
  yvar_col  <- which(dimnames(a1)[[4]]==yvar) 
  cond1_col <- which(dimnames(a1)[[4]]==condvar1) 
  
  # create x-axis and conditioning vectors 
  np <- 2
  sub_fnames <- function(i, fnames ) {
    sep <- c(0, str_locate(fnames,', '), nchar(fnames)+1 )
    substr(fnames, sep[i]+1, sep[i+1]-1 )
  }
  ex_fnames <- function(fnames) vapply(seq(1,2*np,2), sub_fnames, character(1) , fnames=fnames )
  fo <- t(vapply(dimnames(a1)[[1]], ex_fnames, character(np) ))
  
  # create xaxis labels
  labs <- c(paste(evar1name,evar2name,sep='\n'), paste(evar1,evar2,sep='\n') )
  
  # function names conditioning vector
  f1  <- rep(fo[,1], prod(dim(a1)[2:3]) )
  f2  <- rep(fo[,2], prod(dim(a1)[2:3]) )

  # environment vector
  envint <- 1:dim(a1)[2]
  e      <- rep(rep(envint, each=dim(a1)[1] ), dim(a1)[3] )

  # make plots
  # variability against environmental conditions
  p1 <- plot_SA_modvar(yv=a1[,,,yvar], xv=e, xlab='', scales_x=list(at=c(0,envint),labels=labs), ... )

  # variability against environmental conditions broken out by model?
  if(dim(fo)[1]>5) {
    p2 <- plot_SA_modvar(yv=a1[,,,yvar], xv=e, cv1=f2, xlab='', scales_x=list(at=c(0,envint),labels=labs), ... )

  } else {
    p2 <- plot_SA_modvar(yv=a1[,,,yvar], xv=e, cv1=f1, cv2=f2, xlab='', scales_x=list(at=c(0,envint),labels=labs), ... )
    p2 <-
      useOuterStrips(p2,
                     strip=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)),
                     strip.left=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)))
  }
  
  # variability against environmental conditions broken out by first process in SA and condvar 
  if(!is.null(condvar1)) {
    p3 <- plot_SA_modvar(yv=a1[,,,yvar], xv=e, cv1=as.factor(a1[,,,cond1_col]), cv2=f1, xlab='', scales_x=list(at=c(0,envint),labels=labs), ... )
    p3 <-
      useOuterStrips(p3,
                     strip=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)),
                     strip.left=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)))
    
    return(list(p1,p2,p3))
    
  } else {
    
    return(list(p1,p2))
  }

}



# model output vs parameter plots
############################################

# plot data density on bivariate plots using hexbins (bivariate histograms)
plot_parDens <- function(ov, pv, cond1=NULL, cond2=NULL, gls, ... ) {

  pclos <- if(is.null(cond2)) {
    if(is.null(cond1)) ov ~ pv
    else             { ov ~ pv | as.factor(cond1) }              
  } else               ov ~ pv | as.factor(cond1) * as.factor(cond2)
  
  p1 <-
    hexbinplot(
      pclos,
      trans=sqrt, inv=function(x) x^2,
      cex.labels=0.5,cex.title=0.5,
      border = F,colorkey=F,aspect=1,
      ... ,
      colramp=viridis,
      type='l',lwd=2,
      as.table=T ,
      scales=list(cex=0.8,alternating=F,tck=c(-1,0)),
      panel=function(...){
        panel.abline(h=gls,col='gray90')
        panel.hexbinplot(...)
      })
  
  if(!is.null(cond1)&!is.null(cond2)) {
    p1 <- 
      useOuterStrips(p1,
                     strip=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)),
                     strip.left=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)))
  }
  
  rm(ov,pv,cond1,cond2)
  p1
}


# quantile regressions of model output vs parameter values using gamlss - SLOW!
plot_parSens_list <- function(list1,yvar='A',par_col=1,...) {
  
  yvar_col <- which(colnames(list1[[1]]$out_mat)==yvar) 
  
  # values
  ov <- vapply(list1,function(l) l$out_mat[,yvar_col], FUN.VALUE=list1[[1]]$out_mat[,1] )
  p  <- vapply(list1,function(l) l$parsB[,par_col],    FUN.VALUE=list1[[1]]$parsB[,1]   ) 
  f1 <- vapply(list1,function(l) l$fnames,             FUN.VALUE=list1[[1]]$fnames      ) 
  f2 <- vapply(list1,function(l) l$fnamesB,            FUN.VALUE=list1[[1]]$fnamesB     )
  e  <- vapply(list1,function(l) l$env,                FUN.VALUE=list1[[1]]$env         )
  
  # key
  
  p1 <-
    xyplot(
      as.vector(ov) ~ as.vector(p) | as.vector(f1) * as.vector(f2) , 
      groups = as.vector(e) ,
      ... ,
      xlab=colnames(list1[[1]]$parsB)[par_col] ,
      type='l',col=viridis(3),lwd=2,border=F ,
      as.table=T ,
      scales=list(cex=0.8,alternating=F,tck=c(-1,0)) ,
      panel=panel.superpose ,
      panel.groups=function(x,y,...){
        gdata    <- data.frame(x=x,y=y)
        glss_mod <- gamlss(y~cs(x),data=gdata)
        mm       <- c(min(x),max(x))
        cent     <- centiles.pred(glss_mod,xname='x',cent=c(5,50,95),xvalues=seq(mm[1],mm[2],(mm[2]-mm[1])/100),data=gdata ) 
        panel.polygon(x=c(cent$x,rev(cent$x)),y=c(cent$C5,rev(cent$C95)),alpha=0.4,...)
        panel.lines(cent$C50~cent$x,...)
      }
    )
  
  p1 <- 
    useOuterStrips(p1,
                   strip=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)),
                   strip.left=strip.custom(bg='grey90',par.strip.text=list(cex=0.8)))
  
  rm(ov,p,f1,f2,e)
  p1
}



### END ###