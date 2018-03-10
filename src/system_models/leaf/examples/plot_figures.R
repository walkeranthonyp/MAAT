##################################
#
# Process MAAT output and plot for GMD paper
#
# AWalker
# Nov 2017
#
##################################

rm(list=ls())
library(lattice)
library(latticeExtra)
library(stringr)
library(grid)
library('viridis')


wd    <- '#PROJECTDIR#'
ddate <- '#OUTPUTDATE#' 
wdd   <- paste0(wd,'/results/',ddate)

ifile <- 'out.csv'
an_id <- 'AnvNum'
tc_id <- 'Tcurves'


# read data
setwd(wdd)
df1   <- read.csv(paste(an_id,ifile,sep='_'))
head(df1)

df2   <- read.csv(paste(tc_id,ifile,sep='_'))
head(df2)
df2$Tres <- df2$vcmaxlt / df2$vcmax


# Figure labels & parameters
fs     <- 10
ylab   <- list(expression('A ['*mu*mol*' '*m^-2*s^-1*']'), fontsize=fs ) 
xlabCa <- list(expression(C[a]*' [Pa]'), fontsize=fs )
xlabCi <- list(expression(C[i]*' [Pa]'), fontsize=fs )
ylabT  <- list('T scalar [-]' , fontsize=fs )
xlabT  <- list(expression('Leaf T ['^o*'C]'), fontsize=fs )

width  <- 5.5118
col    <- c('blue','red','black')
col2   <- rev( viridis(4)[2:4] )
cex    <- 1
lwd    <- 0.2

deftpar <- trellis.par.get()
gmdtpar <- deftpar
gmdtpar$axis.text$cex   <- cex
gmdtpar$fontsize$text   <- 6
gmdtpar$fontsize$points <- 3


### Analytical vs. numerical
# A/Ca curves
p1 <- 
  xyplot(A~ca|leaf.rs*as.factor(leaf.g0), df1, groups=leaf.solver,
         par.settings=simpleTheme(pch=c(20,20,4), col=col, cex=cex, lwd=lwd ),
         ylab=ylab, xlab=xlabCa, scales=list(alternating=F, tck=c(-0.5,0)),  
         as.table=T, auto.key=list(text=c('Analytical','Analytical Quad','Numerical'), columns=3 ),
         panel=panel.superpose,
         panel.groups=function(group.number, ... ) {        
           type <- c('b','b','p')
           if(group.number==1) panel.abline(h=0, col='grey80')
           panel.xyplot(t=type[group.number],...) 
         })

p1so <- useOuterStrips(p1,
                       strip =      strip.custom(bg='grey90', factor.levels=sapply(str_split(levels(df1$leaf.rs),'_'), function(v) v[length(v)] )),
                       strip.left = strip.custom(bg='grey90', factor.levels=c(expression(g[0]=='0.00' ), expression(g[0]=='0.01'))) )
p1so

# A/Ci curves
p2 <- 
  xyplot(A~ci|leaf.rs*as.factor(leaf.g0), df1, groups=leaf.solver, subset=leaf.solver!='f_A_r_leaf_analytical',
         par.settings=simpleTheme(pch=c(20,20,4), cex=cex, col=col, lwd=lwd ),
         ylab=ylab, xlab=xlabCi, scales=list(alternating=F, tck=c(-0.5,0)),  
         as.table=T,
         panel=panel.superpose,
         panel.groups=function(group.number, subscripts, ... ) {        
           type <- c('b','b','p')
           if(group.number==2) {
             panel.abline(h=0, col='grey80' )
             panel.abline(v=df1$transition[subscripts][1], col='grey80' )
           }
           panel.xyplot(t=type[group.number],...) 
         })

p2so <- useOuterStrips(p2,
                       strip      = strip.custom(bg='grey90', factor.levels=sapply(str_split(levels(df1$leaf.rs),'_'), function(v) v[length(v)] )),
                       strip.left = strip.custom(bg='grey90', factor.levels=c(expression(g[0]=='0.00' ), expression(g[0]=='0.01'))) )
p2so

# write Figure
rx      <- c(0.0)
ry      <- c(0.48, 0.0)
pwidth  <- 1.0
pheight <- c(0.52, 0.50)

pdf('figAC.pdf', width=width, height=width*0.9 ) 

trellis.par.set(gmdtpar)

pushViewport( viewport(y=unit(ry[1],'npc'), x=unit(rx[1],'npc'), width=pwidth, height=pheight[1], just=c('left','bottom')) )
print(p1so,newpage=F)
popViewport()

pushViewport( viewport(y=unit(ry[2],'npc'), x=unit(rx[1],'npc'), width=pwidth, height=pheight[2], just=c('left','bottom')) )
print(p2so,newpage=F)
popViewport()

pushViewport( viewport(y=unit(0.95,'npc'), x=unit(0.00,'npc'), width=0.05, height=0.05, just=c('left','bottom')) )
grid.text('(a)', gp=gpar(fontsize=fs) )
popViewport()

pushViewport( viewport(y=unit(0.45,'npc'), x=unit(0.00,'npc'), width=0.05, height=0.05, just=c('left','bottom')) )
grid.text('(b)', gp=gpar(fontsize=fs) )
popViewport()

dev.off()


# Stomatal limitation Figure
df1s    <- subset(df1,leaf.solver=='f_R_Brent_solver'&leaf.g0==0.01)
df1s_gs <- rbind(data.frame(A=df1s$A,leaf.rs=df1s$leaf.rs,ca=df1s$ca,gs='Stomatal Limitation'), data.frame(A=df1s$A_noR,leaf.rs=df1s$leaf.rs,ca=df1s$ca,gs='No Stomatal Limitation') )
p3 <- 
  xyplot(A~ca|leaf.rs, df1s_gs, groups=gs,
         par.settings=simpleTheme(pch=c(4,20), cex=cex, col=col[c(3,1)] ),
         ylab=ylab, xlab=xlabCa, scales=list(alternating=F, tck=c(-0.5,0)),
         layout=c(5,1), auto.key=list(columns=2),
         strip=strip.custom(bg='grey90', factor.levels=sapply(str_split(levels(df1$leaf.rs),'_'), function(v) v[length(v)] )),
         panel=function(...) {
           panel.abline(h=0, col='grey80')
           panel.xyplot(...)
         })
p3

# write Figure
pdf('fig_gslim.pdf', width=width, height=width/3); trellis.par.set(gmdtpar); p3; dev.off()



### T response Figure

txlim <- c(-3,49)
df2s1 <- subset(df2,leaf.vcmax_tcor_asc!='f_scalar_none'&leaf.vcmax_tcor_des=='f_scalar_none')
df2s2 <- subset(df2,leaf.vcmax_tcor_asc=='f_scalar_none'&leaf.vcmax_tcor_des!='f_scalar_none')
df2s3 <- subset(df2,leaf.vcmax_tcor_asc!='f_scalar_none'&leaf.vcmax_tcor_des!='f_scalar_none')
df2s1$leaf.vcmax_tcor_asc <- as.factor(as.character(df2s1$leaf.vcmax_tcor_asc))
df2s2$leaf.vcmax_tcor_des <- as.factor(as.character(df2s2$leaf.vcmax_tcor_des))
df2s3$leaf.vcmax_tcor_asc <- as.factor(as.character(df2s3$leaf.vcmax_tcor_asc))
df2s3$leaf.vcmax_tcor_des <- as.factor(as.character(df2s3$leaf.vcmax_tcor_des))

des_labs <- sapply(str_split(levels(df2$leaf.vcmax_tcor_des)[2:4],'_'), function(v) v[length(v)-1] )

p4 <- 
  xyplot(Tres~leaf.temp|leaf.vcmax_tcor_des, df2s1, groups=leaf.vcmax_tcor_asc, 
         par.settings=simpleTheme(pch=c(0,20), cex=cex, col=col[1:2] ),
         auto.key=list(text=c('Arrhenius',expression(Q[10])), columns=2),
         ylab=ylabT, xlab=NULL, scales=list(alternating=F, tck=c(-0.5,0)), layout=c(1,1),
         xlim=txlim,ylim=c(0.1,4),
         strip=strip.custom(bg='grey90', factor.levels= 'ascending' ),
         panel=function(...) {
           panel.abline(h=0:1, col='grey80' )
           panel.abline(v=25,  col='grey80' )
           panel.xyplot(...)
         })
p4

p5 <- 
  xyplot(Tres~leaf.temp|leaf.vcmax_tcor_asc, df2s2, groups=leaf.vcmax_tcor_des, 
         par.settings=simpleTheme(pch=c(15,1,2), cex=cex, col=col2 ),
         ylab=NULL, xlab=xlabT, scales=list(alternating=F, tck=c(-0.5,0)), layout=c(1,1),
         xlim=txlim,ylim=c(0,1.35),
         strip=strip.custom(bg='grey90', factor.levels= 'descending' ),
         panel=function(...) {
           panel.abline(h=0:1, col='grey80' )
           panel.abline(v=25,  col='grey80' )
           panel.xyplot(...)
         })
p5

p6 <- 
  xyplot(Tres~leaf.temp|leaf.vcmax_tcor_asc, df2s3, groups=leaf.vcmax_tcor_des, 
         par.settings=simpleTheme(pch=c(15,1,2), cex=cex, col=col2 ),
         auto.key=list(text=c(des_labs), columns=3),
         ylab=NULL, xlab=NULL, scales=list(alternating=F, tck=c(-0.5,0)), layout=c(2,1),
         xlim=txlim,ylim=c(0,1.35),
         strip=strip.custom(bg='grey90', factor.levels=c('Arrhenius',expression(Q[10])) ),
         panel=function(...) {
           panel.abline(h=0:1, col='grey80' )
           panel.abline(v=c(25,32),  col='grey80' )
           panel.xyplot(...)
         })
p6


# write Figure
rx      <- c(0.0,0.26,0.52)
ry      <- c(0.08,0.0,0.08)
pwidth  <- c(0.3,0.3,0.5)
pheight <- c(0.9,0.93,0.9)

pdf('figTres.pdf', width=width, height=width/3 ) 
trellis.par.set(gmdtpar)

pushViewport( viewport(y=unit(ry[1],'npc'), x=unit(rx[1],'npc'), width=pwidth[1], height=pheight[1], just=c('left','bottom')) )
print(p4,newpage=F)
popViewport()

pushViewport( viewport(y=unit(ry[2],'npc'), x=unit(rx[2],'npc'), width=pwidth[2], height=pheight[2], just=c('left','bottom')) )
print(p5,newpage=F)
popViewport()

pushViewport( viewport(y=unit(ry[3],'npc'), x=unit(rx[3],'npc'), width=pwidth[3], height=pheight[3], just=c('left','bottom')) )
print(p6,newpage=F)
popViewport()

pushViewport( viewport(y=unit(0.92,'npc'), x=unit(rx[1],'npc'), width=0.05, height=0.05, just=c('left','bottom')) )
grid.text('(a)', gp=gpar(fontsize=fs) )
popViewport()

pushViewport( viewport(y=unit(0.92,'npc'), x=unit(rx[2],'npc'), width=0.05, height=0.05, just=c('left','bottom')) )
grid.text('(b)', gp=gpar(fontsize=fs) )
popViewport()

pushViewport( viewport(y=unit(0.92,'npc'), x=unit(rx[3],'npc'), width=0.05, height=0.05, just=c('left','bottom')) )
grid.text('(c)', gp=gpar(fontsize=fs) )
popViewport()

dev.off()









