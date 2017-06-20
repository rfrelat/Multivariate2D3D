# From 2D PCA to 3D PTA : introduction to multivariate analysis
# Graphical functions
# author: Romain Frelat
# date: 28 april 2017

source("script_Multivariate2D3D.R")
source("tools.R")
nkeep <- length(keep)

par(mfrow=c(2,2))
#Map of time x area
par(mar = c(3,3,3,1), oma=c(1,1,0,0))
for (i in seq_along(keep)){
  lab <- labkeep[i]
  temp <- pta[[3]]$v[keep[i],] %o% pta[[2]]$v[keep[i],]
  dimnames(temp) <- list(dimnames(IBTS_logscale)[[3]], dimnames(IBTS_logscale)[[2]])
  myHeatmap(temp, pal="BrBG", title=lab, colscale = FALSE, 
            mary = -3,cex.x = 1, cex.y = 1, rm.empty = TRUE,
            tck.x = -0.04, tck.y = -0.04, padj.x = -0.5, hadj.y = 0.5) 
  mtext("Roundfish areas", side = 2, line = 2, xpd=NA, cex=0.8)
}

#Species weight on each PT
par(mar = c(3,3,2,1))
nSp <- dimnames(IBTS_logscale)[[1]]
cooSp <- pta[[1]]$v
limx <- range(cooSp[keep,])
for (i in seq_along(keep)){
  dotchart(sort(cooSp[keep[i],]), xlim=limx, pch=16, cex=0.5,
           labels = nSp[order(cooSp[keep[i],])])
  abline(v=0, lty=2, col="grey")
}
mtext("Weight PC", side = 1, line = 3, xpd=NA, cex=0.7)

#Temporal time series
year <- dimnames(IBTS_logscale)[[2]]
cooTi <- pta[[2]]$v
limy <- range(cooTi[keep,])
for (i in seq_along(keep)){
par(mar=c(3,3,2,0))
plot(year, pta[[2]]$v[keep[i],], type="l", xaxt="n", yaxt="n", ylim=limy)
axis(1, tck=-0.04,  padj = -0.5)
axis(2, tck=-0.04, las=1, hadj=0.7)
mtext("Year", side = 1, line = 2, xpd=NA, cex=0.7)
mtext("Anomaly", side = 2, line = 2, xpd=NA, cex=0.7)
abline(h=0, lty=2, col="grey")
}

#Spatial map
lat <- rnorm(7,mean = 55, sd=2)
lon <- rnorm(7,mean =5, sd=1)
cooSp <- pta[[3]]$v[keep,]
colpal<-brewer.pal(5,"BrBG")
colo <- matrix(colscale(dat = cooSp, col = colpal, met="quant"), ncol = nkeep)
require(maps) #make nicer maps
require(mapdata) #add country background
require(Cairo)
X11(type="cairo")

for (i in seq_along(keep)){
map("worldHires", xlim=c(-10, 10), ylim=c(50, 60), col="gray90", fill=TRUE, border="gray70")
points(lon,lat, pch=21, cex=2, bg=colo[,i], col="grey50")
map.axes()
title(labkeep[i])
}
