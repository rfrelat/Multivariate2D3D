## Define a function for plotting a matrix ----------------
#x : the matrix
#horiz : add a horizontal line 
#pal : name of palette used from RColorBrewer
myHeatmap <- function(x, horiz=NULL, verti =NULL, pal="BuGn", colscale=TRUE, 
                      mary=0, col.l="black", breaks=NULL, labcol= "", 
                      rm.empty=FALSE, cex.x=0.7, cex.y=0.7, 
                      tck.x=-0.03, tck.y=-0.05, padj.x = -0.5, hadj.y=1, ...){
  require(RColorBrewer)
  min <- min(x, na.rm=TRUE)
  max <- max(x, na.rm=TRUE)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  if(colscale) {
    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  }
  # Set color palette
  if (length(pal)==1){
    if (is.null(breaks)){
      ColorRamp<- brewer.pal(9,pal)
      ColorLevels <- seq(min, max, length=10)
    } else {
      ColorRamp<- brewer.pal(length(breaks)-1,pal)
      ColorLevels <- breaks
    }
  } else {
    if (is.null(breaks)){
      ColorRamp<- pal
      ColorLevels <- seq(min, max, length=length(pal)+1)
    } else {
      ColorRamp<- pal
      ColorLevels <- breaks
    }
  }
  
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(2,5+mary,1.5,1))
  #par(mar = c(0.1,0.1,0.1,0.1))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max), breaks=ColorLevels)
  
  if( !is.null(title) ){
    title(main=title)
  }
  if (!is.null(horiz) ){
    abline(h=horiz, lwd=3, col=col.l)
  }
  if (!is.null(verti) ){
    abline(v=verti, lwd=2)
  }
  if (rm.empty){
    axis(BELOW<-1, at=(1:length(xLabels))[xLabels!=""], 
         labels=xLabels[xLabels!=""], cex.axis=cex.x, tck=tck.x, padj = padj.x)
    axis(LEFT <-2, at=(1:length(yLabels))[yLabels!=""], 
         labels=yLabels[yLabels!=""], las= 1, cex.axis= cex.y, 
         tck=tck.y, hadj = hadj.y)
  } else {
    axis(BELOW<-1,  at=1:length(xLabels), labels=xLabels, 
         cex.axis=cex.x, tck=tck.x, padj = padj.x)
    axis(LEFT <-2,  at=1:length(yLabels), labels=yLabels, 
         las= 1, cex.axis=cex.y, tck=tck.y, hadj = hadj.y)
  }
  
  
  # Color Scale
  if (colscale) {
    par(mar = c(3,2.5,2.5,2))
    nColorLevels<-ColorLevels[-1]-(diff(ColorLevels)/2)
    image(1, nColorLevels,
          matrix(data=nColorLevels, ncol=length(nColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    axis(2, at = median(nColorLevels), labels = labcol)
  }
  if(colscale) {
    layout(1)
  } else {
    return(ColorLevels)
  }
}


colscale <- function(dat,col, met=c("quant", "range", "abs"), ...){
  nbreaks <- length(col)
  if (met=="quant"){
    colcut<-quantile(dat, probs = seq(0, 1, length.out = nbreaks+1), ...)
  }
  if (met=="range"){
    colcut<-seq(min(dat, ...), max(dat, ...), length.out = nbreaks+1)
  }
  if (met=="abs"){
    colcut<-seq(-max(abs(dat), ...), max(abs(dat), ...), length.out = nbreaks+1)
  } 
  colord<-cut(dat, breaks = colcut, include.lowest = TRUE)
  res <- as.character(col[colord])
  return(res)
}
