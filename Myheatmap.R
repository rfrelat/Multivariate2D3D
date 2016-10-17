# Myheatmap.R
# Author : Romain Frelat
# Last modified : 16.10.2016
# Inspired from http://www.phaget4.org/R/image_matrix.html

# Define a function for plotting a 2D matrix
# Require RColorBrewer package
# x : the matrix
# pal : name of palette used from RColorBrewer
# colscale : if it adds a color scale on the side (TRUE), or not (FALSE)
# breaks : define the breaks for the color scale. If NULL, the breaks are defined automatically
# mary : define the extra margin on the left side
# horiz, verti : add a horizontal/vertical line to the plot
# col.l : color of the added line
# labcol : define the labels in the x axis
# rm.empty : remove the empty labels
# cex.x, cex.y : size of the labels in x and y-axis
myHeatmap <- function(x, pal="BuGn", colscale=TRUE, breaks=NULL, mary=0,
                      horiz=NULL, verti =NULL, col.l="black", labcol= "", 
                      rm.empty=FALSE, cex.x=0.7, cex.y=0.7, ...){
  require(RColorBrewer)
  min <- min(x)
  max <- max(x)
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
         labels=xLabels[xLabels!=""], cex.axis=cex.x)
    axis(LEFT <-2, at=(1:length(yLabels))[yLabels!=""], 
         labels=yLabels[yLabels!=""], las= 1, cex.axis= cex.y)
  } else {
    axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=cex.x)
    axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, 
         las= 1, cex.axis=cex.y)
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
