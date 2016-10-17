# From 2D PCA to 3D PTA : introduction to multivariate analysis
# author: Romain Frelat
# date: 20 october 2016

## A. Getting ready ---------------------------------------

#A.1 Load package
library(ade4)
library(PTAk)
library(RColorBrewer)
source("Myheatmap.R") #load a homemade function to plot a heatmap

#A.2 Load the dataset
load("IBTS_Tensor.Rdata")
dim(IBTS_tensor)
dimnames(IBTS_tensor)

#A.3 Understanding the variables
#List
length(dimnames(IBTS_tensor))
dimnames(IBTS_tensor)[[2]] #show the second element of the list

#Array
dim(IBTS_tensor)
#Select one element, e.g. Abundance of Cod, in 1998, in RA6
IBTS_tensor[18,14,6]
#Select one vector, e.g. abundance of Cod in 1998
IBTS_tensor[18,,6]
#Select one matrix, e.g. abundance of Cod
IBTS_tensor[18,,]

# Your turn : 
# 1. What is the index of Hake (Merluccius merluccius) in the dataset ?
# 2. What is the abondance of Hake (Merluccius merluccius) in 1988 in RA 1 ?
# 3. What is the abondance of Hake between 2010 and 2015 in RA 1?
# 4. Can you show the evolution of Hake abondance between 1985 and 2015 in RA 1?

#Solution
#1
which(dimnames(IBTS_tensor)[[1]]=="Merluccius merluccius")
#2
IBTS_tensor[33,4,1]
#3
IBTS_tensor[33,26:31,1]
#4
plot(dimnames(IBTS_tensor)[[2]], IBTS_tensor[33,,1], type="l", 
     xlab="Time", ylab="CPUE (n/h)", main="Abundance of Hake in RA1")

## B. Two dimension : Introduction to PCA -----------------

#B.1 Tranformation of the data from 3D to 2D
IBTS_space <- apply(IBTS_tensor,c(3,1),mean)
dim(IBTS_space)

#B.2 Checking the distribution of the data
#boxplot() is used to look at the distribution
boxplot(as.vector(IBTS_space), main="raw CPUE")
#The CPUE is very skewed, so data should be log transformed
IBTS_logspace <- log(IBTS_space+1)
#The new distribution of the log tranformed CPUE
boxplot(as.vector(IBTS_logspace), main="log CPUE")

par(mfrow=c(1,2), mar=c(3,3,3,1))
boxplot(as.vector(IBTS_space), main="raw CPUE")
boxplot(as.vector(IBTS_logspace), main="log CPUE")

#B.3
#PCA normalize your data
par(mfrow=c(1,2), mar=c(4,4,3,1))
boxplot(IBTS_logspace, main="", ylab="CPUE in log", 
        xaxt="n", xlab="species")
boxplot(scale(IBTS_logspace), main="normalized per species", 
        ylab="Anomaly", xaxt="n", xlab="species")
myHeatmap(t(IBTS_logspace), mary=3, cex.y=0.5, title="log CPUE",pal="BrBG")
myHeatmap(t(scale(IBTS_logspace)), mary=3, cex.y=0.5, title="scaled CPUE",pal="BrBG")


#B.4 Run the PCA and choosing the number of PC.
pca_space=dudi.pca(IBTS_logspace, scale = TRUE, center = TRUE)
#**Select the number of axes: **
#Please type 2

#To see how much variance the axes explain :
inertia.dudi(pca_space)$TOT

#B.5 Interpretation of the PC.
par(mfrow=c(1,2))
#Show the weight of the variables :
s.corcircle(pca_space$co, clabel = 0.4, ="Species") #on PC 1 and 2
s.label(pca_space$li, xax=1, yax=2) #on PC 1 and 2

#B.6 Clustering the species
#Compute the distance between species from their projection on PC
dist_species=dist(pca_space$co, method = "euclidean")

#Hierarchical clustering with Ward method
den=hclust(dist_species,method = "ward.D2")

#plot the dendogram
par(mfrow=c(1,1))
plot(den, hang=-1, ax = T, ann=T, xlab="", sub="")

#Choosing the number of cluster
nclust<-5

#Visualize the cutting
rect.hclust(den, k=nclust, border="dimgrey")

#Create the clusters
clusters <- as.factor(cutree(den, k=nclust))

#Visaualize the cluster in the PC axis
par(mfrow=c(1,2))
s.class(pca_space$co,fac=clusters, col=rainbow(nclust),xax=1,yax=2)

##C. Three dimension : Tensor Decomposition ---------------
##C.1 Transforming the data
boxplot(IBTS_tensor, main="raw CPUE")
#Data should be log transformed
IBTS_logtensor <- log(IBTS_tensor+1)
#The new distribution of the log tranformed CPUE
boxplot(IBTS_logtensor, main="log CPUE")
```

```{r, echo=FALSE}
par(mfrow=c(1,2), mar=c(3,3,3,1))
boxplot(as.vector(IBTS_tensor), main="raw CPUE")
boxplot(as.vector(IBTS_logtensor), main="log CPUE")
```

#### Normalize your data
Contrary to PCA, the normalization of a tensor is not straight forward, and have to be done manually before running a PTA.

```{r, collapse=TRUE, comment =""}
#Scaling per species
#Create a new empty array
IBTS_logscale<-array(0,dim=dim(IBTS_tensor))
#Loop scanning each species
for (i in 1:dim(IBTS_tensor)[1]){
  #Calculating the mean and sd of the log CPUE for species i
  ma<-mean(IBTS_logtensor[i,,])
  sa<-sd(IBTS_logtensor[i,,])
  #Saving the anomaly in the array
  IBTS_logscale[i,,]<-(IBTS_logtensor[i,,]-ma)/sa
}
#Copy the labels to the new array
dimnames(IBTS_logscale)<-dimnames(IBTS_tensor)
```

#### Run the PTA and choosing the number of PC.
The PCA is run with the function `dudi.pca`. The data is normalized so the options `scale` and `center` are set to `TRUE`. The function is interactive, showing you the plot of the variance explained by the successive Principal Components (PC)

```{r, comment =""}
PTA<-PTA3(IBTS_logscale, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(PTA,testvar = 0)
```

```{r, comment =""}
#Create the scree plot
out <- !substr(PTA[[3]]$vsnam, 1, 1) == "*"
gct<-(PTA[[3]]$pct*PTA[[3]]$ssX/PTA[[3]]$ssX[1])[out]
barplot(sort(gct, decreasing = TRUE), names.arg = order(gct, decreasing = TRUE),
        ylab="Percentage of variance", xlab="PC")
```


```{r, comment =""}
plot(PTA, mod=c(2,3), nb1 = 1, nb2 = 11)
```

```{r, comment =""}
keep <- c(1, 6, 7, 11)

coo<-t(PTA[[1]]$v[c(keep),])
#labkeep <- paste0(PTA[[3]]$vsnam[keep], " - ", round((100 * (PTA[[3]]$d[keep])^2)/PTA[[3]]$ssX[1],1), "%")
labkeep <- paste0(paste0("PT", 1:4), " - ", round((100 * (PTA[[3]]$d[keep])^2)/PTA[[3]]$ssX[1],1), "%")


par(mar = c(3,3,3,1))
for (i in seq_along(keep)){
  lab <- labkeep[i]
  temp <- PTA[[3]]$v[keep[i],] %o% PTA[[2]]$v[keep[i],]
  dimnames(temp) <- list(dimnames(IBTS_tensor)[[3]], dimnames(IBTS_tensor)[[2]])
  myHeatmap(temp, pal="BrBG", title=lab, colscale = FALSE, mary = -3)
}
```

```{r, comment =""}
dist1=dist(coo, method = "euclidean")
den=hclust(dist1,method = "ward.D2")
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="",labels = FALSE)
nclust<-6
clust.3D <- as.factor(cutree(den, k=nclust))
par(mar=c(1,3,1,1))
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="",labels = FALSE)
rect.hclust(den, k=nclust, border=rainbow(nclust)[c(6,5,2,4,3,1)])
s.class(coo, fac = clust.3D, xax=1, yax=2, col=rainbow(nclust), clabel = 2)
text(min(coo[,1])-0.03,0, labkeep[2], srt=90, xpd=NA, cex=1.5)
text(0,max(coo[,2])+0.02, labkeep[1], xpd=NA, cex=1.5)
s.class(coo, fac = clust.3D, xax=1, yax=3, col=rainbow(nclust), clabel = 2)
text(min(coo[,1])-0.05,0, labkeep[3], srt=90, xpd=NA, cex=1.5)
text(0,max(coo[,3])+0.02, labkeep[1], xpd=NA, cex=1.5)
s.class(coo, fac = clust.3D, xax=1, yax=4, col=rainbow(nclust), clabel = 2)
text(min(coo[,1])-0.03,0, labkeep[4], srt=90, xpd=NA, cex=1.5)
text(0,max(coo[,4])+0.02, labkeep[1], xpd=NA, cex=1.5)
```





####Preview:
png("A.png", width = 1500, height = 1500, res=300)
myHeatmap(t(scale(IBTS_logspace)), mary=3, cex.y=0.5, pal="BrBG", 
          title="Average CPUE per area")
dev.off()
png("B.png", width = 1500, height = 1500, res=300)
layout(matrix(c(1,2,1,3), ncol=2))
par(mar=c(2,2,2,2))
plot(den, hang=-1, ax = T, ann=T, xlab="", sub="",labels = FALSE)
rect.hclust(den, k=nclust, border="dimgrey")
s.label(pca_space$li, xax=1, yax=2) #on PC 1 and 2
s.class(pca_space$co,fac=clusters, col=rainbow(nclust),xax=1,yax=2)
dev.off()