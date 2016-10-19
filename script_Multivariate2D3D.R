# From 2D PCA to 3D PTA : introduction to multivariate analysis
# author: Romain Frelat
# date: 20 october 2016

## A. Getting ready ---------------------------------------

#A.1 Load package
library(ade4)
library(PTAk)
#If there is an error message, type : install.packages(c("ade4", "PTAk")

#A.2 Load the dataset
load("IBTS_Tensor.Rdata")
dim(IBTS_tensor)
dimnames(IBTS_tensor)

#A.3 Understanding the variables
#Array
dim(IBTS_tensor)
IBTS_tensor[18,14,6]#Select one element, e.g. Abundance of Cod, in 1998, in RA6
IBTS_tensor[18,,6]#Select one vector, e.g. abundance of Cod in 1998
IBTS_tensor[18,,]#Select one matrix, e.g. abundance of Cod

#List
names_tensor <- dimnames(IBTS_tensor) # the list of names is stored in a new variable
length(names_tensor) #there are three elements in the list, one element for each dimension
names_tensor[[2]] #show the second element of the list
names_tensor[[1]][18]#show the 18th element of the first element of the list

# Your turn: 
#1. What is the index of Hake (Merluccius merluccius) in the dataset ?
#2. What is the abundance of Hake (Merluccius merluccius) in 1988 in RA 1 ?
#3. What is the abundance of Hake between 2010 and 2015 in RA 1?
#4. Can you show the evolution of Hake abundance between 1985 and 2015 in RA 1?

# Solution:
#1
which(names_tensor[[1]]=="Merluccius merluccius")
#2
IBTS_tensor[33,4,1]
#3
IBTS_tensor[33,26:31,1]
#4
plot(names_tensor[[2]], IBTS_tensor[33,,1], type="l", 
     xlab="Time", ylab="CPUE (n/h)", main="Abundance of Hake in RA1")

## B. Two dimension : Introduction to PCA -----------------

#B.1 Preparing the dataset

# From 3D to 2D
IBTS_space <- apply(IBTS_tensor,c(3,1),mean)
dim(IBTS_space) #the new matrix has 7 rows (areas) and 65 columns (species)

# Checking the distribution of the data
#boxplot() is used to look at the distribution
boxplot(as.vector(IBTS_space), main="raw CPUE")
#The CPUE is very skewed, one can not see the difference between the 1st quarter, the median and the 3rd quarter
#So data should be log transformed
IBTS_logspace <- log(IBTS_space+1)
#The new distribution of the log transformed CPUE
boxplot(as.vector(IBTS_logspace), main="log CPUE")

# Scaling the data
par(mfrow=c(2,1), mar=c(2,4,3,1))
boxplot(IBTS_logspace, main="Abundance", ylab="CPUE in log", 
        xaxt="n")
mtext("Species", side = 1, line = 0)
boxplot(scale(IBTS_logspace), main="Anomaly (= abundance normalized per species)", 
        ylab="Anomaly", xaxt="n")
mtext("Species", side = 1, line = 0)

#B.2 Principal Component Analysis

# Run a PCA and choose the correct number of PC.
pca_space=dudi.pca(IBTS_logspace, scale = TRUE, center = TRUE)
#**Select the number of axes: **
# Please type 2

#To see how much variance the axes explain:
inertia.dudi(pca_space)$TOT

#B.3 Interpretation of the PC.
pca_space$li # or pca_space$co
par(mfrow=c(1,2), mar=c(0,0,0,0))
#Show the weight of the variables:
s.label(pca_space$li, xax=1, yax=2)
s.label(pca_space$co, xax=1, yax=2, clabel = 0.4)

#B.4 Clustering the species

#Compute the distance between species
dist_species=dist(pca_space$co, method = "euclidean")

#Build a tree with Ward method
den=hclust(dist_species,method = "ward.D2")

#Plot the dendogram
par(mar=c(2,3,3,1))
plot(den, hang=-1, ax = T, ann=T, xlab="", sub="", cex=0.6)

#Choosing the number of cluster
nclust<-5

#Visualize the cutting
rect.hclust(den, k=nclust, border="dimgrey")

#Create the clusters
clust_space <- as.factor(cutree(den, k=nclust))

#Visualize the cluster in the PC axis
s.class(pca_space$co,fac=clust_space, col=rainbow(nclust),xax=1,yax=2)

# Your turn: 
#1. How would you interpret the clusters 4 and 5 ?
#2. How many species are located mainly in the south-west of the North Sea (RA 5, i.e. grouped in cluster 1) ?
#3. In which cluster is grouped Saithe (*Pollachius virens*)?

# Solution
#2.
table(clust_space) #There are 14 species in cluster 1
#3.
clust_space[names_tensor[[1]]=="Pollachius virens"] #Saithe is in cluster 4

##C. Three dimension : Tensor Decomposition --------------- 

#C1. Preparing the dataset

# Checking the distribution of the data
#boxplot() is used to look at the distribution
boxplot(IBTS_tensor, main="raw CPUE")
#Data should be log transformed
IBTS_logtensor <- log(IBTS_tensor+1)
#The new distribution of the log tranformed CPUE
boxplot(IBTS_logtensor, main="log CPUE")

par(mfrow=c(1,2), mar=c(3,3,3,1))
boxplot(as.vector(IBTS_tensor), main="raw CPUE")
boxplot(as.vector(IBTS_logtensor), main="log CPUE")

#C2. Scaling the data
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

#C3. Run a PTA and choosing the number of PT.
# Run a PTA 
pta<-PTA3(IBTS_logscale, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta,testvar = 0)

# Choosing the number of PT.
#Create the scree plot
out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct<-(pta[[3]]$pct*pta[[3]]$ssX/pta[[3]]$ssX[1])[out]
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance")

#C4. Interpretation of PT.

par(mfrow=c(1,2))
plot(pta, mod=c(2,3), nb1 = 1, nb2 = 11, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 11, lengthlabels = 3)

par(mfrow=c(2,2))
plot(pta, mod=c(2,3), nb1 = 1, nb2 = 6, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 6, xpd=NA, lengthlabels = 4)
plot(pta, mod=c(2,3), nb1 = 1, nb2 = 7, xpd=NA, lengthlabels = 4)
plot(pta, mod=1, nb1 = 1, nb2 = 7, xpd=NA, lengthlabels = 4)

#C5. Clustering
#Create the matrix with the projection of species on the 4 PT
keep <- c(1, 6, 7, 11) # PT that are kept in the analysis
coo<-t(pta[[1]]$v[c(keep),])
labkeep <- paste0(pta[[3]]$vsnam[keep], " - ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%")

#Compute the distance between species
dist1=dist(coo, method = "euclidean")

#Build a tree with Ward linkage
den=hclust(dist1,method = "ward.D2")

#Plot the dendogram
par(mar=c(1,3,1,1))
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="",labels = FALSE)

#Choose the number of clusters
nclust<-6

#Visualize the cutting
rect.hclust(den, k=nclust, border=rainbow(nclust)[c(6,5,2,4,3,1)])

#Create the clusters
clust_3D <- as.factor(cutree(den, k=nclust))

#Visualize them 
par(mfrow=c(1,3))
s.class(coo, fac = clust_3D, xax=1, yax=2, col=rainbow(nclust), clabel = 2)
text(min(coo[,1])-0.03,0, labkeep[2], srt=90, xpd=NA, cex=1.5)
text(0,max(coo[,2])+0.02, labkeep[1], xpd=NA, cex=1.5)
s.class(coo, fac = clust_3D, xax=1, yax=3, col=rainbow(nclust), clabel = 2)
text(min(coo[,1])-0.03,0, labkeep[3], srt=90, xpd=NA, cex=1.5)
text(0,max(coo[,3])+0.02, labkeep[1], xpd=NA, cex=1.5)
s.class(coo, fac = clust_3D, xax=1, yax=4, col=rainbow(nclust), clabel = 2)
text(min(coo[,1])-0.03,0, labkeep[4], srt=90, xpd=NA, cex=1.5)
text(0,max(coo[,4])+0.02, labkeep[1], xpd=NA, cex=1.5)
