#Conducting various analyses on crustacean community data. 
#Specifically, 16S sequencing data from water samples. 
#Pesticide mixed effects

#Package installing and loading. 
install.packages('lme4')
install.packages('lmerTest')
install.packages('MuMIn')
library(lme4)
library(lmerTest)
library(MuMIn)

#Reading and processing the data stored in '16S_crust_watersamples.txt'
#JA's code to determine the correct multiplier 
x = as.matrix(read.delim('data/data:modified/16S_crust_watersamples.txt'))
x
x = na.omit(matrix(as.numeric(x),nrow(x),ncol(x),dimnames=list(rownames(x),colnames(x))))

y = x %% 1

m = array()

for (i in 1:ncol(y))	{
  for (j in 1:100)	{
    z = (y[,i] * j) %% 1
    z[z > 1 - 1e-4 & z < 1 + 1e-4] = 0
    if (sum(z) < 1e-4)	{
      m[i] = j
      cat(i,j,'\n')
      break
    }
  }
}

#Column 9 looks like 24, and 30 like 18
m[9] = 24
m[30] = 18

z = round(t(t(x) * m))

table(z)

#sadrad
source('r/sadrad.R')

#Defining a function
source('r/inverse_x_distribution.R')

#Richness
?richness
?halfpower
?ipower
?sodds
?logd
?par
#Estimates richness of crustacean communities using halfpower, ipower, sodds, and logd. 
#Calculates AICc values for each richness estimation. 
est = array()
inv_est = array()

dim(z)
a = matrix(NA, 30, 4)

source("r/r:functions:/ipower.R")
for (i in 1:ncol(z))	{
  n = z[,i]
  n = n[n > 0]
  #k = isce(n)
  #	plot(k$fitted.curve,type='l')
  #	points(k$subsampled.richness,cex=0.3,col='red')
  #est[i] = k$asymptotic.richness
  #inv_est[i] = ipower(n)$richness
  inv_est[i] = sodds(n)$richness
  a[i,1] = halfpower(n)$AICc
  a[i,2] = ipower(n)$AICc
  a[i,3] = sodds(n)$AICc
  a[i,4] = logd(n)$AICc
  #	readline()
}
#Calculating fit statistics
#'inv_est[i] = sodds(n)$richness' estimates richness using sodds()
#a[i,] lines calculate AICc for different richness estimation methods and stores them in matrix 'a'

colnames(a)=c("HP", "IP", "SO", "LS")

#Rank abundance distribution plot 
#Plots the rank abundance distribution of the crustacean communities based on the fit statistics 
#Ordination of fit statistics 
a
a[,1]
a[a > 1e4] = NA
#PCA of fit statistics stored in 'a'
plot(princomp(na.omit(a))$scores)
b = array()
for (i in which(! is.na(a[,1]))){ 
  b[i]=colnames(a)[which(a[i,]==min(a[i,]))]}
b

b = na.omit(b)
#Plots the PCA scores 'ps' with labelled points representing the fit statistics
ps = princomp(na.omit(a))$scores
plot(ps, cex = 0)
text(ps, labels = b, cex = exp(ps[,3])^0.05, col=hsv(h=as.numeric(as.factor(b))/4))
dim(ps)
length(b)
exp(ps[,3])^0.05

plot(ps, cex = 0)
text(ps, labels = b, cex = 1 - 0 * exp(ps[,3])^0.05, col=hsv(h=as.numeric(as.factor(b))/4))
#Each point in the plot represents a richness estimation method, and its position is determined by its scores 
#on the principal components. 
#Each point on the plot corresponds to one of the richness estimation methods. The labels near each point typically
#indicate the method used for richness estimation (e.g. HP, IP, SO, and LS)
#Points that are closer together represent richness estimation methods that have similar fit statistics across the different estimation methods 
#Points that are farther apart represent methods with dissimilar fit statistics 
#Methods that cluster together may indicate similarities or differences in their performance across different 
#richness estimation methods. 
#Methods that are outliers or distant from the main cluster may indicate unique or extreme performance characteristics. 

#What distributions are working for the dataset
est[est < 0] = NA
names(inv_est) = colnames(x)

#n = colnames(x)
#for (i in 1:length(n))
#substr(n[i],3,3) = '-'
#names(est) = n

#crust_12S_estimates = est
#Save(crust_12S_estimates,file='crust_12S_estimates')

location = array()

for (i in 1:30)
  location[i] = strsplit(colnames(x)[i],'_')[[1]][1]
location = as.factor(location)

#Fitting a generalised linear model to examine the relationship between the log-transformed richness estimates and location. 
summary(glm(log(inv_est) ~ location))
#The model indicates that location might have a weak effect on log-transformed richness estimates

#Import pesticide data
#Load pesticide data and align with crustacean community data
pesticide <- read.delim('data/data:raw:/pesticide_water_tcp20_without_detection_threshold.txt',row.names = 1)
pesticide

#Figure out how many columns and rows 
dim(pesticide)

#Lining up column names with OTU data and pesticide data 
pesticide2 <- pesticide[rownames(pesticide) %in% colnames(z),]
dim(pesticide2)
rownames(pesticide)
colnames(pesticide)
pesticide2 = pesticide2[, -1]
pesticide2 = t(pesticide2)
pesticide2

#Collapse pesticide data into less variables 
#Principle component analysis on the pesticide data and fits a linear model between the first principal
#component and the log-transformed inverse estimated richness of crustaceans
p <- princomp(scale(log(t(pesticide2))))$loadings
summary(glm(log(inv_est) ~ p[,1:1]))
p

#Mixed-effects
#Fit a mixed-effects model considering the pesticide data along with the first five principal components
#as fixed effects and the location as a random effect 
summary(lmer(log(inv_est) ~ p[,1:5] + (1|location)))

#Evaluates the model using r.squaredGLMM
r.squaredGLMM(lmer(log(inv_est) ~ p[,1:5] + (1|location)))
p[,1:5]
location
is.factor(location)

#Making a dataframe
#Combine inverse estimated richness, principal component scores, and location in a dataframe
#Fit another mixed-effects model
#Exploratory analyses and plotting
xy = cbind(inv_est,p,location)
xy = data.frame(xy)
head(xy)
summary(lmer(log(inv_est) ~ Comp.1 + (1|location),data = xy))
r.squaredGLMM(lmer(log(inv_est) ~ p[,1:5] + (1|location)))
summary(lmer(log(inv_est) ~p[,1:5] + (1|location),data = xy))
p
colnames(p)
dim(p)
cbind(rownames(p),rownames(p3))

#Plotting rank abundance distributions for each of the crustacean communities represented 
#in the dataset 'z'
#Overlays fitted rank abundance distributions using ipower and logd to compare their fits to
#the observed data 
#Series of rank-abundance distribution plots 
for (i in 1:30){
n = z[,i]
n = n[n > 0]
n = rev(sort(n))
plot(n, type = 'l', bty = 'l', log = 'y', xlab = 'Rank of Count', ylab = "Count of Reads")
lines(rev(ipower(n)$fitted.RAD),col = 'violet')
lines(rev(logd(n)$fitted.RAD),col = 'cyan')
legend('topright',c('data', 'IP', 'LS'),col = c('black', 'violet', 'cyan'),lwd = 1, cex = 0.8)
}
#The x-axis represents the rank of each count (species)
#The y-axis represents the count of reads, plotted on a logarithmic scale
#Fitted models
#ipower model(violet line) represents the fitted RAD and attempts to fit a power law to the 
#observed counts 
#logd model (cyan line) represents the fitted RAD and attempts to fit a log-normal distribution 
#to the observed counts 
#Distribution shape 
#The black line (observed counts) typically shows a steep drop, indicating that a few 
#species are very abundance while most are rare
#Model fit
#The closer the ipower and logd models follow the black line, the better the models fit the observed data
#If the lines diverge significantly from the black line, it suggests that the model might not be suitable 
#for describing the observed distribution
#Comparison between models
#The fit of the two models can be compared to determine which better describes the observed data
#A better fit indicates a more appropriate model for the underlying ecological process
#Good fit
#If either the violet or the cyan line follows closely the black line across most ranks, it suggests that
#the corresponding model provides a good description of the observed abundance pattern 
#Poor fit 
#If both lines deviate significantly from the black line, it suggests that neither model adequately captures 
#the distribution of counts. This might prompt consideration of alternative models or a reassessment of the data. 
#Differences between models
#If one model fits better in the high-rank region (most abundant species) and another in the low-rank region (rare species), it suggests that different
#processes might govern the abundance of common and rare species 

 

