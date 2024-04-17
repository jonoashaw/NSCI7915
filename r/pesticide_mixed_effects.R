#Pesticide mixed effects

install.packages('lme4')
install.packages('lmerTest')
install.packages('MuMIn')

library(lme4)
library(lmerTest)
library(MuMIn)

#JA's code to determine the correct multiplier 
x = as.matrix(read.delim('data/data:modified/16S_crust_watersamples.txt'))

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

# column 9 looks like 24, and 30 like 18

m[9] = 24
m[30] = 18

z = round(t(t(x) * m))

table(z)

#sadrad
source('r/sadrad.R')

#Defining a function
source('r/inverse_x_distribution.R')

#Richness
est = array()
inv_est = array()

for (i in 1:ncol(z))	{
  n = z[,i]
  n = n[n > 0]
  #k = isce(n)
  #	plot(k$fitted.curve,type='l')
  #	points(k$subsampled.richness,cex=0.3,col='red')
  #est[i] = k$asymptotic.richness
  inv_est[i] = inv(n)$richness
  #	readline()
}

est[est < 0] = NA
names(inv_est) = colnames(x)

#n = colnames(x)
#for (i in 1:length(n))
#substr(n[i],3,3) = '-'
#names(est) = n

#crust_12S_estimates = est
#save(crust_12S_estimates,file='crust_12S_estimates')

location = array()

for (i in 1:30)
  location[i] = strsplit(colnames(x)[i],'_')[[1]][1]
location = as.factor(location)

summary(glm(log(inv_est) ~ location))

#Import pesticide data
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
pesticde2 = t(pesticde2)
pesticide2

#Collapse pesticide data into less variables 
p <- princomp(t(pesticide2))$loadings
summary(glm(log(inv_est) ~ p[,1:1]))

#Mixed-effects
summary(lmer(log(inv_est) ~ p[,1:5] + (1|location)))
r.squaredGLMM(lmer(log(inv_est) ~ p[,1:5] + (1|location)))
p[,1:5]
location
is.factor(location)

#Making a dataframe
xy = cbind(inv_est,p,location)
xy = data.frame(xy)
head(xy)
summary(lmer(log(inv_est) ~ Comp.1 + (1|location),data = xy))
r.squaredGLMM(lmer(log(inv_est) ~ p[,1:5] + (1|location)))
summary(lmer(log(inv_est) ~p[,1:5] + (1|location),data = xy))
p
