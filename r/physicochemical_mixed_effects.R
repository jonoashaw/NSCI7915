#Physicochemical mixed effects

install.packages('lme4')
install.packages('lmerTest')
install.packages('MuMIn')

library(lme4)
library(lmerTest)
library(MuMIn)

#JA's code to determine the correct multiplier 
x = as.matrix(read.delim('16S_crust_watersamples.txt'))

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
source('sadrad.R')

#Defining a function
source('inverse_x_distribution.R')

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
names(est) = colnames(x)

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

#Import physicochemical data
chemical <- read.delim('Natural_stressors_environmental_data_TCP20.txt',row.names = 1)
chemical

#Figure out how many columns and rows 
dim(chemical)

#Lining up column names with OTU data and physicochemical data 
chemical2 <- chemical[rownames(chemical) %in% colnames(z),]
dim(chemical2)
rownames(chemical)
colnames(chemical)
chemical2 = chemical2[, -1]
chemical2 = t(chemical2)
chemical2

#Collapse physicochemical data into less variables 
p2 <- princomp(t(chemical2))$loadings
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