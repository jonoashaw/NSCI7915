#Physicochemical mixed effects
#Conducting analyses on crustacean community data in the context of physicochemical data
#Package installation and loading 
install.packages('lme4')
install.packages('nlme')
install.packages('lmerTest')
install.packages('MuMIn')
install.packages('bbmle')
install.packages('poilog')
install.packages('sads')
install.packages('r/r:functions:/richness-2.3', repos = NULL, type = "source")
install.packages("minpack.lm")
library(lme4)
library(lmerTest)
library(MuMIn)
library(richness)
library(minpack.lm)

#JA's code to determine the correct multiplier 
x2 = as.matrix(read.delim('data/data:modified/16S_crust_watersamples.txt'))

#Error replaces unknown data with NA's
x2 = na.omit(matrix(as.numeric(x2),nrow(x2),ncol(x2),dimnames=list(rownames(x2),colnames(x2))))

y2 = x %% 1

m2 = array()

for (i in 1:ncol(y2))	{
  for (j in 1:100)	{
    z2 = (y2[,i] * j) %% 1
    z2[z2 > 1 - 1e-4 & z2 < 1 + 1e-4] = 0
    if (sum(z2) < 1e-4)	{
      m2[i] = j
      cat(i,j,'\n')
      break
    }
  }
}

#Column 9 looks like 24, and 30 like 18

m2[9] = 24
m2[30] = 18

z2 = round(t(t(x2) * m2))

table(z2)

#sadrad
source('r/sadrad.R')

#Defining a function
source('r/inverse_x_distribution.R')

#Richness estimation 
est = array()
inv_est = array()

for (i in 1:ncol(z))	{
  n2 = z2[,i]
  n2 = n2[n2 > 0]
  k2 = isce(n2)
	#plot(k$fitted.curve,type='l')
  #points(k$subsampled.richness,cex=0.3,col='red')
  est[i] = k$asymptotic.richness
  inv_est[i] = inv(n)$richness
  #readline()
}
est
est[est < 0] = NA
plot(est,inv_est)
cor.test(est,inv_est,method = "s")
abline(0,1)
names(inv_est) = colnames(x)
dim(x2)
length(inv_est)
inv_est

#n = colnames(x)
#for (i in 1:length(n))
#substr(n[i],3,3) = '-'
#names(est) = n

#crust_12S_estimates = est
#save(crust_12S_estimates,file='crust_12S_estimates')

#Create factor variable 'location' based on the sample names
location = array()

for (i in 1:30)
  location[i] = strsplit(colnames(x2)[i],'_')[[1]][1]
location = as.factor(location)

#Fit generalised linear model to analyse the relationship between the log-transformed inverse estimated richness and location
summary(glm(log(inv_est) ~ location))

#Import physicochemical data
chemical <- read.delim('data/data:raw:/natural_stressors_environmental_data_tcp20.txt',row.names = 1)
chemical

#Removing uneccessary coordinate data 
cbind(colnames(chemical), 1:ncol(chemical))
chemical = chemical[,-c(9,20:21,23:33)]
chemical

#Figure out how many columns and rows 
dim(chemical)

#Lining up column names with OTU data and physicochemical data 
z2
rownames(chemical)
n = colnames(z2)
for (i in 1:length(n))
  substr(n[i],3,3) = '-'
colnames(z2) = n
chemical2 <- chemical[rownames(chemical) %in% colnames(z2),]
dim(chemical2)
rownames(chemical)
colnames(chemical)
chemical2 = chemical2[, -1]
chemical2 = t(chemical2)
chemical2

#Collapse physicochemical data into less variables 
#Conduct principal component analysis on the scaled physicochemical data and retain the scores
p3 <- princomp(scale(t(chemical2)))$scores
p3

#Fit generalised linear model to analyse the relationship between the log-transformed inverse estimated richness
#and the first five principal components of the physicochemical data 
#Water conductivity/Salinity - how are crustaceans affected? 
#This is salinity - look at Parikh et al. 2024 which shows similar results
summary(glm(log(inv_est) ~ p3[,1:5]))
#None of the first five principal components are significantly associated with 'log(inv_est)'

#GLM with location as a variable
#Diversity of MA is really low 
summary(glm(log(inv_est) ~ location))
chemical2[6,]

#Fit mixed-effects model to analyse the relationship between the log-transformed inverse estimated richness 
#and the first five principal components of the physicochemical data, with location as a random effect
summary(lmer(log(inv_est) ~ p3[,1:5] + (1|location)))
#Some variability in 'log(inv_est)' between different locations, but the fixed effects do not explain a significant
#portion of the variance

#Calculate marginal and conditional R-squared values
r.squaredGLMM(lmer(log(inv_est) ~ p3[,1:5] + (1|location)))

#An idea of what distribution we are dealing with 
for (i in 1:ncol(x))	{
  readline()
plot(rev(sort(x2[x2[,i] > 0,i])), log = 'xy',type = 'l')
} 
p3
colnames(p3)
rownames(p3)
dim(p3)
str(p3)

