#Physicochemical mixed effects

install.packages('lme4')
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
  k = isce(n)
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
dim(x)
length(inv_est)
inv_est

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
chemical <- read.delim('data/data:raw:/natural_stressors_environmental_data_tcp20.txt',row.names = 1)
chemical
chemical = chemical[,-c(9,23:33)]

#Figure out how many columns and rows 
dim(chemical)

#Lining up column names with OTU data and physicochemical data 
z
rownames(chemical)
n = colnames(z)
for (i in 1:length(n))
  substr(n[i],3,3) = '-'
colnames(z) = n
chemical2 <- chemical[rownames(chemical) %in% colnames(z),]
dim(chemical2)
rownames(chemical)
colnames(chemical)
chemical2 = chemical2[, -1]
chemical2 = t(chemical2)
chemical2

#Collapse physicochemical data into less variables 
p3 <- princomp(t(chemical2))$scores
princomp(t(chemical2))$loadings
#Water conductivity - how are crustaceans affected? 
#This is salinity - look at Aashi's recently published paper which shows that ->
#fish and shark communities are influenced by salinity in the same estuaries 
summary(glm(log(inv_est) ~ p3[,1:5]))

#glm with location as a variable
#Diversity of MA is really low 
summary(glm(log(inv_est) ~ location))
chemical2[6,]

#Mixed-effects
summary(lmer(log(inv_est) ~ p3[,1:5] + (1|location)))
r.squaredGLMM(lmer(log(inv_est) ~ p3[,1:5] + (1|location)))

#An idea of what distribution we are dealing with 
for (i in 1:ncol(x))	{
  readline()
plot(rev(sort(x[x[,i] > 0,i])), log = 'xy',type = 'l')
} 


