#Master model to analyse the relationship between various explanatory variables (derived from the PCA)
#and the log-transformed inverse estimated richness
#Physicochemical princomp -> p3
#Pesticide princomp -> p 
#Data is as close to normal as possible 
#lmer is just not working
#No signal in the data between anything (look at lm) - this is all you need to talk about for your report

#Package installation and loading 
install.packages('paran')
library(lme4)
library(paran)
library(nlme)

#Merging PCA scores from two datasets 
masterpca <- cbind(p, p3)
masterpca

#Using Horn's parallel analysis - removing the variance that is not important 
#Horn's parallel analysis performed on the scaled data sets 'pesticide2' and 'chemical2' to determine 
#the significant principal components
#This analysis helps in identifying the meaningful variance in the data
#Horn's parallel analysis, the 'paran()' function, compares the eigenvalues obtained from the actual data to those 
#obtained from randomly generated data sets of the same size. By comparing the actual eigenvalues to the random
#eigen values, Horn's parallel analysis helps to identify the significant principal components
paran(scale(log(t(pesticide2))))
paran(scale(t(chemical2)))
#The result of Horn's parallel analysis for component retention on the 'pesticide2' and 'chemical2' datasets suggest that there 
#are two and four components that should be retained, respectively. 
#These components likely represent underlying patterns or structures in the data that explain significant variance. 
#expl_variable dataframe created combining specific columns from 'p' and 'p3' 
expl_variable <- cbind(p[,1:2], p3[,1:4]) 

#Fixed effects
#Assigning column names and data preparation
colnames(expl_variable) = c('pt1', 'pt2', 'ch1', 'ch2', 'ch3', 'ch4')
head(expl_variable)
something = na.omit(data.frame(inv_est, expl_variable, location))

#A LMM is constructed using the explanatory variables 'pt1', 'pt2', etc. and a random effect for 'location'
#The summary and the calculated pesudo-R squared are provided to assess the model fit
mastermodel <- lmer(log(inv_est) ~ pt1 + pt2 + ch1 + ch2 + ch3 + ch4 + (1 | location)
                    , data = something)
summary(mastermodel)
#Master linear mixed-effects model summary output
#The random effect 'location' has extremely small estimated variance 
#None of the fixed effects are statistically significant as indicated by the p-values
#The linear mixed-effects model did not find statistically significant relationships between the predictor variables 
#('pt1', 'pt2', 'ch1', 'ch2', 'ch3', 'ch4') and the response variable (log-transformed inverse estimated crustacean richness)
r.squaredGLMM(mastermodel)
ranef(mastermodel)
sum(chemical2 == 0)
sum(pesticide2 < 0)

#Q-Q plot
#Extract the residuals
residuals <- residuals(mastermodel)
#Create a Q-Q plot
png("qqplot.png")
qqnorm(residuals)
qqline(residuals)
dev.off()

#Save the mastermodel summary output 
summary_text <- capture.output(summary(mastermodel))
writeLines(summary_text, "summary_output.txt")

#Histograms of pesticide data are plotted for visual inspection 
for (i in 1:nrow(pesticide2)) { 
  hist(as.numeric(pesticide2[i,]))
  readline()
  }
i

#Factor analysis on the 'expl_variable' data, aiming to identify underlying latent factors 
#that explain correlation structure between observed variables 
fact <- factanal(expl_variable, factors = 2, scores = "regression")$scores
fact

plot(p3[,1], p3[,2])

#Non-linear mixed-effects models the relationship between the log-transformed inverse estimated richness 
#and the explanatory variables, considering the location as a random effect
fm <- nlme(log(inv_est) ~ expl_variable + location
           , data = data.frame(inv_est, expl_variable, location)
           , fixed = expl_variable ~ 1, random = location ~ 1, start = c(expl_variable = 1))

#Histogram of the log-transformed inverse estimated richness
hist(log(inv_est))

rand(mastermodel)
#ANOVA-like table for the random-effects in the linear mixed-effects model 'mastermodel' fitted to the data
#Pr(>Chisq) shows the p-value is 1, suggesting that the random effect 'location' is not significant, as its 
#inclusion does not significantly improve model fit. 
sort(inv_est)

#Fit a linear regression model
summary(lm(log(inv_est) ~ pt1 + pt2 + ch1 + ch2 + ch3 + ch4 + location, data = something))
#Summary suggests that the linear regression model does not provide a good fit to the data, as indicated 
#by the low R^2 value, and the non-significant F-statistic. 
#Additionally, none of the individual predictor variables appear to be statistically significant
#, as indicated by their p-values. 

chemical2
pesticide2
