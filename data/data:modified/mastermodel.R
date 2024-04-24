#Master model
#physicochemical princomp -> p3
#pesticide princomp -> p 

install.packages('paran')

library(lme4)
library(paran)

#merge
masterpca <- cbind(p, p3)
masterpca

#Using Horn's parallel analysis - removing the variance that is not important 
paran(scale(log(t(pesticide2))))
paran(scale(t(chemical2)))
expl_variable <- cbind(p[,1:2], p3[,1:4]) 

#fixed effects
colnames(expl_variable) = c('pt1', 'pt2', 'ch1', 'ch2', 'ch3', 'ch4')
head(expl_variable)
mastermodel <- lmer(log(inv_est) ~ expl_variable + (1 | location), data = data.frame(inv_est, expl_variable, location))
summary(mastermodel)
r.squaredGLMM(mastermodel)
ranef(mastermodel)
sum(chemical2 == 0)
sum(pesticide2 < 0)

for (i in 1:nrow(pesticide2)) { 
  hist(as.numeric(pesticide2[i,]))
  readline()
  }
i

fact <- factanal(expl_variable, factors = 2, scores = "regression")$scores
fact

plot(p3[,1], p3[,2])
nlme 