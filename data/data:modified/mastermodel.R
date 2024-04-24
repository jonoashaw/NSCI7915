#Master model
#physicochemical princomp -> p3
#pesticide princomp -> p 
#data is as close to normal as possible 
#lmer is just not working
#no signal in the data between anything (look at lm) - this is all you need to talk about for your report

install.packages('paran')

library(lme4)
library(paran)
library(nlme)

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
something = na.omit(data.frame(inv_est, expl_variable, location))
mastermodel <- lmer(log(inv_est) ~ pt1 + pt2 + ch1 + ch2 + ch3 + ch4 + (1 | location)
                    , data = something)
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

fm <- nlme(log(inv_est) ~ expl_variable + location
           , data = data.frame(inv_est, expl_variable, location)
           , fixed = expl_variable ~ 1, random = location ~ 1, start = c(expl_variable = 1))

hist(log(inv_est))

rand(mastermodel)
sort(inv_est)

summary(lm(log(inv_est) ~ pt1 + pt2 + ch1 + ch2 + ch3 + ch4 + location, data = something))
