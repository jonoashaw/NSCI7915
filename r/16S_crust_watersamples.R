# test data

x = matrix(sample(1:20,100,replace=T),10,10)
x = t(t(x) / 1:10)


x = as.matrix(read.delim('/Users/alroy/Data/16S_crust_watersamples.txt'))

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

library(richness)

for (i in 1:ncol(z))	{
  n = z[,i]
  n = rev(sort(n[n > 0]))
  plot(n,type='l',log='y')
  lines(rev(logd(n)$fitted.RAD),col='red')
  #	lines(rev(pln(n)$fitted.RAD),col='orange')
  lines(rev(sodds(n)$fitted.RAD),col='dodgerblue')
  lines(rev(inv(n)$fitted.RAD),col='green3')
  #	lines(rev(nbin(n)$fitted.RAD),col='green3')
  readline()
}


est = array()
inv_est = array()

for (i in 1:ncol(z))	{
  n = z[,i]
  n = n[n > 0]
  k = isce(n)
  #	plot(k$fitted.curve,type='l')
  #	points(k$subsampled.richness,cex=0.3,col='red')
  est[i] = k$asymptotic.richness
  inv_est[i] = inv(n)$richness
  #	readline()
}

est[est < 0] = NA
names(est) = colnames(x)

n = colnames(x)
for (i in 1:length(n))
  substr(n[i],3,3) = '-'
names(est) = n

crust_12S_estimates = est
save(crust_12S_estimates,file='crust_12S_estimates')

location = array()

for (i in 1:30)
  location[i] = strsplit(colnames(x)[i],'_')[[1]][1]

summary(glm(log(inv_est) ~ location))



