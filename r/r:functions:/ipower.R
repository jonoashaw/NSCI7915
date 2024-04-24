ipower<-function(n)	{
	if (length(n) < 3 || max(n) < 3)
		return(list('richness' = NA, 'beta' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	library(stats4)
	S <- length(n)
	n2 <- n[n <= 2^14]
	S2 <- length(n2)
	s <- array(dim=2^14,data=0)
	t <- table(n2)
	s[as.numeric(names(t))] <- t
	u <- unique(n2)
	x <- 2:(2^14 + 2)
	like<-function(b)	{
		if (b <= 0)
			return(1e10)
		p <- -diff(1 / x^b) / exp(-b * log(2))
		if (is.nan(p[1]) || p[n2[S2]] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
			return(1e10)
		ll <- -sum(s[u] * log(p[u]))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
	b <- coef(stats4::mle(like,lower=list(b=0),upper=list(b=10),start=list(b=0.5)))
	if (b == 0 || b == 10)
		return(list('richness' = NA, 'beta' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
	aicc <- 2 * like(b) + 2 + 4 / (S2 - 2)
	x <- 2:(2^20 + 2)
	p <- -diff(1 / x^b) / exp(-b * log(2))
	return(list('richness' = as.numeric(S / exp(-b * log(2))), 'beta' = as.numeric(b), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
