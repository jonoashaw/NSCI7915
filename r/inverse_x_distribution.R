inv<-function(n)	{
  if (length(n) < 3 || max(n) < 3)
    return(list('richness' = NA, 'lambda' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
  library(stats4)
  S <- length(n)
  n2 <- n[n <= 2^14]
  S2 <- length(n2)
  s <- array(dim=2^14,data=0)
  t <- table(n2)
  s[as.numeric(names(t))] <- t
  u <- unique(n2)
  x <- (1:(2^14 + 1))^0.5
  like<-function(l)	{
    if (l <= -1)
      return(1e10)
    p <- -diff(0.5^(1 / x^l))
    p <- p / sum(p)
    if (is.nan(p[1]) || p[n2[S2]] < 1e-100 || min(p[u]) == 0 || sum(p[u]) > 1 - 1e-120)
      return(1e10)
    ll <- -sum(s[u] * log(p[u]))
    if (is.infinite(ll) || is.nan(ll))
      return(1e10)
    ll
  }
  l <- coef(stats4::mle(like,lower=list(l=-1),upper=list(l=1e6),start=list(l=0.01)))
  if (l == -1 || l == 1e6)
    return(list('richness' = NA, 'lambda' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))
  aicc <- 2 * like(l) + 2 + 4 / (S2 - 2)
  x <- (1:(2^20 + 1))^0.5
  p <- -diff(0.5^(1 / x^l))
  p <- p / sum(p)
  return(list('richness' = as.numeric(S / exp(-l)), 'lambda' = as.numeric(l), 'AICc' = aicc, 'fitted.RAD' = sadrad(length(n),p), 'fitted.SAD' = p[1:2^12]))
}
