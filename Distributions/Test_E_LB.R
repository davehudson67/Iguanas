r <- 2
lowerBound <- 12

## test non-truncated samples
sims <- replicate(10000, rExponentialLB(1, r, lowerBound))
hist(sims, freq = FALSE)

## add PDF
lines(0:100, dExponentialLB(0:100, r, lowerBound), col = "blue", lwd = 2)


## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- replicate(10000, rExponentialLB(1, r, lowerBound))
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pExponentialLB(1:80, r, lowerBound, lower.tail = 0), col = "blue")
