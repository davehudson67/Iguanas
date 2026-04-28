a <- 0.05
b <- 0.02
c1 <- exp(-5)
lowerBound <- 0

## test non-truncated samples
sims <- replicate(10000, rGompertzMakehamLB(1, a, b, c1, 3))
hist(sims, freq = FALSE)

## add PDF
lines(0:100, dGompertzMakehamLB(0:100, a, b, c1, 3), col = "blue", lwd = 2)


## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- replicate(10000, rGompertzMakehamLB(1, a, b, c1, 3))
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pGompertzMakehamLB(1:80, a, b, c1, 3, lower.tail = 0), col = "blue")
      