a1 <- exp(-2)
a2 <- exp(-11)
b1 <- 1.2
b2 <- 0.11
c <- exp(-5) 
lowerBound <- 0

## test non-truncated samples
sims <- replicate(10000, rSilerLB(1, a1, a2, b1, b2, c, 3))
hist(sims, freq = FALSE)

## add PDF
lines(0:100, dSilerLB(0:100, a1, a2, b1, b2, c, 3), col = "blue", lwd = 2)


## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- replicate(10000, rSilerLB(1, a1, a2, b1, b2, c, 3))
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pSilerLB(1:80, a1, a2, b1, b2, c, 3, lower.tail = 0), col = "blue")
      