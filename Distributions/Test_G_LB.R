a <- 0.05
b <- 0.02
lowerBound <- 12

## test non-truncated samples
sims <- replicate(10000, rGompertzLB(1, a, b, lowerBound))
hist(sims, freq = FALSE)

## add PDF
lines(0:100, dGompertzLB(0:100, a, b, lowerBound), col = "blue", lwd = 2)


## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- replicate(10000, rGompertzLB(1, a, b, lowerBound))
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pGompertzLB(1:80, a, b, lowerBound, lower.tail = 0), col = "blue")
