source("local.R")
library(msm)

## New Pade matrix exponential

A <- matrix(c(-0.11, 0.01, 0.001,  0.2, -0.2, 0,  0, 0, 0), nrow=3, byrow=TRUE)
MatrixExp(A, method="pade")
MatrixExp(A, method="series")


### New d-p-q-r functions

## Truncated normal distribution
set.seed(220676)
(rl <- rnorm(10))
dtnorm(rl)
dnorm(rl)
dtnorm(rl, mean=2, sd=1.2)
dnorm(rl, mean=2, sd=1.2)
dtnorm(rl, mean=2, sd=1.2, lower=seq(-4,5))
ptnorm(c(-1000, 0, 1000))
ptnorm(c(-1, 0, 1), mean=c(0,1,2), sd=c(1,2,3), lower=c(-2,1,0))
qtnorm(ptnorm(rl))
try(qtnorm(c(-1, 0, 1, 2)))
rtnorm(10, lower=0, upper=2)

dnorm(2)
dmenorm(2)
dmenorm(c(-2, 0, 2), mean=c(0,1,2))
dmenorm(c(-2, 0, 2), mean=c(0,1,2), lower=c(-3,-2,-1), sderr=c(2,3,4))
pmenorm(c(-2, 0, 2))
pnorm(c(-2, 0, 2))
pmenorm(c(-2, 0, 2), sderr=5)
pmenorm(c(-2, 0, 2), sderr=5, meanerr=10)
qmenorm(pmenorm(c(-2, 0, 2), sderr=5, lower=0), sderr=5, lower=0)
qmenorm(pmenorm(c(-2, 0, 2)))
rmenorm(10, lower=0, upper=2, sderr=0.2)

dmeunif(c(-2, 0, 0.7, 1, 2))
dunif(c(-2, 0, 0.7, 1, 2))
dmeunif(c(-2, 0, 0.7, 1, 2))
dunif(c(-2, 0, 0.7, 1, 2), min=-3:1, max=4:8)
dmeunif(c(-2, 0, 0.7, 1, 2), lower=-3:1, upper=4:8)
dmeunif(c(-2, 0, 0.7, 1, 2), lower=-3:1, upper=4:8, sderr=1)
pmeunif(c(0.1, 0.5, 0.9))
punif(c(0.1, 0.5, 0.9))
pmeunif(c(0.1, 0.5, 0.9), sderr=5)
pmeunif(c(0.1, 0.5, 0.9), sderr=5, meanerr=10)
qmeunif(pmeunif(c(0.1, 0.5, 0.9), sderr=5, lower=-1), sderr=5, lower=-1)
qmeunif(pmeunif(c(0.1, 0.5, 0.9)))

## Example in help(deltamethod)
## Simple linear regression, E(y) = alpha + beta x 
x <- 1:100
set.seed(220676)
y <- rnorm(100, 4*x, 5)
toy.lm <- lm(y ~ x)
(estmean <- coef(toy.lm))
estvar <- summary(toy.lm)$cov.unscaled

## Estimate of (1 / (alphahat + betahat))
1 / (estmean[1] + estmean[2])
## Approximate standard error
deltamethod (~ 1 / (x1 + x2), estmean, estvar) 
