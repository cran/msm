source("local.R")
library(msm)

## New Pade matrix exponential

A <- matrix(c(-0.11, 0.01, 0.001,  0.2, -0.2, 0,  0, 0, 0), nrow=3, byrow=TRUE)
me <- MatrixExp(A, method="pade")
stopifnot(isTRUE(all.equal(c(0.896703832431769, 0.171397960992687, 0, 0.00856989804963433, 0.81957474998506, 0, 0.00094726269518597, 9.0272890222537e-05, 1), as.numeric(me), tol=1e-06)))
me <- MatrixExp(A, method="series")
stopifnot(isTRUE(all.equal(c(0.896703832431769, 0.171397960992687, 0, 0.00856989804963433, 0.81957474998506, 0, 0.00094726269518597, 9.0272890222537e-05, 1), as.numeric(me), tol=1e-06)))
ev <- eigen(A)
me2 <- ev$vectors %*% diag(exp(ev$values)) %*% solve(ev$vectors)
stopifnot(isTRUE(all.equal(me2, me, tol=1e-06)))

### New d-p-q-r functions

## Truncated normal distribution
set.seed(220676)
(rl <- rnorm(10))

stopifnot(isTRUE(all.equal(dtnorm(rl), dnorm(rl), tol=1e-06)))
stopifnot(isTRUE(all.equal(dtnorm(rl, mean=2, sd=1.2), dnorm(rl, mean=2, sd=1.2), tol=1e-06)))

d <- dtnorm(rl, mean=2, sd=1.2, lower=seq(-4,5))
stopifnot(isTRUE(all.equal(c(0.260110259383406, 0.108097895222820, 0.0558659556833655, 0.160829438765247, 0.343919966894772, 0, 0, 0, 0, 0), d, tol=1e-06)))

stopifnot(isTRUE(all.equal(c(0, 0.5, 1), ptnorm(c(-1000, 0, 1000)), tol=1e-06)))

stopifnot(isTRUE(all.equal(c(0.139068959153926, 0, 0.156451685781240), ptnorm(c(-1, 0, 1), mean=c(0,1,2), sd=c(1,2,3), lower=c(-2,1,0)), tol=1e-06)))

stopifnot(isTRUE(all.equal(rl, qtnorm(ptnorm(rl)), tol=1e-03)))

try(qtnorm(c(-1, 0, 1, 2)))
set.seed(220676)
rtnorm(10, lower=0, upper=2)
set.seed(220676)
rt <- rtnorm(10)
set.seed(220676)
r <- rnorm(10)
stopifnot(isTRUE(all.equal(rt, r, tol=1e-06)))

## Measurement error distributions 

stopifnot(isTRUE(all.equal(dnorm(2), dmenorm(2), tol=1e-06)))
stopifnot(isTRUE(all.equal(dnorm(2, log=TRUE), dmenorm(2, log=TRUE), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.0539909665131881, 0.241970724519143, 0.398942280401433), dmenorm(c(-2, 0, 2), mean=c(0,1,2)), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.119536494085260, 0.120031723608082, 0.0967922982964366), dmenorm(c(-2, 0, 2), mean=c(0,1,2), lower=c(-3,-2,-1), sderr=c(2,3,4)), tol=1e-06)))
stopifnot(isTRUE(all.equal(pmenorm(c(-2, 0, 2)), pnorm(c(-2, 0, 2)), tol=1e-06)))
stopifnot(isTRUE(all.equal(pmenorm(c(-2, 0, 2), log.p=TRUE), pnorm(c(-2, 0, 2), log.p=TRUE), tol=1e-06)))
stopifnot(isTRUE(all.equal(pmenorm(c(-2, 0, 2), lower.tail=FALSE), pnorm(c(-2, 0, 2), lower.tail=FALSE), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.347443301205908, 0.500000000140865, 0.652556698813763), pmenorm(c(-2, 0, 2), sderr=5), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.00930146266876999, 0.0249300921973760, 0.0583322325986182), pmenorm(c(-2, 0, 2), sderr=5, meanerr=10), tol=1e-06)))

stopifnot(isTRUE(all.equal(qmenorm(pmenorm(c(-2, 0, 2), sderr=5, lower=0), sderr=5, lower=0), qmenorm(pmenorm(c(-2, 0, 2))), tol=1e-03)))

set.seed(220676)
rmenorm(10, lower=0, upper=2, sderr=0.2)

stopifnot(isTRUE(all.equal(c(0,1,1,0,0), dmeunif(c(-2, 0, 0.7, 1, 2)))))
stopifnot(isTRUE(all.equal(dunif(c(-2, 0, 0.7, 1, 2), min=-3:1, max=4:8), dmeunif(c(-2, 0, 0.7, 1, 2), lower=-3:1, upper=4:8), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.120192106440279, 0.139607083057178, 0.136490639905731, 0.120192106440279, 0.120192106440279), dmeunif(c(-2, 0, 0.7, 1, 2), lower=-3:1, upper=4:8, sderr=1), tol=1e-06)))
stopifnot(isTRUE(all.equal(pmeunif(c(0.1, 0.5, 0.9)), punif(c(0.1, 0.5, 0.9)), tol=1e-04)))
stopifnot(isTRUE(all.equal(c(0.468171571157871, 0.500000000120507, 0.531828429094026), pmeunif(c(0.1, 0.5, 0.9), sderr=5), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.0189218497312070, 0.0229301821964305, 0.0276311076816442), pmeunif(c(0.1, 0.5, 0.9), sderr=5, meanerr=10), tol=1e-06)))
stopifnot(isTRUE(all.equal(c(0.1, 0.5, 0.9), qmeunif(pmeunif(c(0.1, 0.5, 0.9), sderr=5, lower=-1), sderr=5, lower=-1), tol=1e-03)))
stopifnot(isTRUE(all.equal(c(0.1, 0.5, 0.9), qmeunif(pmeunif(c(0.1, 0.5, 0.9))), tol=1e-03)))


## Exponential distribution with piecewise constant hazard

stopifnot(isTRUE(all.equal(1, integrate(dpexp, 0, Inf)$value)))
rate <- c(0.1, 0.2, 0.05, 0.3)
t <- c(0, 10, 20, 30)
stopifnot(isTRUE(all.equal(1, integrate(dpexp, 0, Inf, rate=rate, t=t)$value, tol=1e-04)))
x <- rexp(10)
stopifnot(isTRUE(all.equal(dpexp(x), dexp(x))))
stopifnot(isTRUE(all.equal(dpexp(x, log=TRUE), log(dpexp(x)))))
stopifnot(isTRUE(all.equal(dpexp(x, log=TRUE), dexp(x, log=TRUE))))


stopifnot(ppexp(-5) == 0)
stopifnot(ppexp(0) == 0)
stopifnot(ppexp(Inf) == 1)
set.seed(22061976)
q <- rexp(10)
stopifnot(isTRUE(all.equal(pexp(q), ppexp(q))))
stopifnot(isTRUE(all.equal(pexp(q, log.p=TRUE), ppexp(q, log.p=TRUE))))
rate <- c(0.1, 0.2, 0.05, 0.3)
t <- c(0, 10, 20, 30)
stopifnot(ppexp(-5, rate, t) == 0)
stopifnot(ppexp(0, rate, t) == 0)
stopifnot(isTRUE(all.equal(1, ppexp(Inf, rate, t))))
stopifnot(isTRUE(all.equal(1, ppexp(9999999, rate, t))))
stopifnot(isTRUE(all.equal(pexp(c(5, 6, 7), rate[1]), ppexp(c(5, 6, 7), rate, t))))
try(ppexp(q, rate=c(1,2,3), t=c(1,2))) # rate and t different lengths
try(ppexp(q, rate=-4)) # negative rates, NaN
try(ppexp(q, rate=c(1,2,3), t=c(-1, 4, 6))) # first elt of t not zero 


set.seed(22061976)
p <- runif(10)
stopifnot(isTRUE(all.equal(qpexp(p), qexp(p), tol=1e-03)))
stopifnot(isTRUE(all.equal(qpexp(p, lower.tail=FALSE), qexp(p, lower.tail=FALSE), tol=1e-03)))
stopifnot(isTRUE(all.equal(qpexp(log(p), log.p=TRUE), qexp(log(p), log.p=TRUE), tol=1e-03)))
stopifnot(isTRUE(all.equal(p, ppexp(qpexp(p)), tol=1e-03)))
set.seed(22061976)
q <- rexp(10)
stopifnot(isTRUE(all.equal(q, qpexp(ppexp(q)), tol=1e-03)))


set.seed(220676)
rt <- rpexp(10)
set.seed(220676)
r <- rexp(10)
stopifnot(isTRUE(all.equal(rt, r, tol=1e-06)))


## Example in help(deltamethod)
## Simple linear regression, E(y) = alpha + beta x 
x <- 1:100
set.seed(220676)
y <- rnorm(100, 4*x, 5)
toy.lm <- lm(y ~ x)
(estmean <- coef(toy.lm))
estvar <- summary(toy.lm)$cov.unscaled

## Estimate of (1 / (alphahat + betahat))
stopifnot(isTRUE(all.equal(0.206982798128202, as.numeric(1 / (estmean[1] + estmean[2])))))
## Approximate standard error
stopifnot(isTRUE(all.equal(0.00850451118374906, deltamethod (~ 1 / (x1 + x2), estmean, estvar))))

cat("utils.R: ALL TESTS PASSED\n")
