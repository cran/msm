source("local.R")
library(msm)
# for (i in list.files("../R", ".+\\.R$")) source(paste("../R/",i,sep=""))

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))

## all Inf
Q <- twoway4.q
efpt.msm(qmatrix=Q, tostate=3)
## none Inf
Q <- twoway4.q[1:3,1:3]; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
efpt.msm(qmatrix=Q, tostate=3)
## calculate by hand
Q <- rbind(c(-0.25,0.25,0), c(0.166, -0.332, 0.166), c(0, 0.25, -0.25))
stopifnot(all.equal(efpt.msm(qmatrix=Q, tostate=3)[c(1,2)], solve(-Q[1:2,1:2], c(1,1)),tol=1e-06))
## some Inf and some not.
Q <- twoway4.q; Q[2,4] <- Q[2,1] <- 0; diag(Q) <- 0; diag(Q) <- -rowSums(Q)
efpt.msm(qmatrix=Q, tostate=3)
## test against sim
if (0) { 
n <- 10000
t <- numeric(n)
for (i in 1:n){
sim <- sim.msm(qmatrix=Q, maxtime=200)
t[i] <- min(sim$times[sim$states==3])
}
mean(t)
}

## from a fitted model.  
data(psor)
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
efpt.msm(psor.msm, tostate=c(2))
efpt.msm(psor.msm, tostate=c(3))
efpt.msm(psor.msm, tostate=c(2,3))

if (0) { 
efpt.msm(psor.msm, tostate=c(3), ci="normal", B=1000)
efpt.msm(psor.msm, tostate=c(3), ci="boot", B=100)

## Test boot ci for totlos accepts integrate options
set.seed(12082012)
totlos.msm(psor.msm, ci="normal", B=10, subdivisions=50)
totlos.msm(psor.msm, subdivisions=5, rel.tol=10, abs.tol=10)
set.seed(12082012)
totlos.msm(psor.msm, ci="normal", B=10, subdivisions=5, rel.tol=10, abs.tol=10)
set.seed(12082012)
totlos.msm(psor.msm, ci="boot", B=10)
set.seed(12082012)
totlos.msm(psor.msm, ci="boot", B=10, subdivisions=5, rel.tol=10, abs.tol=10)
}
