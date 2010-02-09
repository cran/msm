source("local.R")
library(msm)
#library(msm, lib.loc="~/lib/R")
data(psor)

psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))

## all panel data 

psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, 
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), control = list(REPORT=1,trace=2),
                method="BFGS")
pears <- pearson.msm(psor.msm)

stopifnot(isTRUE(all.equal(4027, pears$test$stat, tol=1)))

## exact death times 

psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, death=c(4), 
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), control = list(REPORT=1,trace=2),
                method="BFGS")
pears <- pearson.msm(psor.msm)

stopifnot(isTRUE(all.equal(1799, pears$test$stat, tol=1)))

## multiple death times 

psor2.q <- rbind(c(0,0.1,0,0,0),c(0,0,0.1,0,0),c(0,0,0,0.1,0.1),c(0,0,0,0,0),c(0,0,0,0,0))
psor2 <- psor; psor2$state[psor2$state==4][1:10] <- 5 
psor.msm <- msm(state ~ months, subject=ptnum, data=psor2, qmatrix = psor2.q, covariates = ~ollwsdrt+hieffusn, death=c(4,5), 
                constraint = list(hieffusn=c(1,1,1,1),ollwsdrt=c(1,1,2,1)), control = list(REPORT=1,trace=2),
                method="BFGS")
pears <- pearson.msm(psor.msm)

stopifnot(isTRUE(all.equal(2047, pears$test$stat, tol=1)))

cat("pearson.R: ALL TESTS PASSED\n")
