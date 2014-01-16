## Test integrating over covariates in expected prevalences.  New feature in 1.2.

if (0) {
    
library(msm)
data(psor)

for (i in list.files("../R", ".+\\.R$")) source(paste("../R/",i,sep=""))

psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,  covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)))
x <- psor.msm
times=NULL; interp="start"; censtime=Inf

## No covariates
psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q)
expected.msm(psor.msm)

## One covariate
psor.msm <- msm(state ~ months, covariates=~ollwsdrt, subject=ptnum, data=psor, qmatrix = psor.q)
get.covhist(psor.msm)
observed.msm(psor.msm)
expected.msm(psor.msm)
expected.msm(psor.msm, covariates="mean")

## Two covariates
psor.msm <- msm(state ~ months, covariates=~ollwsdrt+hieffusn, subject=ptnum, data=psor, qmatrix = psor.q)
get.covhist(psor.msm)
observed.msm(psor.msm)
expected.msm(psor.msm)
expected.msm(psor.msm, covariates="mean")
## just a bit different from the mean, as expected

## PCI
psor.msm <- msm(state ~ months, covariates=~ollwsdrt+hieffusn, subject=ptnum, data=psor, qmatrix = psor.q, pci=c(5,10))
x <- psor.msm
(covhist <- get.covhist(psor.msm))
observed.msm(psor.msm)
expected.msm(psor.msm)
expected.msm(psor.msm, covariates="mean")
## just a bit different from the mean, as expected

## Manual time-inhomog, also continuous covs
psor2 <- psor; psor2$ollwsdrt <- rnorm(nrow(psor2), 0, 1)
psor.msm <- msm(state ~ months, covariates=~ollwsdrt+hieffusn, subject=ptnum, data=psor2, qmatrix = psor.q, pci=c(5,10))
psor.msm
(covhist <- get.covhist(psor.msm))
observed.msm(psor.msm)
expected.msm(psor.msm)
expected.msm(psor.msm, covariates="mean")
## just a bit different from the mean, as expected
plot.prevalence.msm(psor.msm)  ## Takes a minute or so with continuous covs. 

## Misclassification models
misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, covariates = ~ cumrej, 
                qmatrix = oneway4.q, ematrix=ematrix, death = 4)
misc.msm
observed.msm(misc.msm)
expected.msm(misc.msm)
expected.msm(misc.msm, covariates="mean")

## CIs
psor.msm <- msm(state ~ months, covariates=~ollwsdrt, subject=ptnum, data=psor, qmatrix = psor.q)
prevalence.msm(psor.msm, ci="normal", B=10)
prevalence.msm(psor.msm, covariates="mean")

## Plotting - is this slow? 
plot.prevalence.msm(psor.msm)
plot.prevalence.msm(psor.msm, covariates="mean")

## summary.msm
summary(psor.msm)

## Test subset arg to observed.
## Simulate a model with a significant covariate effect
sim.df <- na.omit(psor); colnames(sim.df)[colnames(sim.df)%in%c("ptnum","months")] <- c("subject","time")
sim.df <- simmulti.msm(sim.df, qmatrix=psor.q, covariates=list(ollwsdrt=c(3,3,3)))
psor.msm <- msm(state ~ time, covariates=~ollwsdrt, subject=subject, data=sim.df, qmatrix = psor.q)
plot.prevalence.msm(psor.msm, interp="midpoint")
plot.prevalence.msm(psor.msm, covariates="mean", interp="midpoint")

subs <- psor.msm$data$subject[!duplicated(psor.msm$data$subject) & psor.msm$data$cov.orig$ollwsdrt==0]
observed.msm(psor.msm, subset=subs)$obsperc
plot.prevalence.msm(psor.msm, covariates=list(ollwsdrt=0), subset=subs, interp="midpoint")
plot.prevalence.msm(psor.msm, covariates=list(ollwsdrt=1), subset=subs, interp="midpoint")

subs <- psor.msm$data$subject[!duplicated(psor.msm$data$subject) & psor.msm$data$cov.orig$ollwsdrt==1]
observed.msm(psor.msm, subset=subs)$obsperc
plot.prevalence.msm(psor.msm, covariates=list(ollwsdrt=0), subset=subs, interp="midpoint")
plot.prevalence.msm(psor.msm, covariates=list(ollwsdrt=1), subset=subs, interp="midpoint")

} 
