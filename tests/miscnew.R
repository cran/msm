### TESTS OF MISCLASSIFICATION MODELS IMPLEMENTED WITH hmodel OBJECTS
### CHECK THESE against misc.R. 
source("local.R")
library(msm)
data(heart)

## Plain misc model using hmmCat

oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))

miscnew.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, death = 4, fixedpars=TRUE,
                hmodel=list(
                  hmmCat(prob=c(0.9, 0.1, 0, 0)),
                  hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                  hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent())
                )
miscnew.msm

if (developer.local) { 
    system.time(miscnew.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                                   qmatrix = oneway4.q, death = 4,  
                                   hmodel=list(
                                     hmmCat(prob=c(0.9, 0.1, 0, 0)),
                                     hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                                     hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                                   control = list(trace=1, REPORT=1), method="BFGS"
                                   ))
    print(miscnew.msm)
    if(interactive()) save(miscnew.msm, file="~/msm/devel/models/miscnew.msm.rda")
}

## Does misc model with no misc reduce to simple 
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
miscnew.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = twoway4.q,
                   hmodel=list(hmmCat(prob=c(1, 0, 0, 0)), hmmCat(prob=c(0, 1, 0, 0)), hmmCat(prob=c(0, 0, 1, 0)), hmmIdent()),
                   death = 4, fixedpars=TRUE)
miscnew.msm
miscnew.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = twoway4.q, death = 4, fixedpars=TRUE)
miscnew.msm

## Covs on misc probs, new way, with hmodel 

misccovnew.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, death = 4, fixedpars=TRUE,
                hmodel=list(
                  hmmCat(prob=c(0.9, 0.1, 0, 0)),
                  hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                  hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1),
                hcovinits = list(c(0.01,-0.013), c(0.02,-0.014,0.03,-0.015), c(0.04,-0.016), NULL),
                )
misccovnew.msm

if (developer.local) { 
    system.time(misccovnew.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                                      qmatrix = oneway4.q, death = 4, fixedpars=FALSE,
                                      hmodel=list(
                                        hmmCat(prob=c(0.9, 0.1, 0, 0)),
                                        hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                                        hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                                      hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1),
                                      hcovinits = list(c(0.01,-0.013), c(0.02,-0.014,0.03,-0.015), c(0.04,-0.016), NULL),
                                      control = list(trace=1, REPORT=1), method="BFGS"
                                      ))
    print(misccovnew.msm)
    if(interactive()) save(misccovnew.msm, file="~/msm/devel/models/misccovnew.msm.rda")

    ## Misclassification-specific output functions with misc model specified new way

    if(interactive()) load("~/msm/devel/models/miscnew.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccovnew.msm.rda")

    ## Don't allow ematrix 
    try(ematrix.msm(miscnew.msm))
                                        #ematrix.msm(miscnew.msm)[c("estimates","SE")]
                                        #print(ematrix.msm(miscnew.msm), digits=2)
    ## Instead extract parameters like so
    print(miscnew.msm$hmodel)
    print(viterbi.msm(miscnew.msm)[1:50,])
    print(viterbi.msm(miscnew.msm)[viterbi.msm(miscnew.msm)$subject==100063,])
    ## Don't alllow odds.msm. need model fitted with ematrix.
    try(odds.msm(misccovnew.msm))
    
    ## Non misclassification-specific output functions.  All OK. 

    print(qmatrix.msm(misccovnew.msm))
    print(sojourn.msm(misccovnew.msm))
    print(pmatrix.msm(misccovnew.msm, 10))
    print(qratio.msm(misccovnew.msm, c(1,2), c(2,3), cl=0.99))
    print(prevalence.msm(misccovnew.msm))
    print(summary.msm(misccovnew.msm))
    if (interactive()) plot.msm(misccovnew.msm)
                                        #coef.msm(misccovbothnew.msm)
                                        #hazard.msm(misccovbothnew.msm)
    print(transient.msm(misccovnew.msm))
    print(absorbing.msm(misccovnew.msm))
    print(totlos.msm(misccovnew.msm))
    print(logLik.msm(misccovnew.msm))
}

## Covariate initial values defaulting to 0

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, death = 4, fixedpars=TRUE,
                hmodel=list(
                  hmmCat(prob=c(0.9, 0.1, 0, 0)),
                  hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                  hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent()),
                hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1),
                )

misc.msm

## Covs on misc probs, new way, with ematrix 
## Don't allow, since hcovariates doesn't logically correspond to ematrix. 
try ( misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                         qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=1:17,
                         hcovariates=list(~dage + sex, ~dage + sex, ~dage + sex, ~1),
                         hcovinits = list(c(0.01,0.013), c(0.01,0.013,0.01,0.013), c(0.01,0.013), NULL),
                         control = list(trace=1, REPORT=1), method="BFGS") )

## initprobs

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, death = 4, fixedpars=TRUE,
                initprobs=c(0.7, 0.1, 0.1, 0.1),
                hmodel=list(
                  hmmCat(prob=c(0.9, 0.1, 0, 0)),
                  hmmCat(prob=c(0.1, 0.8, 0.1, 0)),
                  hmmCat(prob=c(0, 0.1, 0.9, 0)), hmmIdent())
                )
misc.msm
