### TESTS OF MISCLASSIFICATION MODELS
### SPECIFIED USING THE OLD ematrix SYNTAX
source("local.R")
library(msm)
data(heart)

oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))

## Plain misc model with no covs
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS")
misc.msm

if (developer.local) { 
    system.time(misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                                qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                                control = list(trace=1, REPORT=1), method="BFGS")) # new 43, old 129
    print(misc.msm)
    if(interactive()) save(misc.msm, file="~/msm/devel/models/misc.msm.rda")
}

## Does misc model with no misc reduce to simple 
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
nomisc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), death = 4, fixedpars=TRUE)
nomisc.msm
simple.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = twoway4.q, death = 4, fixedpars=TRUE)
simple.msm

## Covs on transition rates
misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars = TRUE, 
                control = list(trace=1, REPORT=1), method="BFGS",
                covariates = ~ sex, covinits=list(sex=rep(0.1, 5)))
misccov.msm

## Covs on misc probs, old way. 
misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                   qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE,
                   misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                   control = list(trace=1, REPORT=1), method="BFGS")
misccov.msm

if (developer.local) { 
    system.time(misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                                   qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                                   misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                                   control = list(trace=1, REPORT=1), method="BFGS"))
    print(misccov.msm)
    if(interactive()) save(misccov.msm, file="~/msm/devel/models/misccov.msm.rda")
}

## Covs on both
misccovboth.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE, 
                       control = list(trace=1, REPORT=1), method="BFGS",
                       covariates = ~ sex, covinits=list(sex=rep(0.1, 5)),
                       misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016))
                       )
misccovboth.msm

if (developer.local) { 
    system.time(misccovboth.msm <- msm(state ~ years, subject = PTNUM, data = heart, 
                                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=FALSE, 
                                       control = list(trace=1, REPORT=1), method="BFGS",
                                       covariates = ~ sex, covinits=list(sex=rep(0.1, 5)),
                                       misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016))
                                       ))
    if(interactive()) save(misccovboth.msm, file="~/msm/devel/models/misccovboth.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccovboth.msm.rda")
    print(misccovboth.msm)

##########    OUTPUT FUNCTIONS    ###################

    if(interactive()) load("~/msm/devel/models/misc.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccov.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccovboth.msm.rda")

    print(ematrix.msm(misc.msm)  )
    print(ematrix.msm(misc.msm)[c("estimates","SE")])
    print(ematrix.msm(misc.msm), digits=2)
    print(viterbi.msm(misc.msm)[1:50,])
    print(viterbi.msm(misc.msm)[viterbi.msm(misc.msm)$subject==100063,])
    print(odds.msm(misccov.msm))

    print(ematrix.msm(misccov.msm)  )
    print(ematrix.msm(misccov.msm)[c("estimates","SE")])
    print(ematrix.msm(misccov.msm, covariates=0)  )
    print(ematrix.msm(misccov.msm, covariates=0)[c("estimates","SE")])
    print(ematrix.msm(misccov.msm, covariates=list(dage=50, sex=0))  )
    print(ematrix.msm(misccov.msm, covariates=list(dage=50, sex=0))[c("estimates","SE")])

    print(qmatrix.msm(misccovboth.msm))  # Slightly different model fit in version 0.5
    print(qmatrix.msm(misccovboth.msm)[c("estimates","SE")])
    print(qmatrix.msm(misccovboth.msm, covariates=0)  )
    print(qmatrix.msm(misccovboth.msm, covariates=0)[c("estimates","SE")])
    print(qmatrix.msm(misccovboth.msm, covariates=list(sex=1))  )
    print(qmatrix.msm(misccovboth.msm, covariates=list(sex=1))[c("estimates","SE")])

### Non misclassification-specific output functions 

    print(qmatrix.msm(misccov.msm))
    print(sojourn.msm(misccov.msm))
    print(pmatrix.msm(misccov.msm, 10))
    print(qratio.msm(misccov.msm, c(1,2), c(2,3), cl=0.99))
    print(prevalence.msm(misccov.msm))
    print(summary.msm(misccov.msm))
    if (interactive()) plot.msm(misccov.msm)
    print(coef.msm(misccovboth.msm))
    print(hazard.msm(misccovboth.msm))
    print(transient.msm(misccov.msm))
    print(absorbing.msm(misccov.msm))
    print(totlos.msm(misccov.msm))
    print(logLik.msm(misccov.msm))

}


##########    OTHER FEATURES      ###################

if (developer.local) { 
    ## Baseline intens constraints
    misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=4:7,
                    qconstraint = c(1, 2, 1, 2, 3), 
                    control = list(trace=1, REPORT=1), method="BFGS")
    print(misc.msm)
    print(qmatrix.msm(misc.msm)[c("estimates","SE")])

    ## Baseline misc constraints
    ematrix2 <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.11, 0),c(0, 0.11, 0, 0),c(0, 0, 0, 0))
    misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                    qmatrix = oneway4.q, ematrix=ematrix2, death = 4, fixedpars=1:5,
                    econstraint = c(1, 1, 2, 2), 
                    control = list(trace=1, REPORT=1), method="BFGS")
    print(misc.msm)
    ematrix.msm(misc.msm)[c("estimates","SE")]

    ## intens covariate constraints
    ## Give replicated inits for replicated cov effs, as that's consistent with constraints on q, e and h.

    misc.msm <- msm(state ~ years, subject = PTNUM, data = heart, fixedpars=c(1:5, 9:12), 
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                    control = list(trace=1, REPORT=1), method="BFGS",
                    covariates = ~ sex, covinits=list(sex=c(0, 0, 0.1, 0, 0)), 
                    constraint = list(sex = c(1, 2, 1, 2, 3))  )
    print(misc.msm)
    print(qmatrix.msm(misc.msm, covariates=0))

    ## misc covariate constraints. Version 0.4.1 breaks when calculating SEs here (singular in dgesv).

    misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=c(1:5), 
                       misccovariates = ~dage + sex,
                       misccovinits = list(dage=c(0.01,0.01,0.001,0.001), sex=c(0.0131,0.0132,0.0133,0.0134)),
                       miscconstraint = list(dage = c(1, 1, 2, 2)), 
                       control = list(trace=1, REPORT=1), method="BFGS")
    print(misccov.msm)
    print(ematrix.msm(misccov.msm)[c("estimates","SE")])
    print(ematrix.msm(misccov.msm, covariates=0)[c("estimates","SE")])

    ## fixedpars for misc covariates.  Parameters are ordered within covariate, within parameter, within state.

    misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                       misccovariates = ~dage + sex,
                       misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                       fixedpars = c(10, 11, 12, 15),
                       control = list(trace=1, REPORT=1), method="BFGS")
    print(misccov.msm)
}

## multiple death states (Jean-Luc's data)
## Misclassification between states 2 and 3 

if (developer.local) { 
  c2.df <- read.table("~/msm/tests/jeanluc/donneesaveccancerPT.txt", header=TRUE)
  print(statetable.msm(state, PTNUM, c2.df))
  qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
  ex <- rbind( c(0, 0, 0, 0, 0), c(0, 0, 0.1, 0, 0), c(0, 0.1, 0, 0, 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0) ) 
  c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx, ematrix=ex, death=c(4, 5), method="BFGS", fixedpars = TRUE, 
                control=list(trace=2, REPORT=1, fnscale=100000))
  print(c2.msm)
  print(logLik(c2.msm))
  
  ## multiple death states specified using an obstype vector 
  d45 <- rep(1, nrow(c2.df)); d45[c2.df$state %in% c(4,5)] <- 3
  c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx, ematrix=ex, obstype=d45, method="BFGS", fixedpars = TRUE, 
                control=list(trace=2, REPORT=1, fnscale=100000))
  print(c2.msm)
  print(logLik(c2.msm))
}

## exact times
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, exacttimes=TRUE, fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS") # should warn about redundant death argument
misc.msm

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, exacttimes=TRUE, fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS")
misc.msm

misc.msm$minus2loglik

## exact times specified using an obstype vector
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, obstype=rep(2, nrow(heart)), fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS")
misc.msm

## initprobs
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE, initprobs=c(0.7, 0.1, 0.1, 0.1),
                control = list(trace=1, REPORT=1), method="BFGS")
misc.msm

## Censored states
## Different in 0.4 or less, censored states are not subject to misclassification in >= 0.4.1

heart.cens <- heart
heart.cens$state[heart$state==4][1:50] <- 99
heart.cens2 <- heart
heart.cens2$state[heart$state==4][1:50] <- 99
heart.cens2$state[heart$state==4][51:100] <- 999
heart.cens3 <- heart
ns <- c(heart$state[2:nrow(heart)], 0)
heart.cens3$state[heart$state==4][1:50] <- 99
heart.cens3$state[ns==4][1:50] <- 999

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens,
                qmatrix = oneway4.q, ematrix=ematrix, death=TRUE, censor=99, fixedpars=TRUE)
misc.msm

## Two types of censoring
misc.msm <- msm(state ~ years, subject=PTNUM, data=heart.cens2, qmatrix=oneway4.q, ematrix=ematrix, censor=c(99, 999), death=4, censor.states=list(c(1,2,3), c(2,3)), fixedpars=TRUE)
misc.msm

## Does misc model with no misc reduce to simple, with censoring 

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens,
                qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), censor=99, death=TRUE, fixedpars=TRUE)
misc.msm
simple.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens, qmatrix = twoway4.q, death=TRUE, censor=99, fixedpars=TRUE)
simple.msm

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens,
                qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), censor=99, fixedpars=TRUE)
misc.msm
simple.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens, qmatrix = twoway4.q, censor=99, fixedpars=TRUE)
simple.msm

## Viterbi with non-HMM model
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=TRUE)
viterbi.msm(heart.msm)[1:50,] # no error, returns observed states. 
