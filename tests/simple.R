source("local.R")
library(msm)
data(heart)

### TESTS FOR SIMPLE NON-HIDDEN MARKOV MODELS

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
twoway3.q <- rbind(c(-0.5, 0.25, 0), c(0.166, -0.498, 0.166), c(0, 0.25, -0.5))
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well","Mild","Severe","Death")
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
twoway4.i <- twoway4.q; twoway4.i[twoway4.i!=0] <- 1
oneway4.i <- oneway4.q; oneway4.i[oneway4.i!=0] <- 1

### HEART DATA
statetable.msm(state, PTNUM, data=heart)
(cinits <- crudeinits.msm(state ~ years, PTNUM, data=heart, qmatrix=twoway4.q))

## Simple model
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                 qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1)) 
heart.msm
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                 qmatrix = cinits, death = TRUE, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1)) 
heart.msm
if (developer.local) {
    system.time(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                                 qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,
                                 method="BFGS", control=list(trace=5, REPORT=1)) ) # 9.81 on new, 55.14 on old! 
    print(heart.msm)
}

## auto-generated initial values. 
state.g <- heart$state; time.g <- heart$years; subj.g <- heart$PTNUM
heart.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = twoway4.i, gen.inits=TRUE, fixedpars=TRUE)
heart.msm

if (developer.local) {
    system.time(heart.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = twoway4.i, gen.inits=TRUE, fixedpars=TRUE))
    print(heart.msm)
    heart.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = crudeinits.msm(state ~ years, PTNUM, twoway4.i, data=heart), fixedpars=TRUE)
    print(heart.msm)
}

## Covariates
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                 qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE, 
                 covariates = ~ sex, covinits = list(sex=rep(0.01, 7)), # , dage=rep(0, 7)),
                 method="BFGS", control=list(trace=5, REPORT=1))
heart.msm$minus2loglik 
heart.msm

if (developer.local) {
    system.time(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                                 qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, 
                                 covariates = ~ sex, method="BFGS", control=list(trace=5, REPORT=1))) # 44.13 on new, 260.14 on old
    print(heart.msm) # Close enough. 
    print(qmatrix.msm(heart.msm)[c("estimates","SE")])
    print(sojourn.msm(heart.msm))
}

if (developer.local) {
    ## Baseline constraints
    system.time(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                                 qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,
                                 qconstraint = c(1,1,2,2,2,3,3),
                                 method="BFGS", control=list(trace=2, REPORT=1)
                                 )) # 3.22 on new, 17.55 on old
    print(heart.msm) # 4116.226864
    print(qmatrix.msm(heart.msm)[c("estimates","SE")])
    print(sojourn.msm(heart.msm))

    ## Covariate constraints. 
    system.time(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, 
                                 covariates = ~ sex, covinits = list(sex=rep(0.01, 7)), constraint=list(sex=c(1,2,3,1,2,3,2)), 
                                 method="BFGS", control=list(trace=1, REPORT=1))) # 21.2 on new, 86.58 on old
    print(heart.msm)
    print(qmatrix.msm(heart.msm)[c("estimates","SE")])
}

## Constraints with psoriatic arthritis data
data(psor)
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)),
                fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
psor.msm
qmatrix.msm(psor.msm)[c("estimates","SE")]

## No death state
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, 
                 qmatrix = twoway4.q, death = FALSE, fixedpars=TRUE)
heart.msm

## Multiple death states (Jean-Luc's data)
if (developer.local) { 
    c2.df <- read.table("~/msm/tests/jeanluc/donneesaveccancerPT.txt", header=TRUE)
    qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                  qmatrix=qx, death=c(4, 5), method="BFGS", fixedpars = 1:5, 
                  control=list(trace=2, REPORT=1, fnscale=100000))
    print(c2.msm$minus2loglik)
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                  qmatrix=qx, method="BFGS", fixedpars = 1:5, 
                  control=list(trace=2, REPORT=1, fnscale=100000))
    print(c2.msm$minus2loglik)

    ## Same using an "obstype" vector. 
    obstype <- ifelse(c2.df$state %in% c(4,5), 3, 1)
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx,
                  obstype=obstype, method="BFGS", fixedpars = 1:5, 
                  control=list(trace=2, REPORT=1, fnscale=100000))
    print(c2.msm$minus2loglik)
    obstype <- rep(1, length(c2.df$state))
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, obstype=obstype, 
                  qmatrix=qx, method="BFGS", fixedpars = 1:5, 
                  control=list(trace=2, REPORT=1, fnscale=100000))
    print(c2.msm$minus2loglik)

### G Marshall's diabetic retinopathy data
    marsh.df <- read.table("~/msm/tests/markov/test.dat", col.names=c("subject","eyes","time","duration","hba1"))
    marsh.df$hba1 <- marsh.df$hba1 - mean(marsh.df$hba1)
    marsh.msm <-
      msm(eyes ~ time, subject=subject, qmatrix = rbind(c(0,0.02039,0,0), c(0.007874,0,0.01012,0), c(0,0.01393,0,0.01045), c(0,0,0,0)),
          covariates = ~ hba1, data = marsh.df, fixedpars=TRUE)
    marsh.msm <-
      msm(eyes ~ time, subject=subject, qmatrix = rbind(c(0,0.02039,0,0), c(0.007874,0,0.01012,0), c(0,0.01393,0,0.01045), c(0,0,0,0)),
          covariates = ~ hba1, data = marsh.df, control=list(trace=1, REPORT=1))
    print(qmatrix.msm(marsh.msm, covariates=0))
    print(marsh.msm)

}

## Exact times (BOS data)
fiveq <- rbind(c(0,0.01,0,0,0.002), c(0,0,0.07,0,0.01), c(0,0,0,0.07,0.02), c(0,0,0,0,0.03), c(0,0,0,0,0))
data(bos)
(msmtest5 <- msm(state ~ time, qmatrix = fiveq,  subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=1:7))
(msmtest5 <- msm(state ~ time, qmatrix = fiveq,  subject = ptnum, data = bos, obstype=rep(2, nrow(bos)), fixedpars=1:7))
(msmtest5 <- msm(state ~ time, qmatrix = fiveq,  subject = ptnum, data = bos, obstype=rep(1, nrow(bos)), fixedpars=1:7))
## Death and exact times (should be same!) 
(msmtest5 <- msm(state ~ time, qmatrix = fiveq,  subject = ptnum, data = bos, death=5, obstype=rep(2, nrow(bos)), exacttimes=TRUE, fixedpars=1:7))
(msmtest5 <- msm(state ~ time, qmatrix = fiveq,  subject = ptnum, data = bos, death=5, obstype=rep(2, nrow(bos)), fixedpars=1:7))
## Autogenerated inits (start at MLE, shouldn't need any optimisation)
fiveq.i <- fiveq; fiveq.i[fiveq.i!=0] <- 1
(msmtest5 <- msm(state ~ time, qmatrix = fiveq.i, gen.inits=TRUE, subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=1:7))
(msmtest5 <- msm(state ~ time, qmatrix = fiveq.i, gen.inits=TRUE, subject = ptnum, data = bos, exacttimes=TRUE))

### Aneurysm dataset different in 0.5 (not fromto, includes imputed initial states)


##########    OUTPUT FUNCTIONS    ###################

qmatrix.msm(psor.msm)
qmatrix.msm(psor.msm, covariates=list(hieffusn=0.1, ollwsdrt=0.4))
try(qmatrix.msm(psor.msm, covariates=list(hieffusn=0.1, foo=0.4))) # deliberate error
qmatrix.msm(psor.msm, covariates=list(hieffusn=0.1, ollwsdrt=0.4), cl=0.99)
qmatrix.msm(psor.msm, covariates=list(hieffusn=0.1, ollwsdrt=0.4), sojourn=TRUE)$sojourn

sojourn.msm(psor.msm, covariates=list(hieffusn=0.1, ollwsdrt=0.4))
sojourn.msm(psor.msm, covariates=list(hieffusn=0.1, ollwsdrt=0.4), cl=0.99)

pmatrix.msm(psor.msm, t=10)
try(pmatrix.msm(psor.msm, t=10, covariates=list(hieffusn=0.1))) # deliberate error
pmatrix.msm(psor.msm, t=10, covariates=list(hieffusn=0.1, ollwsdrt=0.2))

qratio.msm(psor.msm, c(1,2), c(2,3))
qratio.msm(psor.msm, c(1,2), c(2,3), cl=0.99)
qratio.msm(psor.msm, c(1,1), c(2,3))
qratio.msm(psor.msm, c(2,2), c(2,3))

prevalence.msm(psor.msm)
prevalence.msm(psor.msm, times=seq(0,60,5))

summary.msm(psor.msm)

print(interactive())
if (interactive()) plot.msm(psor.msm)
if (interactive()) plot.msm(psor.msm, from=c(1,3), to=4, range=c(10,30))
if (interactive()) plot.msm(psor.msm, from=c(1,2), to=4, range=c(10,80), legend.pos=c(70,0.1))

coef.msm(psor.msm)

hazard.msm(psor.msm)
hazard.msm(psor.msm, hazard.scale=2)
hazard.msm(psor.msm, hazard.scale=c(1,2))

transient.msm(psor.msm)

absorbing.msm(psor.msm)

totlos.msm(psor.msm)
totlos.msm(psor.msm, fromt=1, tot=30)
totlos.msm(psor.msm, start=2, fromt=10, tot=30)

logLik.msm(psor.msm)

### pmatrix.piecewise.msm
pmatrix.msm(psor.msm, 10)
times <- c(5, 10, 15)
covariates <- list(list(hieffusn=0, ollwsdrt=0),
                   list(hieffusn=0, ollwsdrt=1),
                   list(hieffusn=1, ollwsdrt=0),
                   list(hieffusn=1, ollwsdrt=1)
                   )
pmatrix.msm(psor.msm, 3, covariates=covariates[[1]])
pmatrix.piecewise.msm(psor.msm, 0, 3, times, covariates)
pmatrix.piecewise.msm(psor.msm, 0, 7, times, covariates)
pmatrix.piecewise.msm(psor.msm, 0, 19, times, covariates)
pmatrix.msm(psor.msm, 5, covariates[[1]]) %*% pmatrix.msm(psor.msm, 5, covariates[[2]]) %*% pmatrix.msm(psor.msm, 5, covariates[[3]]) %*% pmatrix.msm(psor.msm, 4, covariates[[4]])



#######  MISCELLANEOUS FEATURES  ##########

### Variables in global environment, not a data frame.
state.g <- heart$state; time.g <- heart$years; subj.g <- heart$PTNUM
heart.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = twoway4.q, fixedpars=TRUE)
heart.msm <- msm(state.g ~ time.g, subject=PTNUM, data=heart, qmatrix = twoway4.q, fixedpars=TRUE)
heart.msm

### Factor covariates

heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ factor(pdiag), covinits=list(sex=rep(0.1,7)), fixedpars=TRUE)
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ factor(pdiag), covinits=list("factor(pdiag)Nonexistentlevel"=rep(0.1,7)), fixedpars=TRUE)
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ factor(pdiag), covinits=list("factor(pdiag)Hyper"=rep(0.1,7)), fixedpars=TRUE) # OK 
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ factor(pdiag), covinits=list(sex=rep(0.1,7)), fixedpars=TRUE)
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ pdiag, covinits=list(pdiag=rep(0.1,7)), fixedpars=TRUE)
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ pdiag, covinits=list(pdiagNonexistentlevel=rep(0.1,7)), fixedpars=TRUE)
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q,
                       covariates = ~ pdiag, covinits=list(pdiagHyper=rep(0.1,7)), fixedpars=TRUE) # OK 
heartfaccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ pdiag,
                       covinits=list(pdiagHyper=rep(0.1,7),pdiagIDC=rep(0.1,7),pdiagIHD=rep(0.1,7),pdiagOther=rep(0.1,7),pdiagRestr=rep(0.1,7),), fixedpars=TRUE) # OK 
heartfaccov.msm


## Some data with censored states
## Replace first few death states by censorings 

heart.cens <- heart
heart.cens$state[heart$state==4][1:50] <- 99
heart.cens2 <- heart
heart.cens2$state[heart$state==4][1:50] <- 99
heart.cens2$state[heart$state==4][51:100] <- 999
heart.cens3 <- heart
ns <- c(heart$state[2:nrow(heart)], 0)
heart.cens3$state[heart$state==4][1:50] <- 99
heart.cens3$state[ns==4][1:50] <- 999

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
### Censored observations - final state censored
statetable.msm(state, PTNUM, heart.cens)
heartcens.msm <- msm(state ~ years, subject=PTNUM, data=heart.cens, qmatrix=twoway4.q, censor=99, fixedpars=TRUE)
heartcens.msm

### Censored observations - two kinds of censoring 
statetable.msm(state, PTNUM, heart.cens2)
heartcens.msm <- msm(state ~ years, subject=PTNUM, data=heart.cens2, qmatrix=twoway4.q, censor=c(99, 999), censor.states=list(c(1,2,3), c(2,3)), fixedpars=TRUE)
heartcens.msm

### Censored observations - intermediate state censored 
statetable.msm(state, PTNUM, heart.cens3)
heartcens.msm <- msm(state ~ years, subject=PTNUM, data=heart.cens3, qmatrix=twoway4.q, censor=c(99, 999), censor.states=list(c(2,3), c(1,2,3)), fixedpars=TRUE)
heartcens.msm

### crudeinits with censoring
try(crudeinits.msm(state~ years, PTNUM, twoway4.q, heart.cens))
crudeinits.msm(state~ years, PTNUM, twoway4.q, heart.cens, censor=99)
crudeinits.msm(state~ years, PTNUM, twoway4.q, heart.cens2, censor=c(99,999), censor.states=list(c(1,2,3),c(2,3)))
crudeinits.msm(state~ years, PTNUM, twoway4.q, heart.cens3, censor=c(99,999), censor.states=list(c(2,3),c(1,2,3)))

### Death with state at previous instant known - HIV model 
heart.dp <- heart
ns <- c(heart.dp$state[2:nrow(heart.dp)], 0)
heart.dp$years[ns==4][1:50] <- heart.dp$years[heart.dp$state==4][1:50]
### Observations at identical times not allowed in <= 0.4.1. 
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart.dp, qmatrix = twoway4.q, death = 4, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1))  # Works, with -2L of 2 less than baseline.
heart.msm

### Use "exacttimes" instead of "death" observation schemes for those obs of death with state at previous instant known
### Should be just the same.    
### Lik contrib for death, sum_r p(prev, r, t=0), q(r, death) =  q(prev, death)
### Lik contrib for exacttimes, p(prev, prev, t=0), q(prev, death) =  q(prev, death)
### since p(r,r,0) = 1. 
obstype <- rep(1, nrow(heart))
obstype[heart$state==4] <- 3
obstype[heart$state==4][1:50] <- 2
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart.dp, qmatrix = twoway4.q, obstype=obstype, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1))
heart.msm
obstype[heart$state==4][1:50] <- 3 # just to make sure this is the same 
heart.msm <- msm( state ~ years, subject=PTNUM, data = heart.dp, qmatrix = twoway4.q, obstype=obstype, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1))
heart.msm

