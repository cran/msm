### TESTS OF ERROR CHECKING
source("local.R")
print(developer.local)
library(msm)
data(heart)

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well","Mild","Severe","Death")
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
oneway4.i <- oneway4.q; oneway4.i[oneway4.i!=0] <- 1
ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))

### Rubbish in formula 

try(heart.msm <- msm(state, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))
try(heart.msm <- msm(~1, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))
try(heart.msm <- msm("foo", subject=PTNUM, data = heart, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))

### Rubbish in qmatrix 
wrong.q <- cbind(c(0,1,2), c(0,1,2))
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, death = TRUE, fixedpars=TRUE))
wrong.q <- cbind(c(0,1), c(0,1))
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, death = TRUE, fixedpars=TRUE))
wrong.q <- "foo"
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, death = TRUE, fixedpars=TRUE))
wrong.q <- 1
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, death = TRUE, fixedpars=TRUE))

### Rubbish in ematrix 
wrong.e <- "foo"
try(misc.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = oneway4.q, ematrix=wrong.e, death = 4, fixedpars=TRUE))
wrong.e <- 1
try(misc.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = oneway4.q, ematrix=wrong.e, death = 4, fixedpars=TRUE))
wrong.e <- cbind(c(0,1,2), c(0,1,2))
try(misc.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = oneway4.q, ematrix=wrong.e, death = 4, fixedpars=TRUE))
wrong.e <- cbind(c(0,1), c(0,2))
try(misc.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = oneway4.q, ematrix=wrong.e, death = 4, fixedpars=TRUE))

### Rubbish in subject
try(heart.msm <- msm(state~years, subject="foo", data = heart, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))
try(heart.msm <- msm(state~years, subject=foo, data = heart, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))

### Rubbish in obstype 
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype="foo", death = TRUE, fixedpars=TRUE))
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype=rep(1,10), death = TRUE, fixedpars=TRUE))
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype=rep(4, nrow(heart)), death = TRUE, fixedpars=TRUE))

### covariates not in data
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = "wibble"))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sux))

### misccovariates not in data
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=TRUE, misccovariates = "wobble"))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sox))

### covinits not in data, wrong length, rubbish 
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits="foo", fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits=list(sex="foo", age="bar"), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits=list(sex=c(1,2,3), age="bar"), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits=list(age=rep(0.1, 7)), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits=list(age=rep(0.1, 7), foo=1, bar=2), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits=list(sex=1), fixedpars=TRUE))
### misccovinits not in data, wrong length, rubbish 
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, misccovinits="foo", fixedpars=TRUE))
if (developer.local) { 
    try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, misccovinits=list(sex="foo", age="bar"), fixedpars=TRUE)) # covariates but misccovinits,  no error
    heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, misccovinits=list(sex="foo", age="bar"), fixedpars=TRUE)
}
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, misccovariates = ~ sex, misccovinits=list(sex=1, age="bar"), fixedpars=TRUE))
## misccovs specified but misc false 
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, misccovariates = ~ sex, misccovinits=list(sex=1, age="bar"), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, ematrix=ematrix, misccovariates = ~ sex, misccovinits=list(sex=1, age="bar"), fixedpars=TRUE))

### constraint not in data, wrong length, rubbish 
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, constraint="foo"))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, constraint=list(foo="bar")))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, constraint=list(foo=1)))

### miscconstraint
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, constraint="foo", fixedpars=TRUE)) # constraint but no covs 
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, miscconstraint=list(foo="bar")))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, miscconstraint=list(foo=1)))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, miscconstraint=list(sex=1)))

### initprobs
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, initprobs="poo", fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, initprobs=c(1,2), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, initprobs=c(2,1,1,1), fixedpars=TRUE)) # correct - scaled to sum to 1.

### qconstraint
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, qconstraint="foo", fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, qconstraint=list(c(1,1,2)), fixedpars=TRUE))
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, qconstraint=c(1,1,2), fixedpars=TRUE))

### econstraint
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, econstraint="foo", fixedpars=TRUE) )
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, econstraint=list(c(1,1,2)), fixedpars=TRUE) )
try(heart.msm <- msm( state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, econstraint=c(1,1,2), fixedpars=TRUE))

### States in data not in state set specified with qmatrix 
wrong.q <- cbind(c(0,1), c(0,1)) # extra states in data 
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, fixedpars=TRUE))
wrong.q <- rbind(c(0,1,2,3,1), c(0,1,3,4,1), c(0,1,2,3,2), c(0,1,2,3,4), c(0,0,0,0,0))
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, fixedpars=TRUE)) # too few states in data, warning only 
wrong.q <- rbind(c(0,1,2,3,1,0), c(0,1,3,4,1,0), c(0,1,2,3,2,0), c(0,1,2,3,4,0), c(0,0,0,0,0,0), c(0,0,0,0,0,0))
try(heart.msm <- msm(state~years, subject=PTNUM, data = heart, qmatrix = wrong.q, fixedpars=TRUE))

### Observations not increasing with time
heart.wrong <- heart
heart.wrong$years[3:5] <- 4:2
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))

### Observations from a subject not all adjacent
heart.wrong <- heart
heart.wrong$PTNUM[4:5] <- 100003
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE))

### Data inconsistent with P matrix (snapshot) 
heart.wrong <- heart
heart.wrong$state[4] <- 1
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = oneway4.q, death = TRUE, fixedpars=TRUE))
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.i, gen.inits=TRUE, fixedpars=TRUE))

heart.cens <- heart
heart.cens$state[heart$state==4][1:50] <- 99
### censored states inconsistent (currently gives Inf) 
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.cens, qmatrix = oneway4.q, censor=99, censor.states=list(c(1,2)), fixedpars=TRUE))
heart.msm
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.cens, qmatrix = oneway4.q, censor=99, censor.states=list(c(1,2,3)), fixedpars=TRUE))
heart.msm
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.cens, qmatrix = oneway4.q, censor=99, censor.states=list(c(1,2)), exacttimes=TRUE, fixedpars=TRUE))
heart.msm

### Data inconsistent with Q matrix (exacttimes) 
heart.wrong <- heart
heart.wrong$state[4] <- 1
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = twoway4.q, exacttimes=TRUE, fixedpars=TRUE))
obstype <- rep(2, nrow(heart))
obstype[10] <- 1
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = twoway4.q, obstype=obstype, fixedpars=TRUE))
heart.msm

### rubbish in death
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = "foo", fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = 5, fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = 1:5, fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, death = 3, fixedpars=TRUE) ) # non-absorbing death state

### rubbish in censor 
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, censor="rubbish", fixedpars=TRUE))
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, censor=1, fixedpars=TRUE)) # warning, censor=actual state
# heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, censor=1, fixedpars=FALSE, method="BFGS", control=list(trace=5, REPORT=1)) # works but extreme estimates 
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, censor.states="rubbish", fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.cens, qmatrix = twoway4.q, censor=99, censor.states="rubbish", fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.cens, qmatrix = twoway4.q, censor=99, censor.states=list(c(1,2,3), "rubbish"), fixedpars=TRUE) )

### obstype
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype="rubbish", fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype=c(1,2,3), fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype=4, fixedpars=TRUE) )
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, obstype=rep(4,nrow(heart)), fixedpars=TRUE) )

### fixedpars
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars="foo"))
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=list(c(1,3,4))))
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=1:8))
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=0.5))

### crudeinits
heart.wrong <- heart
heart.wrong$state[4] <- 1
crudeinits.msm(state ~ years, PTNUM, twoway4.q, heart.wrong) # no error. ignores inconsistent transitions 
crudeinits.msm(state ~ years, PTNUM, oneway4.q, heart) # no error. ignores inconsistent transitions 
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = twoway4.q, exacttimes=TRUE, fixedpars=TRUE))
obstype <- rep(2, nrow(heart))
obstype[10] <- 1
try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart.wrong, qmatrix = twoway4.q, obstype=obstype, fixedpars=TRUE))
heart.msm


#### OUTPUTS #####

try(heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=TRUE))
try(heartfit.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q))
try(heartcov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, covariates = ~ sex, covinits=list(sex=rep(0.1,7)), fixedpars=TRUE))
try(heartmisc.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, fixedpars=TRUE))
try(heartmisccov.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = oneway4.q, ematrix=ematrix, misccovariates = ~ sex, misccovinits=list(sex=rep(0.1,4)), fixedpars=TRUE))

### qematrix 
try(qmatrix.msm("foo"))
try(qmatrix.msm(heart.msm))
try(qmatrix.msm(heart.msm, covariates="foo"))
try(qmatrix.msm(heart.msm, covariates=list(foo=1)))
try(qmatrix.msm(heartcov.msm, covariates=list(foo=1)))

### ematrix 
try(ematrix.msm("Foo"))
try(ematrix.msm(heart.msm)) # not a misc model, returns NULL.
try(ematrix.msm(heartmisc.msm, covariates="foo"))
try(ematrix.msm(heartmisc.msm, covariates=list(foo=1)))
try(ematrix.msm(heartmisccov.msm, covariates=list(foo=1)))

### sojourn
try(sojourn.msm("foo"))
try(sojourn.msm(heart.msm, covariates="foo"))
try(sojourn.msm(heart.msm, covariates=list(foo=1)))

### pmatrix
try(pmatrix.msm("foo"))
try(pmatrix.msm(heart.msm, t="foo"))
try(pmatrix.msm(heart.msm, -9))
try(pmatrix.msm(heart.msm, 0))
try(pmatrix.msm(heart.msm, 1, covariates=list(foo=1)))

### pmatrix.piecewise
try(pmatrix.piecewise.msm("foo", 1,2, c(1, 2), c(1,2)))
try(pmatrix.piecewise.msm("foo", 1, 2, c(1, 0.5, 2), list(0, 1, 0, 1)))
try(pmatrix.piecewise.msm(heart.msm, 1, 2, c(1, 0.5, 2), list(0, 1, 0, 1)))
try(pmatrix.piecewise.msm(heart.msm, 1, 2, "rubbish", list(0, 1, 0, 1)))
try(pmatrix.piecewise.msm(heartcov.msm, 1, 2, c(1, 1.5, 2), "rubbish"))
try(pmatrix.piecewise.msm(heartcov.msm, 1, 2, c(1, 1.5, 2), list("rubbish", "foo","bar","boing")))
try(pmatrix.piecewise.msm(heartcov.msm, 1, 2, c(1, 1.5, 2), list(0, 1, 0, 1)))
pmatrix.piecewise.msm(heartcov.msm, 1, 2, c(1, 1.5, 2), list(list(0), list(1), list(0), list(1))) # OK
pmatrix.msm(heartcov.msm, 1) # OK

### qratio
try(qratio.msm("foo"))
try(qratio.msm(heart.msm, "foo"))
try(qratio.msm(heart.msm, c(1,8), c(1,0)))
try(qratio.msm(heart.msm, c(1,2), c(1,0)))
qratio.msm(heart.msm, c(1,2), c(2,3)) # OK
qratio.msm(heartfit.msm, c(1,2), c(2,3)) # OK 
try(qratio.msm(heartfit.msm, c(1,2), c(2,3), cl="foo"))
try(qratio.msm(heartfit.msm, c(1,2), c(2,3), cl=2))

### summary, prevalence, observed
try(prevalence.msm("foo"))
try(summary.msm("foo"))

### plot
if (interactive()) { 
    try(plot.msm("foo"))
    try(plot.msm(heart.msm, from="foo"))
    try(plot.msm(heart.msm, to="foo"))
    try(plot.msm(heart.msm, from = 1:8, to=3))
    try(plot.msm(heart.msm, to = 3))
    try(plot.msm(heart.msm, range="foo"))
    try(plot.msm(heart.msm, range=1:6))
    try(plot.msm(heart.msm))
    try(plot.msm(heartfit.msm))
}

### coef, hazard, odds
try(coef.msm("foo"))
try(hazard.msm("foo"))
try(hazard.msm(heartcov.msm, hazard.scale="foo"))
try(hazard.msm(heartcov.msm))
try(hazard.msm(heartcov.msm, hazard.scale=c(1,2,3,4)))
try(hazard.msm(heartcov.msm, hazard.scale=c(1,2,3)))
try(hazard.msm(heartcov.msm, hazard.scale=2))

try(odds.msm("foo"))
try(odds.msm(heartcov.msm, odds.scale="foo"))
try(odds.msm(heartmisccov.msm, odds.scale="foo"))
odds.msm(heartmisccov.msm) # OK
try(odds.msm(heartmisccov.msm, odds.scale=c(1,2,3,4)))
try(odds.msm(heartmisccov.msm, odds.scale=c(1,2,3)))
odds.msm(heartmisccov.msm, odds.scale=2) # OK 

### transient, absorbing
try(transient.msm("foo"))
try(absorbing.msm("foo"))
try(transient.msm(qmatrix="foo"))
try(absorbing.msm(qmatrix="foo"))
try(transient.msm(qmatrix=c(1,4,5,6)))
try(absorbing.msm(qmatrix=c(1,4,5,6)))
try(transient.msm(qmatrix=cbind(c(1,2,3),c(1,3,2))))
try(absorbing.msm(qmatrix=cbind(c(1,2,3),c(1,3,2))))
try(transient.msm())
try(absorbing.msm())

### totlos
try(totlos.msm("foo"))
try(totlos.msm(heart.msm, start="foo"))
try(totlos.msm(heart.msm, start=-1))
try(totlos.msm(heart.msm, start=1, fromt=1, tot=0))
try(totlos.msm(heart.msm, start=1, fromt=c(1,2), tot=c(3,4)))
try(totlos.msm(heart.msm, start=1, fromt="foo", tot=2))
try(totlos.msm(heart.msm, start=1, fromt=-3, tot=-2))
totlos.msm(heart.msm, start=1, fromt=1, tot=2) # OK

### logLik
try(logLik.msm("foo"))
logLik.msm(heart.msm)

### viterbi
try(viterbi.msm("foo"))
