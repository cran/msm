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
stopifnot(isTRUE(all.equal(4296.9155995778, misc.msm$minus2loglik, tol=1e-06)))

if (developer.local) { 
    system.time(misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                                qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                                control = list(trace=1, REPORT=1), method="BFGS")) # new 43, old 129
    stopifnot(isTRUE(all.equal(3951.82919869367, misc.msm$minus2loglik, tol=1e-06)))
    if(interactive()) save(misc.msm, file="~/msm/devel/models/misc.msm.rda")
}

## Does misc model with no misc reduce to simple 
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
nomisc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), death = 4, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4908.81676837903, nomisc.msm$minus2loglik, tol=1e-06)))
simple.msm <- msm(state ~ years, subject = PTNUM, data = heart, qmatrix = twoway4.q, death = 4, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4908.81676837903, simple.msm$minus2loglik, tol=1e-06)))

## Covs on transition rates
misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars = TRUE, 
                control = list(trace=1, REPORT=1), method="BFGS",
                covariates = ~ sex, covinits=list(sex=rep(0.1, 5)))
stopifnot(isTRUE(all.equal(4299.35653620144, misccov.msm$minus2loglik, tol=1e-06)))

## Covs on misc probs, old way. 
misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                   qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE,
                   misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                   control = list(trace=1, REPORT=1), method="BFGS")
stopifnot(isTRUE(all.equal(4306.82007050922, misccov.msm$minus2loglik, tol=1e-06)))

if (developer.local) { 
    system.time(misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                                   qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                                   misccovariates = ~dage + sex,
                                   control = list(trace=1, REPORT=1), method="BFGS"))
    stopifnot(isTRUE(all.equal(3929.39438312539, misccov.msm$minus2loglik, tol=1e-06)))
    if(interactive()) save(misccov.msm, file="~/msm/devel/models/misccov.msm.rda")
}

## Covs on both
misccovboth.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE, 
                       control = list(trace=1, REPORT=1), method="BFGS",
                       covariates = ~ sex, covinits=list(sex=rep(0.1, 5)),
                       misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016))
                       )
stopifnot(isTRUE(all.equal(4309.26368021750, misccovboth.msm$minus2loglik, tol=1e-06))) ## warnings here: 15: is.na() applied to non-(list or vector) in: is.na(i) 16: no non-missing arguments to max; returning -Inf. between first two msm.form.covdata's? 
#library(msm, lib.loc="d:/work/msm/lib/windows/0.7")


if (developer.local) { 
    system.time(misccovboth.msm <- msm(state ~ years, subject = PTNUM, data = heart, 
                                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=FALSE, 
                                       control = list(trace=1, REPORT=1), method="BFGS",
                                       covariates = ~ sex, covinits=list(sex=rep(0.1, 5)),
                                       misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016))
                                       ))
    stopifnot(isTRUE(all.equal(3921.40046811911, misccovboth.msm$minus2loglik, tol=1e-06)))
    if(interactive()) save(misccovboth.msm, file="~/msm/devel/models/misccovboth.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccovboth.msm.rda")
    print(misccovboth.msm)

##########    OUTPUT FUNCTIONS    ###################

    if(interactive()) load("~/msm/devel/models/misc.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccov.msm.rda")
    if(interactive()) load("~/msm/devel/models/misccovboth.msm.rda")

    e <- ematrix.msm(misc.msm)
    stopifnot(isTRUE(all.equal(0.00766195389341105, e$estimates[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.00334013763683592, e$SE[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.00325333223799247, e$L[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.0179372312008091, e$U[1,2], tol=1e-06)))
    
    print(ematrix.msm(misc.msm), digits=2)
    print(viterbi.msm(misc.msm)[1:50,])
    vit <- viterbi.msm(misc.msm)[viterbi.msm(misc.msm)$subject==100063,]
    stopifnot(isTRUE(all.equal(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2), vit$fitted, tol=1e-06)))

    odds <- odds.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.924920277547759, odds$dage[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.9350691888108385, odds$dage[2,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(24.76686951404301, odds$sex[1,3], tol=1e-04)))
    stopifnot(isTRUE(all.equal(30.90173037227264, odds$sex[3,3], tol=1e-04)))

    e <- ematrix.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.003848160553285282, e$estimates[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.01036733815389283, e$SE[1,2], tol=1e-03)))

    e <- ematrix.msm(misccov.msm, covariates=0)
    stopifnot(isTRUE(all.equal(0.005168493472109885, e$estimates[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.007347868833274329, e$SE[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.0003155487782038569, e$L[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.07877543777698336, e$U[1,2], tol=1e-04)))

    e <- ematrix.msm(misccov.msm, covariates=list(dage=50, sex=0))  
    stopifnot(isTRUE(all.equal(0.0119933869179702, e$estimates[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.01233677079484181, e$SE[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.001575057755385304, e$L[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.08542810558239693, e$U[1,2], tol=1e-04)))

### Non misclassification-specific output functions 

    q <- qmatrix.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.2359981283678808, q$estimates[2,3], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.03911018205606053, q$SE[2,3], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.1705475068597598, q$L[2,3], tol=1e-04)))
    stopifnot(isTRUE(all.equal(0.3265665832273966, q$U[2,3], tol=1e-04)))

    soj <- sojourn.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(c(6.7952181814242, 3.82452445494119, 3.30260275866043, 0.507093246772922, 
0.418641351699219, 0.395262006422924, 5.87059985108007, 3.08604986663779, 
2.61205827140167, 7.86546371827108, 4.73971190957426, 4.17570507554508
), as.numeric(unlist(soj)), tol=1e-06)))
    
    p <- pmatrix.msm(misccov.msm, 10)
    stopifnot(isTRUE(all.equal(0.1238626028171057, p[1,3], tol=1e-06)))
    
    q <- qratio.msm(misccov.msm, c(1,2), c(2,3), cl=0.99)
    stopifnot(isTRUE(all.equal(c(0.450437513842305, 0.0950161215095277, 0.261613747651346, 0.775547751973777), as.numeric(q), tol=1e-04)))

    p <- prevalence.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(158, p$Observed[5,4], tol=1e-06)))
    stopifnot(isTRUE(all.equal(134.6193887170075, p$Expected[5,4], tol=1e-06)))
    stopifnot(isTRUE(all.equal(31.43564356435644, p$"Observed percentages"[4,4], tol=1e-06)))
    stopifnot(isTRUE(all.equal(27.55723041956071, p$"Expected percentages"[4,4], tol=1e-06)))

    summ <- summary.msm(misccov.msm)
    p <- summ$prevalences
    stopifnot(isTRUE(all.equal(158, p$Observed[5,4], tol=1e-06)))
    stopifnot(isTRUE(all.equal(134.6193887170075, p$Expected[5,4], tol=1e-06)))
    stopifnot(isTRUE(all.equal(31.43564356435644, p$"Observed percentages"[4,4], tol=1e-06)))
    stopifnot(isTRUE(all.equal(27.55723041956071, p$"Expected percentages"[4,4], tol=1e-06)))

    if (interactive()) plot.msm(misccov.msm)

    cf <- coef.msm(misccovboth.msm)
    stopifnot(isTRUE(all.equal(-0.531938307861184, cf$Qmatrices$sex[1,2], tol=1e-04)))
    stopifnot(isTRUE(all.equal(-7.6573606063652, cf$Ematrices$sex[1,2], tol=1e-04)))

    stopifnot(isTRUE(all.equal(c(1,2,3), transient.msm(misccov.msm), tol=1e-06)))
    
    stopifnot(isTRUE(all.equal(4, absorbing.msm(misccov.msm), tol=1e-06)))

    tot <- totlos.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(c(6.79521818195318, 2.76263786152668, 2.15322224353698), as.numeric(tot), tol=1e-06)))

    stopifnot(isTRUE(all.equal(1964.697191562695, as.numeric(logLik.msm(misccov.msm)), tol=1e-06)))

}


##########    OTHER FEATURES      ###################

if (developer.local) { 
    ## Baseline intens constraints
    misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=4:7,
                    qconstraint = c(1, 2, 1, 2, 3), 
                    control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(4209.65938095232, misc.msm$minus2loglik, tol=1e-06)))
    q <- qmatrix.msm(misc.msm)
    stopifnot(isTRUE(all.equal(-0.145819054714827, q$estimates[1,1], tol=1e-06)))

    ## Baseline misc constraints
    ematrix2 <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.11, 0),c(0, 0.11, 0, 0),c(0, 0, 0, 0))
    misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                    qmatrix = oneway4.q, ematrix=ematrix2, death = 4, fixedpars=1:5,
                    econstraint = c(1, 1, 2, 2), 
                    control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(4163.71591850274, misc.msm$minus2loglik, tol=1e-06)))
    e <- ematrix.msm(misc.msm)
    stopifnot(isTRUE(all.equal(0.96982246176299, e$estimates[1,1], tol=1e-06)))

    ## intens covariate constraints
    ## Give replicated inits for replicated cov effs, as that's consistent with constraints on q, e and h.

    misc.msm <- msm(state ~ years, subject = PTNUM, data = heart, fixedpars=c(1:5, 9:12), 
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                    control = list(trace=1, REPORT=1), method="BFGS",
                    covariates = ~ sex, covinits=list(sex=c(0, 0, 0.1, 0, 0)), 
                    constraint = list(sex = c(1, 2, 1, 2, 3))  )
    stopifnot(isTRUE(all.equal(4277.77801412343, misc.msm$minus2loglik, tol=1e-06)))
    q <- qmatrix.msm(misc.msm, covariates=0)
    stopifnot(isTRUE(all.equal(-0.17619654476216, q$estimates[1,1], tol=1e-06)))

    ## misc covariate constraints. Version 0.4.1 breaks when calculating SEs here (singular in dgesv).

    misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=c(1:5), 
                       misccovariates = ~dage + sex,
                       misccovinits = list(dage=c(0.01,0.01,0.001,0.001), sex=c(0.0131,0.0132,0.0133,0.0134)),
                       miscconstraint = list(dage = c(1, 1, 2, 2)), 
                       control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(4013.68740622917, misccov.msm$minus2loglik, tol=1e-06)))
    e <- ematrix.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.999568228582394, e$estimates[1,1], tol=1e-06)))
    e <- ematrix.msm(misccov.msm, covariates=0)
    stopifnot(isTRUE(all.equal(0.993214215656684, e$estimates[1,1], tol=1e-06)))
    
    ## fixedpars for misc covariates.  Parameters are ordered within covariate, within parameter, within state.

    misccov.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, 
                       misccovariates = ~dage + sex,
                       misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                       fixedpars = c(10, 11, 12, 15),
                       control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(3945.10580606042, misccov.msm$minus2loglik, tol=1e-06)))
}

## multiple death states (Jean-Luc's data)
## Misclassification between states 2 and 3 

if (developer.local) { 
  c2.df <- read.table("~/msm/tests/jeanluc/donneesaveccancerPT.txt", header=TRUE)
  print(statetable.msm(state, PTNUM, c2.df))
  qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
  ex <- rbind( c(0, 0, 0, 0, 0), c(0, 0, 0.1, 0, 0), c(0, 0.1, 0, 0, 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0) ) 
  c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx, ematrix=ex, death=c(4, 5), method="BFGS", fixedpars = TRUE)
  stopifnot(isTRUE(all.equal(70084.3665626129, c2.msm$minus2loglik, tol=1e-06)))
  
  ## multiple death states specified using an obstype vector 
  d45 <- rep(1, nrow(c2.df)); d45[c2.df$state %in% c(4,5)] <- 3
  c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx, ematrix=ex, obstype=d45, method="BFGS", fixedpars = TRUE)
  stopifnot(isTRUE(all.equal(70084.3665626129, c2.msm$minus2loglik, tol=1e-06)))
}

## exact times
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, exacttimes=TRUE, fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS") # should warn about redundant death argument
stopifnot(isTRUE(all.equal(4864.14764195147, misc.msm$minus2loglik, tol=1e-06)))

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, exacttimes=TRUE, fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS")
stopifnot(isTRUE(all.equal(4864.14764195147, misc.msm$minus2loglik, tol=1e-06)))

## exact times specified using an obstype vector
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, obstype=rep(2, nrow(heart)), fixedpars=TRUE,
                control = list(trace=1, REPORT=1), method="BFGS")
stopifnot(isTRUE(all.equal(4864.14764195147, misc.msm$minus2loglik, tol=1e-06)))

## initprobs
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE, initprobs=c(0.7, 0.1, 0.1, 0.1),
                control = list(trace=1, REPORT=1), method="BFGS")
stopifnot(isTRUE(all.equal(4725.9078185031, misc.msm$minus2loglik, tol=1e-06)))

## initprobs in Viterbi : bug fix for 0.6.1
if (developer.local) { 
    miscinitp.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, initprobs=c(0.6, 0.4, 0, 0), 
                    control = list(trace=1, REPORT=1), method="BFGS")
    if(interactive()) save(miscinitp.msm, file="~/msm/devel/models/miscinitp.msm.rda")
    if(interactive()) load(file="~/msm/devel/models/miscinitp.msm.rda")
    vitinitp <- viterbi.msm(miscinitp.msm)
    table(vitinitp$fitted[vitinitp$time==0]) / nrow(vitinitp[vitinitp$time==0,])
}

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
stopifnot(isTRUE(all.equal(4025.42265024404, misc.msm$minus2loglik, tol=1e-06)))

## Two types of censoring
misc.msm <- msm(state ~ years, subject=PTNUM, data=heart.cens2, qmatrix=oneway4.q, ematrix=ematrix, censor=c(99, 999), death=4, censor.states=list(c(1,2,3), c(2,3)), fixedpars=TRUE)
stopifnot(isTRUE(all.equal(3822.04540210944, misc.msm$minus2loglik, tol=1e-06)))

## Does misc model with no misc reduce to simple, with censoring 

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
misc.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens,
                qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), censor=99, death=TRUE, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4759.28151596975, misc.msm$minus2loglik, tol=1e-06)))

simple.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens, qmatrix = twoway4.q, death=TRUE, censor=99, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4759.28151596975, simple.msm$minus2loglik, tol=1e-06)))

misc.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens,
                qmatrix = twoway4.q, ematrix=matrix(0, nrow=4, ncol=4), censor=99, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4724.26606344485, misc.msm$minus2loglik, tol=1e-06)))

simple.msm <- msm(state ~ years, subject = PTNUM, data = heart.cens, qmatrix = twoway4.q, censor=99, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4724.26606344485, simple.msm$minus2loglik, tol=1e-06)))


## Viterbi with non-HMM model
twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
heart.msm <- msm(state ~ years, subject=PTNUM, data = heart, qmatrix = twoway4.q, fixedpars=TRUE)
vit <- viterbi.msm(heart.msm)[1:50,] # no error, returns observed states. 
stopifnot(all.equal(vit$observed, vit$fitted))


#### Estimating initprobs

if (developer.local)
  misc.msm <- msm(state ~ years, subject = PTNUM, data = heart,
                  qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=FALSE, initprobs=rep(0.25, 4), est.initprobs=TRUE, 
                  control = list(trace=1, REPORT=1), method="BFGS")
### just converges to 1,0,0,0

#### Simulate data with known initprobs. 
nsubj <- 50; nobspt <- 6
sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt), 
                     x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 2 )
(three.q <- msm:::msm.fixdiag.qmatrix(rbind(c(0, exp(-3), exp(-6)), c(0, 0, exp(-3)), c(0, 0, 0))))
ematrix3 <- rbind(c(0, 0.1, 0), c(0.1, 0, 0), c(0,0,0))
initprobs <- c(0.5, 0.5, 0)
set.seed(22061976)
sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, ematrix=ematrix3, start=sample(1:3, 50, prob=initprobs, replace=TRUE))
misc.msm <- msm(obs ~ time, subject = subject, data = sim2.df,
                qmatrix = three.q, ematrix=ematrix3, initprobs=c(0.1, 0.9, 0), fixedpars=7, est.initprobs=TRUE, 
                control = list(trace=1, REPORT=1), method="BFGS")
stopifnot(misc.msm$hmodel$initprobs["State 2","L95"] < 0.5 && 0.5 < misc.msm$hmodel$initprobs["State 2","U95"])

cat("misc.R: ALL TESTS PASSED\n")



