### TESTS OF MISCLASSIFICATION MODELS
### SPECIFIED USING THE OLD ematrix SYNTAX
source("local.R")

library(msm)
# library(msm, lib.loc="~/msm/src/1.2.7")
# data(cav) # not needed since 1.3

oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- rownames(ematrix) <- colnames(ematrix) <- c("Well","Mild","Severe","Death")

## Plain misc model with no covs

misc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE)
misc.msm
stopifnot(isTRUE(all.equal(4296.9155995778, misc.msm$minus2loglik, tol=1e-06)))

if (developer.local) {
    system.time(misc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                                qmatrix = oneway4.q, ematrix=ematrix, death = 4, # pci=5, # fixedpars=1:5,
                                control = list(trace=1, REPORT=1), method="BFGS"))
    stopifnot(isTRUE(all.equal(3951.82919869367, misc.msm$minus2loglik, tol=1e-06)))
    if(interactive()) save(misc.msm, file="~/msm/devel/models/misc.msm.rda")
}

## Covs on transition rates
misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars = TRUE,
                control = list(trace=1, REPORT=1), method="BFGS",
                covariates = ~ sex, covinits=list(sex=rep(0.1, 5)))
#stopifnot(isTRUE(all.equal(4299.35653620144, misccov.msm$minus2loglik, tol=1e-06))) # with last obs included in cov means, as pre 1.2.3
stopifnot(isTRUE(all.equal(4299.38058878142, misccov.msm$minus2loglik, tol=1e-06)))


## Covs on misc probs, old way.
misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                   qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=TRUE,
                   misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                   control = list(trace=1, REPORT=1), method="BFGS")
#stopifnot(isTRUE(all.equal(4304.90609473048, misccov.msm$minus2loglik, tol=1e-06))) # with last obs included in cov means
stopifnot(isTRUE(all.equal(4306.3077053482, misccov.msm$minus2loglik, tol=1e-06)))

if (developer.local) {
    system.time(misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                                   qmatrix = oneway4.q, ematrix=ematrix, death = 4,
                                   misccovariates = ~dage + sex, # fixedpars=c(1:5,11),
                                   control = list(trace=1, REPORT=1), method="BFGS"))
    ##    stopifnot(isTRUE(all.equal(3929.39438312539, misccov.msm$minus2loglik, tol=1e-06))) ## 0.7.1 and earlier
##    stopifnot(isTRUE(all.equal(3929.59504496140, misccov.msm$minus2loglik, tol=1e-06))) # 1.2.2 and earlier
    stopifnot(isTRUE(all.equal(3929.60136975058, misccov.msm$minus2loglik, tol=1e-06)))
    if(interactive()) save(misccov.msm, file="~/work/msm/devel/models/misccov.msm.rda")

    system.time(misccovboth.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=FALSE,
                                       control = list(trace=1, REPORT=1), method="BFGS",
                                       covariates = ~ sex, covinits=list(sex=rep(0.1, 5)),
                                       misccovariates = ~dage + sex, misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016))
                                       ))
    ##    stopifnot(isTRUE(all.equal(3921.40046811911, misccovboth.msm$minus2loglik, tol=1e-06))) # 0.7.1 and earlier
##     stopifnot(isTRUE(all.equal(3921.42240883417, misccovboth.msm$minus2loglik, tol=1e-06))) 1.2.2 and earlier
    stopifnot(isTRUE(all.equal(3921.42117997675, misccovboth.msm$minus2loglik, tol=1e-06)))
    if(interactive()) save(misccovboth.msm, file="~/work/msm/devel/models/misccovboth.msm.rda")
    if(interactive()) load("~/work/msm/devel/models/misccovboth.msm.rda")
    print(misccovboth.msm)

##########    OUTPUT FUNCTIONS    ###################

    if(interactive()) load("~/work/msm/devel/models/misc.msm.rda")
    if(interactive()) load("~/work/msm/devel/models/misccov.msm.rda")
    if(interactive()) load("~/work/msm/devel/models/misccovboth.msm.rda")

    e <- ematrix.msm(misc.msm)
    stopifnot(isTRUE(all.equal(0.00766164690017842, e$estimates[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.00334014173130528, e$SE[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.00325308687951399, e$L[1,2], tol=1e-06)))
    stopifnot(isTRUE(all.equal(0.0179371415700995, e$U[1,2], tol=1e-06)))

    print(ematrix.msm(misc.msm), digits=2)
    print(viterbi.msm(misc.msm)[1:50,])
    vit <- viterbi.msm(misc.msm)[viterbi.msm(misc.msm)$subject==100063,]
    stopifnot(isTRUE(all.equal(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2), vit$fitted, tol=1e-06)))

### TEST VITERBI ON SUBSET WITH FIXEDPARS BUG REPORT

    if(interactive()) load("~/work/msm/devel/models/misc.msm.rda")
    if(interactive()) load("~/work/msm/devel/models/misccov.msm.rda")
    pt <- 100046
    v <- viterbi.msm(misc.msm)
    for (pt in unique(cav$PTNUM)){
#        pt <- 100057
#        pt <- 100059
        subs <- cav[cav$PTNUM==pt,]
        x <- v[v$subject==pt,]
        miscfix.msm <- msm(state ~ years, subject = PTNUM, data = subs, qmatrix = qmatrix.msm(misc.msm)$estimates, ematrix=ematrix.msm(misc.msm)$estimates, death = 4, fixedpars=TRUE)
        y <- viterbi.msm(miscfix.msm)
        stopifnot(all(x$fitted==y$fitted))
    }

    odds <- odds.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.92024511282455, odds$dage[2,2], tol=1e-03))) # tests hadn't been updated since change in order in Version 0.9
    stopifnot(isTRUE(all.equal(0.937622692212605, odds$dage[1,2], tol=1e-03)))
    stopifnot(isTRUE(all.equal(23.6267954503319, odds$sex[2,3], tol=1e-01)))
    stopifnot(isTRUE(all.equal(30.7225710461079, odds$sex[4,3], tol=1e-01)))

    e <- ematrix.msm(misccov.msm) # 1.3 chenge
    stopifnot(isTRUE(all.equal(0.00425, e$estimates[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.0089, e$SE[1,2], tol=1e-02)))

    e <- ematrix.msm(misccov.msm, covariates=0)
    stopifnot(isTRUE(all.equal(0.00523189426375198, e$estimates[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.00717339683841909, e$SE[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.000182524613416697, e$L[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.06105006, e$U[1,2], tol=1e-02)))

    e <- ematrix.msm(misccov.msm, covariates=list(dage=50, sex=0))
    stopifnot(isTRUE(all.equal(0.0126402839581582, e$estimates[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.0127991364840901, e$SE[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.000570, e$L[1,2], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.031, e$U[1,2], tol=1e-02)))

### Non misclassification-specific output functions

    q <- qmatrix.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.234871828782083, q$estimates[2,3], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.0388940534715571, q$SE[2,3], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.169775324070530, q$L[2,3], tol=1e-02)))
    stopifnot(isTRUE(all.equal(0.324928114597637, q$U[2,3], tol=1e-02)))

    soj <- sojourn.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(c(6.83075273024829, 3.81966398452026, 3.29901891409971, 0.499352190044205,
0.41662507532632, 0.395976119914232, 5.91892340913934, 3.08447614101601,
2.60745249850692, 7.88305231146405, 4.73008455492069, 4.17400731243225
), as.numeric(unlist(soj)), tol=1e-02)))

    p <- pmatrix.msm(misccov.msm, 10)
    stopifnot(isTRUE(all.equal(0.122616549949547, p[1,3], tol=1e-03)))

    q <- qratio.msm(misccov.msm, c(1,2), c(2,3), cl=0.99)
    stopifnot(isTRUE(all.equal(c(0.449639898189716, 0.0948158687296418, 0.261198374155947, 0.774032528714554), as.numeric(q), tol=1e-02)))

    p <- prevalence.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(158, p$Observed[5,4], tol=1e-03)))
    stopifnot(isTRUE(all.equal(134.668020428293, p$Expected[5,4], tol=1e-03)))
    stopifnot(isTRUE(all.equal(31.43564356435644, p$"Observed percentages"[4,4], tol=1e-03)))
    stopifnot(isTRUE(all.equal(27.5796638007428, p$"Expected percentages"[4,4], tol=1e-03)))

    summ <- summary.msm(misccov.msm)
    p <- summ$prevalences
    stopifnot(isTRUE(all.equal(158, p$Observed[5,4], tol=1e-03)))
    stopifnot(isTRUE(all.equal(134.668020428293, p$Expected[5,4], tol=1e-03)))
    stopifnot(isTRUE(all.equal(31.43564356435644, p$"Observed percentages"[4,4], tol=1e-03)))
    stopifnot(isTRUE(all.equal(27.5796638007428, p$"Expected percentages"[4,4], tol=1e-03)))

    if (interactive()) plot.msm(misccov.msm)

    cf <- coef.msm(misccovboth.msm)
    stopifnot(isTRUE(all.equal(-0.533, cf$Qmatrices$sex[1,2], tol=1e-03)))
    stopifnot(isTRUE(all.equal(-7.83, cf$Ematrices$sex[1,2], tol=1e-03)))

    stopifnot(isTRUE(all.equal(c(1,2,3), as.numeric(transient.msm(misccov.msm)), tol=1e-03)))

    stopifnot(isTRUE(all.equal(4, as.numeric(absorbing.msm(misccov.msm)), tol=1e-03)))

    tot <- totlos.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(c(6.83075273082287, 2.75543070960308, 2.13503613252349), as.numeric(tot[1:3]), tol=1e-03)))

    stopifnot(isTRUE(all.equal(-1964.79752248070, as.numeric(logLik.msm(misccov.msm)), tol=1e-03)))


##########    OTHER FEATURES      ###################

    ## Baseline intens constraints
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=4:7,
                    qconstraint = c(1, 2, 1, 2, 3), analyticp=TRUE, 
                    control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(4209.65938095232, misc.msm$minus2loglik, tol=1e-06)))

    q <- qmatrix.msm(misc.msm)
    stopifnot(isTRUE(all.equal(-0.145819054714827, q$estimates[1,1], tol=1e-06)))

    ## Baseline misc constraints
    ematrix2 <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.11, 0),c(0, 0.11, 0, 0),c(0, 0, 0, 0))
    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                    qmatrix = oneway4.q, ematrix=ematrix2, death = 4, fixedpars=1:5,
                    econstraint = c(1, 1, 2, 2),
                    control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(4161.05767600322, misc.msm$minus2loglik, tol=1e-06)))
    e <- ematrix.msm(misc.msm)
    stopifnot(isTRUE(all.equal(0.970221499780077, e$estimates[1,1], tol=1e-06)))

    ## intens covariate constraints
    ## Give replicated inits for replicated cov effs, as that's consistent with constraints on q, e and h.

    misc.msm <- msm(state ~ years, subject = PTNUM, data = cav, fixedpars=c(1:5, 9:12),
                    qmatrix = oneway4.q, ematrix=ematrix, death = 4,
                    control = list(trace=1, REPORT=1), method="BFGS",
                    covariates = ~ sex, covinits=list(sex=c(0, 0, 0.1, 0, 0)),
                    constraint = list(sex = c(1, 2, 1, 2, 3))  )
    stopifnot(isTRUE(all.equal(4276.385, misc.msm$minus2loglik, tol=1e-06)))
#    stopifnot(isTRUE(all.equal(4277.77801412343, misc.msm$minus2loglik, tol=1e-06))) # pre 1.3
    q <- qmatrix.msm(misc.msm, covariates=0)
    stopifnot(isTRUE(all.equal(-0.176, q$estimates[1,1], tol=1e-02)))

    ## misc covariate constraints.
    misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars=c(1:5),
                       misccovariates = ~dage + sex,
                       misccovinits = list(dage=c(0.01,0.01,0.001,0.001), sex=c(0.0131,0.0132,0.0133,0.0134)),
                       miscconstraint = list(dage = c(1, 1, 2, 2)),
                       control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(4017.02699645160, misccov.msm$minus2loglik, tol=1e-06)))
    e <- ematrix.msm(misccov.msm)
    stopifnot(isTRUE(all.equal(0.9997, e$estimates[1,1], tol=1e-04)))
    e <- ematrix.msm(misccov.msm, covariates=0)
    stopifnot(isTRUE(all.equal(0.992813688955682, e$estimates[1,1], tol=1e-06)))

    ## fixedpars for misc covariates.  Parameters are ordered within covariate, within parameter, within state.
    ## memory errors with this
    misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = oneway4.q, ematrix=ematrix, death = 4,
                       misccovariates = ~dage + sex,
                       misccovinits = list(dage=c(0.01,0.02,0.03,0.04), sex=c(-0.013,-0.014,-0.015,-0.016)),
                       fixedpars = c(10, 11, 12, 15),
                       control = list(trace=1, REPORT=1), method="BFGS")
    stopifnot(isTRUE(all.equal(3946.96900738597, misccov.msm$minus2loglik, tol=1e-06)))

## multiple death states (Jean-Luc's data)
## Misclassification between states 2 and 3

  c2.df <- read.table("~/work/msm/tests/jeanluc/donneesaveccancerPT.txt", header=TRUE)
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

cat("misc.R: ALL TESTS PASSED\n")
