source("local.R")
library(msm)
data(fev)

### TESTS FOR HIDDEN MARKOV MODELS WITH CONTINUOUS RESPONSES
### USING LUNG TRANSPLANT FEV1 DATA

three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))

four.q <-  rbind(c(0, exp(-6), 0, exp(-9)), c(0, 0, exp(-6.01), exp(-9)), c(0, 0, 0, exp(-6.02)), c(0, 0, 0, 0))

five.q <-  rbind(c(0, exp(-6), 0, 0, exp(-9)),
                 c(0, 0, exp(-6.01), 0, exp(-9)),
                 c(0, 0, 0, exp(-6.02), exp(-6.03)),
                 c(0, 0, 0, 0, exp(-6.04)),
                 c(0, 0, 0, 0, 0))

### Test against old C code 
## Old program on same parameters. 
## When using covariates, likelihoods match when covs are not centred on their means.

oldlik <- function(cmdfile) {
    cmdfile <- paste("../../hid/newcmdfiles/", cmdfile, ".cmd", sep="")
    cmd <- paste("../../hid/hidden", cmdfile)
    output <- readLines(pipe(cmd))
    lik <- output[grep("-2\\*loglikelihood =", output)]
    lik <- as.numeric(gsub("-2\\*loglikelihood =", "", lik))
    lik
}

test.oldlik <- function(msmmod, cmdfile, tol=1e-04) {
    old <- oldlik(cmdfile)
    likdiff <- abs (old - msmmod$minus2loglik)
    if (likdiff < tol)
      cat("PASSED. Likelihood matches old, difference = ", likdiff, "\n")
    else cat("FAILED. Old lik ", old, ", new ", msmmod$minus2loglik, ", difference = ", likdiff, "\n")
    invisible()
}

## One-level model

hmodel3 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
hmodel4 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=72.5, sd=10), hmmNorm(mean=42.5, sd=18), hmmIdent(999))
hmodel5 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=72.5, sd=10), hmmNorm(mean=62.5, sd=10), hmmNorm(mean=42.5, sd=18), hmmIdent(999))

(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
if (developer.local) test.oldlik(fev3.hid, "onelevel3")

(fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, death=4, hmodel=hmodel4, fixedpars=TRUE))
if (developer.local) test.oldlik(fev4.hid, "onelevel4") ## old program on same parameters.

(fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, death=5, hmodel=hmodel5, fixedpars=TRUE))
if (developer.local) test.oldlik(fev5.hid, "onelevel5") ## old program on same parameters.

## One-level model with covariate on response - indicator for acute events within 14 days 

(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3,
                 hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                 hconstraint = list(acute = c(1,1)), 
                 fixedpars=TRUE, center=FALSE))
if (developer.local) test.oldlik(fev3.hid, "onelevel3cov")

(fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, death=4, hmodel=hmodel4,
                 hcovariates=list(~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, NULL),
                 hconstraint = list(acute = c(1,1,1)), 
                 fixedpars=TRUE, center=FALSE))
if (developer.local) test.oldlik(fev4.hid, "onelevel4cov")

(fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, death=5, hmodel=hmodel5,
                 hcovariates=list(~acute, ~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, -8, NULL),
                 hconstraint = list(acute = c(1,1,1,1)), 
                 fixedpars=TRUE, center=FALSE))
if (developer.local) test.oldlik(fev5.hid, "onelevel5cov")

## Two-level Satten and Longini model 

hmodel3 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                hmmIdent(999))
hmodel4 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                hmmMEUnif(sderr=8, lower=65, upper=80, meanerr=0),
                hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=65, meanerr=0),
                hmmIdent(999))
hmodel5 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                hmmMEUnif(sderr=8, lower=65, upper=80, meanerr=0),
                hmmMEUnif(sderr=8, lower=50, upper=65, meanerr=0),
                hmmMETNorm(mean=42, sd=18, sderr=8, lower=0, upper=50, meanerr=0),
                hmmIdent(999))

(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3, fixedpars=TRUE))
if (developer.local) test.oldlik(fev3.hid, "twolevel3")

(fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, death=4, hmodel=hmodel4, fixedpars=TRUE))
if (developer.local) test.oldlik(fev4.hid, "twolevel4")

(fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, death=5, hmodel=hmodel5, fixedpars=TRUE))
if (developer.local) test.oldlik(fev5.hid, "twolevel5")


## Two-level Satten and Longini model with acute event covariate on the measurement error model

(fev3.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3,
                 hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                 fixedpars=TRUE, center=FALSE))
if (developer.local) test.oldlik(fev3.hid, "twolevel3cov")

(fev4.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=four.q, death=4, hmodel=hmodel4,
                 hcovariates=list(~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, NULL),
                 fixedpars=TRUE, center=FALSE))
if (developer.local) test.oldlik(fev4.hid, "twolevel4cov")

(fev5.hid <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=five.q, death=5, hmodel=hmodel5,
                 hcovariates=list(~acute, ~acute, ~acute, ~acute, NULL), hcovinits = list(-8, -8, -8, -8, NULL),
                 fixedpars=TRUE, center=FALSE))
if (developer.local) test.oldlik(fev5.hid, "twolevel5cov")


#### MODELS ACTUALLY FITTED. These are the two models presented in the manual. 

### On some platforms (not bumblebee) doesn't converge (no SEs) with days - use months instead. 

if (developer.local) { 

    three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))
    months <- fev$days * 12 / 365
    hmodel1 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))
    (fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel1,
                     hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                     hconstraint = list(acute = c(1,1)), 
                     fixedpars=FALSE, center=FALSE, method="BFGS", control=list(trace=1, REPORT=1)))
    print(fev1.msm)
    if (interactive()) save(fev1.msm, file="~/msm/devel/models/fev1.msm.rda")
    if (developer.local) system("../../hid/hidden ../../hid/newcmdfiles/onelevel3covfit.cmd") # OK - see fev2.out


    hmodel2 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                    hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                    hmmIdent(999))
    (fev2.msm <- msm(fev ~ months, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel2,
                     hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                     hconstraint = list(sderr = c(1,1), acute = c(1,1)), 
                     method="BFGS", control=list(trace=3, REPORT=1), 
                     fixedpars=FALSE, center=FALSE))
    print(fev2.msm)

    (fev2.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel3,
                     hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                     hconstraint = list(sderr = c(1,1), acute = c(1,1)), 
                     control=list(trace=3, REPORT=1), 
                     method="BFGS", fixedpars=FALSE))
    print(fev2.msm)

    if (interactive()) save(fev1.msm, file="~/msm/devel/models/fev1.msm.rda")
    if (developer.local) system("../../hid/hidden ../../hid/newcmdfiles/twolevel3covfit.cmd") # OK - see fev1.out


#########  OUTPUT FUNCTIONS FOR HIDDEN MARKOV MODELS  ############################

    if (interactive()) load(file="~/msm/devel/models/fev1.msm.rda")
    if (interactive()) load(file="~/msm/devel/models/fev2.msm.rda")

    print(viterbi.msm(fev1.msm)[1:100,])
    print(qmatrix.msm(fev1.msm))
    print(sojourn.msm(fev1.msm))
    print(pmatrix.msm(fev1.msm, 10))
    print(qratio.msm(fev1.msm, c(1,2), c(2,3), cl=0.99))
    print(prevalence.msm(fev1.msm))
    print(summary.msm(fev1.msm))
    if (interactive()) plot.msm(fev1.msm)
    print(coef.msm(fev1.msm))
    print(hazard.msm(fev1.msm))
    print(transient.msm(fev1.msm))
    print(absorbing.msm(fev1.msm))
    print(totlos.msm(fev1.msm))
    print(logLik.msm(fev1.msm))


### example of Viterbi algorithm

    keep <- fev$ptnum==1 & fev$fev<999
    vit <- viterbi.msm(fev1.msm)[keep,]
    print(max1 <- max(vit$time[vit$fitted==1]))
    print(min2 <- min(vit$time[vit$fitted==2]))
    if (interactive())  {
        plot(fev$days[keep], fev$fev[keep], type="l", ylab=expression(paste("% baseline ", FEV[1])), xlab="Days after transplant")
        abline(v = mean(max1,min2), lty=2)
        text(max1 - 500, 50, "STATE 1")
        text(min2 + 500, 50, "STATE 2")
    }
}
