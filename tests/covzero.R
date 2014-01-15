### COVARIATES ON TRANSITION INTENSITIES

if (0) {

data(psor)
psor <- read.table("~/work/msm/msm/data/psor.txt", header=TRUE)
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))

msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, fixedpars=c(4,5,7))
msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn))

### additional fixedpars
msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, fixedpars=c(4,5,7))
msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), fixedpars=6)
try(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), fixedpars=7)) # fixedpars out of range in new, but in range in old

### constraints
ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint=list(hieffusn=c(1,1)))
try(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint=list(hieffusn=c(1,1,1)))) # wrong length
try(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint="foo")) # junk in constraints
try(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), constraint=list(hieffusn=c("foo",1)))) # junk in constraints

### alternative covinits
msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), covinits=list(hieffusn=c(0.1,0.1)), fixedpars=TRUE)
msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, covinits=list(hieffusn=c(0,0.1,0.1)), fixedpars=TRUE)
try(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), covinits=list(hieffusn=c(0, 0.1,0.1)), fixedpars=TRUE)) # covinits wrong length
try(msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn), covinits=list(hieffusn=c("foo")), fixedpars=TRUE)) # junk in covinits

### factors
psor$faccov <- factor(sample(1:3, size=nrow(psor), replace=TRUE))
ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn+faccov))
qmatrix.msm(ps, covariates=0)
qmatrix.msm(ps, covariates=list(faccov=1))
qmatrix.msm(ps, covariates=list(faccov=2))
### constraints on factors
try(ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn+faccov, "2-3" = ~hieffusn+faccov), constraint=list(faccov=c(1,1))))
ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn+faccov, "2-3" = ~hieffusn+faccov), constraint=list(faccov2=c(1,1)))

### ordered factors
psor$faccov <- ordered(sample(1:3, size=nrow(psor), replace=TRUE))
ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt+hieffusn, "2-3" = ~hieffusn+faccov))
qmatrix.msm(ps, covariates=0)
qmatrix.msm(ps, covariates=list(faccov=1))
qmatrix.msm(ps, covariates=list(faccov=2))
## added to package todo, warn that poly contrasts not supported in extractors

### interactions
ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = list("3-4"=~ollwsdrt*hieffusn, "2-3" = ~hieffusn*ollwsdrt))
ps


### with a hmm as well
emat <- rbind(c(0,0.1,0,0),c(0.1,0,0,0),c(0,0,0,0),c(0,0,0,0))
ps <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, ematrix=emat, covariates = list("3-4"=~ollwsdrt*hieffusn, "2-3" = ~hieffusn), initprobs=c(0.5,0.2,0.2,0.1), fixedpars=c(8,9))
ps



### COVARIATES ON MISCLASSIFICATION PROBS
### don't support. see below.
if (0) {
cav <- read.table("~/work/msm/msm/data/cav.txt", header=TRUE)
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))

misccov.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                   qmatrix = oneway4.q, ematrix=ematrix, death = 4, fixedpars = TRUE,
                   control = list(trace=1, REPORT=1), method="BFGS",
                   misccovariates = list("2-1" = ~ sex), misccovinits=list(sex=rep(0.05, 1)))
}

## ecmodel has already been converted into hmodel at this point!
## note     hcovariates <- lapply(ifelse(rowSums(emodel$imatrix)>0, deparse(misccovariates), deparse(~1)), as.formula)
## makes a list of formulae, one for each state
## these are then processed by msm.form.covdata(hcovariates....)
## MN log reg, cov effs are on log(p/pbase). same for all to within from state
## hcovariates has one element for each state.
## note this comment in miscnew.R: "Don't allow, since hcovariates doesn't logically correspond to ematrix."
## just say sod it? don't need it for current work anyway.

### bits removed from msm.R when decided to drop this
if (0) {
        ## if (is.list(misccovariates))
        ##     misccovdata <- msm.form.covdata.byrate(misccovariates, emodel, data, NULL, center)
        ## else if (inherits(misccovariates, "formula"))
        ## else stop(deparse(substitute(misccovariates)), " should be a formula or list of formulae")

        ## ### TESTME same for misc probs, and both.
        ## offset <- qmodel$ndpars + qcmodel$ndpars + emodel$ndpars
        ## inds[offset + ecmodel$constr[!duplicated(ecmodel$constr)]] <-
        ##     ecmodel$cri[!duplicated(ecmodel$constr)]


}

}
