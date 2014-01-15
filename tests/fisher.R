### TODO more testing

source("local.R")

if(developer.local){

    library(msm)
    data(cav)
    twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
    data(psor)
    psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))

#source("../R/msm.R")
#if (is.loaded("msmCEntry")) dyn.unload(lib);
#dyn.load(lib)

## Fisher scoring is over twice as fast as previous best method (BFGS)!!

system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="fisher", qmatrix = twoway4.q, fixedpars=FALSE, control=list(trace=1)))
# stopifnot(isTRUE(all.equal(3968.7978930519, cav.msm$minus2loglik, tol=1e-06)))
cav.msm$minus2loglik
cav.msm$paramdata$estimates
cav.msm$paramdata$ci
#cav.msm$paramdata$deriv
#cav.msm$paramdata$info
#cav.msm$paramdata$opt$hessian
#sqrt(diag(solve(0.5*cav.msm$paramdata$info)))
#sqrt(diag(solve(0.5*cav.msm$paramdata$opt$hessian)))

system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="optim", method="BFGS", qmatrix = twoway4.q, fixedpars=FALSE))
cav.msm$minus2loglik
cav.msm$paramdata$estimates
cav.msm$paramdata$ci

system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="optim", method="Nelder-Mead", qmatrix = twoway4.q, fixedpars=FALSE))
cav.msm$minus2loglik # converges to local max
cav.msm$paramdata$estimates
cav.msm$paramdata$ci

system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="nlm", qmatrix = twoway4.q, fixedpars=FALSE))
cav.msm$paramdata$estimates
cav.msm$paramdata$ci

## expected versus observed information
cav.msm$paramdata$info
cav.msm$opt$hessian

solve(0.5*cav.msm$paramdata$info)
cav.msm$covmat # should be nearly the same


system.time(psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,
                            covariates = ~ollwsdrt+hieffusn,  method="BFGS",
                            constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2))))
psor.msm$paramdata$estimates
psor.msm$paramdata$ci
system.time(psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q,
                            covariates = ~ollwsdrt+hieffusn,  opt.method="fisher",
                            constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2))))
psor.msm$paramdata$estimates
psor.msm$paramdata$ci
## OK.


## optim improves lik but falls over with singular info mat close to the MLE.
## unless damp parameter provided, as in Andrew's code.
system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, fixedpars=FALSE, covariates=~sex, method="BFGS", control=list(trace=1)))
qmatrix.msm(cav.msm)
      cav.msm$minus2loglik

    system.time(cavf.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="fisher", qmatrix = qmatrix.msm(cav.msm,ci="none"), fixedpars=FALSE, covariates=~sex, covinits=list(sex=cav.msm$Qmatrices$sex[cbind(c(1,1,2,2,2,3,3),c(2,4,1,3,4,2,4))]), control=list(trace=1,damp=1)))

    system.time(cavf.msm <- msm( state ~ years, subject=PTNUM, data = cav, opt.method="fisher", qmatrix = twoway4.q, fixedpars=FALSE, covariates=~sex, control=list(trace=1,damp=0.5)))

cavf.msm$minus2loglik




}
