source("local.R")
library(msm)

if (developer.local) {
  data(psor)
  psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")

  set.seed(22061976)
  q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline, file="~/msm/devel/psor.q.boot.rda")
  ## new SEs, if more accurate, should be greater than asymp SE from Hessian.  OK.
  apply(array(unlist(q.list), dim=c(4,4,5)), c(1,2), sd)
  apply(array(unlist(q.list), dim=c(4,4,5)), c(1,2), function(x)quantile(x, c(0.025, 0.975)))
  ## old SEs
  psor.msm$QmatricesSE$baseline

  ## Pmatrix CI utility function
  pmatrix.msm(psor.msm)
  (p.ci <- pmatrix.msm(psor.msm, ci.boot=TRUE, B=5))

  t.list <- boot.msm(psor.msm, totlos.msm, B=5)
  t.ci <- totlos.ci.msm(psor.msm, B=3)
  totlos.msm(psor.msm, ci.boot=TRUE, B=3)

  ## test on a HMM
  load(file="~/msm/devel/models/misc.msm.rda")
  pmatrix.msm(misc.msm, 1)
  oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
  rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
  ematrix <- rbind(c(0, 0.1, 0, 0),c(0.1, 0, 0.1, 0),c(0, 0.1, 0, 0),c(0, 0, 0, 0))
  p.list <- boot.msm(misc.msm, function(x)pmatrix.msm(x, t=1), 3)

  ### models with censoring
  cav.cens <- cav
  cav.cens$state[cav$state==4][1:50] <- 99
  twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
  cavcens.msm <- msm(state ~ years, subject=PTNUM, data=cav.cens, qmatrix=twoway4.q, censor=99, fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
  p.list <- boot.msm(cavcens.msm, B=3)


  ## parallel processing
  ## Using doSMP now since multicore doesn't work on Windows Vista or 7
  ## ?? Passing environment works from package, but not with manual function loading

  data(psor)
  psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix =
                  psor.q, covariates = ~ollwsdrt+hieffusn, constraint =
                  list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), fixedpars=FALSE)

system.time(q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline, B=50, cores=3))
system.time(q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline, B=50, cores=1))

  pmatrix.msm(psor.msm, B=50, ci="boot", cores=3)
  pmatrix.msm(psor.msm, B=50, ci="normal")
  qmatrix.msm(psor.msm, B=50, ci="boot", cores=3)
  qmatrix.msm(psor.msm, B=50, ci="normal")
  pnext.msm(psor.msm, B=50, ci="boot", cores=3)
  pnext.msm(psor.msm, B=50, ci="normal")
  qratio.msm(psor.msm, c(1,2), c(2,3), B=50, ci="boot", cores=3)
  qratio.msm(psor.msm, c(1,2), c(2,3), B=50, ci="normal")
  totlos.msm(psor.msm, B=50, ci="boot", cores=3)
  totlos.msm(psor.msm, B=50, ci="normal")
  expected.msm(psor.msm, B=50, ci="boot", cores=3)
  expected.msm(psor.msm, B=50, ci="normal")
  pmatrix.piecewise.msm(psor.msm, B=50, ci="boot", cores=3)
  pmatrix.piecewise.msm(psor.msm, B=50, ci="normal")

  times <- c(5, 10, 15)
  covariates <- list(list(ollwsdrt=0, hieffusn=0),
                     list(ollwsdrt=0, hieffusn=1),
                     list(ollwsdrt=1, hieffusn=0),
                     list(ollwsdrt=1, hieffusn=1)
                     )
  pmatrix.piecewise.msm(psor.msm, 0, 3, times, covariates, ci="boot", cores=3, B=50)
  pmatrix.piecewise.msm(psor.msm, 0, 3, times, covariates, ci="normal", B=50)

  if (0){
      times <- numeric(10)
      for (i in 1:10)
          times[i] <- system.time(q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline, file="~/msm/devel/psor.q.boot.rda", B=100, cores=i))["elapsed"]
      print(times)
### on 4 core machine
### [1] 80.117 41.936 37.123 26.334 27.725 22.678 24.545 25.094 23.127 25.674
  }

}
