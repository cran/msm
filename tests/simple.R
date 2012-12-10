source("local.R")
library(msm)
#library(msm, lib.loc="~/lib/R")
#library(msm, lib.loc="~/msm/lib/sun/0.7.5")
data(cav)

### TESTS FOR SIMPLE NON-HIDDEN MARKOV MODELS

twoway4.q <- rbind(c(-0.5, 0.25, 0, 0.25), c(0.166, -0.498, 0.166, 0.166), c(0, 0.25, -0.5, 0.25), c(0, 0, 0, 0))
twoway4.q2 <- rbind(c(-0.51, 0.24, 0, 0.25), c(0.162, -0.498, 0.168, 0.166), c(0, 0.26, -0.5, 0.25), c(0, 0, 0, 0))
twoway3.q <- rbind(c(-0.5, 0.25, 0), c(0.166, -0.498, 0.166), c(0, 0.25, -0.5))
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171), c(0, 0, 0.202, 0.081), c(0, 0, 0, 0.126), c(0, 0, 0, 0))
rownames(twoway4.q) <- colnames(twoway4.q) <- c("Well","Mild","Severe","Death")
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well","Mild","Severe","Death")
twoway4.i <- twoway4.q; twoway4.i[twoway4.i!=0] <- 1
oneway4.i <- oneway4.q; oneway4.i[oneway4.i!=0] <- 1

### CAV DATA
(stab <- statetable.msm(state, PTNUM, data=cav))
stopifnot(all(stab == c(1367, 46, 4, 204, 134, 13, 44, 54, 107, 148, 48, 55)))
(cinits <- crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q))
stopifnot(isTRUE(all.equal(c(-0.117314905786477, 0.116817878212849, 0, 0, 0.067989320398981, -0.375848825554382, 0.049084006577444, 0, 0, 0.137134030945518, -0.256747111328168, 0, 0.049325585387496, 0.121896916396016, 0.207663104750724, 0), as.numeric(cinits))))

## Simple model
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                 qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1))
stopifnot(isTRUE(all.equal(4908.81676837903, cav.msm$minus2loglik, tol=1e-06)))
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                 qmatrix = cinits, death = TRUE, fixedpars=TRUE,
                 method="BFGS", control=list(trace=5, REPORT=1))
stopifnot(isTRUE(all.equal(4113.16601901957, cav.msm$minus2loglik, tol=1e-06)))

if (developer.local) {
    system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                               qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, # opt.method="nlm"
                               method="BFGS",
                               control=list(trace=5, REPORT=1, fnscale=1)
                               )
                )
    stopifnot(isTRUE(all.equal(3968.7978930519, cav.msm$minus2loglik, tol=1e-06)))
    
    pnext.msm(cav.msm)
    pnext.msm(cav.msm, ci="normal")
    pnext.msm(cav.msm, ci="bootstrap", B=3)
    if (0)
        pnext.msm(cav.msm, ci="bootstrap")   
}

## No death state.
system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                 qmatrix = twoway4.q, death = FALSE, fixedpars=FALSE,
                 method="BFGS", control=list(trace=1, REPORT=1, fnscale=4000)) )
stopifnot(isTRUE(all.equal(3986.08765893935, cav.msm$minus2loglik, tol=1e-06)))
cav.msm

## auto-generated initial values.
state.g <- cav$state; time.g <- cav$years; subj.g <- cav$PTNUM
cav.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = twoway4.i, gen.inits=TRUE, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4119.9736299032, cav.msm$minus2loglik, tol=1e-06)))

if (developer.local) {
    system.time(cav.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = twoway4.i, gen.inits=TRUE, fixedpars=TRUE))
    stopifnot(isTRUE(all.equal(4119.9736299032, cav.msm$minus2loglik, tol=1e-06)))
    cav.msm <- msm(state.g ~ time.g, subject=subj.g, qmatrix = crudeinits.msm(state ~ years, PTNUM, twoway4.i, data=cav), fixedpars=TRUE)
    stopifnot(isTRUE(all.equal(4119.9736299032, cav.msm$minus2loglik, tol=1e-06)))
}

## Covariates

cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                 qmatrix = twoway4.q, death = TRUE, fixedpars=TRUE,
                 covariates = ~ sex, covinits = list(sex=rep(0.01, 7)), # , dage=rep(0, 7)),
                 method="BFGS", control=list(trace=5, REPORT=1))
stopifnot(isTRUE(all.equal(4909.08442586298, cav.msm$minus2loglik, tol=1e-06)))

if (developer.local) {
    system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                                 qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,
                                 covariates = ~ sex, method="BFGS", control=list(trace=5, REPORT=1))) # 44.13 on new, 260.14 on old
    stopifnot(isTRUE(all.equal(3954.77699876128, cav.msm$minus2loglik, tol=1e-06)))
    qmat <- qmatrix.msm(cav.msm)[c("estimates","SE")]
    stopifnot(isTRUE(all.equal(c(0.226768225214069, -0.583572966218831, 0.337068668450776, 0.0197360725539862), as.numeric(qmat$estimates[2,]), tol=1e-06)))
    stopifnot(isTRUE(all.equal(c(0.0341912476474469, 0.177757314782152, 0.0383135782775986, 0.171176619065207), as.numeric(qmat$SE[2,]), tol=1e-02)))   ### SEs slightly different on different machines.
    stopifnot(isTRUE(all.equal(5.29003443121721, sojourn.msm(cav.msm)[1,3], tol=1e-06)))
}

if (developer.local) {
    ## Baseline constraints
    system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                                 qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,
                                 qconstraint = c(1,1,2,2,2,3,3),
                                 method="BFGS", control=list(trace=2, REPORT=1)
                                 )) # 3.22 on new, 17.55 on old
    stopifnot(isTRUE(all.equal(4116.22686367935, cav.msm$minus2loglik, tol=1e-06)))

    qmat <- qmatrix.msm(cav.msm)[c("estimates","SE")]
    stopifnot(isTRUE(all.equal(0.171691686850913, qmat$estimates[2,1], tol=1e-06)))
    stopifnot(isTRUE(all.equal(2.96494230954565, sojourn.msm(cav.msm)[3,1], tol=1e-06)))

    ## Covariate constraints.
    system.time(cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE,
                                 covariates = ~ sex, covinits = list(sex=rep(0.01, 7)), constraint=list(sex=c(1,2,3,1,2,3,2)),
                                 method="BFGS", control=list(trace=1, REPORT=1))) # 21.2 on new, 86.58 on old
    stopifnot(isTRUE(all.equal(3959.35551766943, cav.msm$minus2loglik, tol=1e-06)))
    qmat <- qmatrix.msm(cav.msm)[c("estimates","SE")]
    stopifnot(isTRUE(all.equal(0.228458428847816, qmat$estimates[2,1], tol=1e-06)))
}

## Constraints with psoriatic arthritis data
data(psor)
psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
system.time(psor.msm <- msm(state ~ months, subject=ptnum, data=psor,
                qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, # covinits=list(hieffusn = c(0.5, 0.1, 0), ollwsdrt=c(0.2, 0.1, -0.1)),
                 constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), 
                fixedpars=FALSE, control = list(REPORT=1,trace=2), method="BFGS"))
stopifnot(isTRUE(all.equal(1114.89946121717, psor.msm$minus2loglik, tol=1e-06)))
stopifnot(isTRUE(all.equal(0.0953882330391683, qmatrix.msm(psor.msm)$estimates[1,2], tol=1e-03)))

pn <- pnext.msm(psor.msm)

## center=FALSE
psor.nocen.msm <- msm(state ~ months, subject=ptnum, data=psor,
                    qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn,
                    constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), 
                    fixedpars=FALSE, center=FALSE, control = list(REPORT=1,trace=2), method="BFGS")
## $baseline and component should obey center
stopifnot(isTRUE(all.equal(psor.nocen.msm$Qmatrices$baseline, qmatrix.msm(psor.nocen.msm, covariates=0, ci="none"))))
stopifnot(isTRUE(all.equal(psor.nocen.msm$Qmatrices$baseline, qmatrix.msm(psor.nocen.msm, covariates=list(ollwsdrt=0), ci="none"))))
stopifnot(isTRUE(all.equal(psor.msm$Qmatrices$baseline, qmatrix.msm(psor.msm, covariates="mean", ci="none"))))
stopifnot(isTRUE(all.equal(psor.nocen.msm$sojourn, sojourn.msm(psor.nocen.msm, covariates=0))))
stopifnot(isTRUE(all.equal(psor.msm$sojourn, sojourn.msm(psor.msm, covariates="mean"))))
stopifnot(isTRUE(all.equal(exp(psor.nocen.msm$Qmatrices$logbaseline[c(5,10,15)]), psor.nocen.msm$Qmatrices$baseline[c(5,10,15)])))
stopifnot(isTRUE(all.equal(exp(psor.msm$Qmatrices$logbaseline[c(5,10,15)]), psor.msm$Qmatrices$baseline[c(5,10,15)])))

## Bug fixed in 1.1.3
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.nocen.msm, covariates=0)$SE[1,2],
                           qmatrix.msm(psor.nocen.msm, covariates=list(hieffusn=0, ollwsdrt=0))$SE[1,2])))
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.nocen.msm, covariates=0)$SE[1,1],
                           qmatrix.msm(psor.nocen.msm, covariates=list(hieffusn=0, ollwsdrt=0))$SE[1,1])))
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.nocen.msm, covariates=0)$L[1,2],
                           qmatrix.msm(psor.nocen.msm, covariates=list(hieffusn=0, ollwsdrt=0))$L[1,2])))
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.nocen.msm, covariates=0)$L[1,2],
                           qmatrix.msm(psor.nocen.msm, covariates=list(hieffusn=0))$L[1,2])))
cm <- psor.nocen.msm$qcmodel$covmeans
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.nocen.msm, covariates="mean")$SE[1,2],
                           qmatrix.msm(psor.nocen.msm, covariates=list(hieffusn=cm["hieffusn"], ollwsdrt=cm["ollwsdrt"]))$SE[1,2])))
                           
## values not supplied are set to zero, even if center=TRUE (undocumented)
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.msm, covariates=list(hieffusn=0)),
                           qmatrix.msm(psor.msm, covariates=list(ollwsdrt=0, hieffusn=0)))))
stopifnot(isTRUE(all.equal(qmatrix.msm(psor.nocen.msm, covariates=list(hieffusn=0)),
                           qmatrix.msm(psor.nocen.msm, covariates=list(ollwsdrt=0, hieffusn=0)))))


## No death state
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                 qmatrix = twoway4.q, death = FALSE, fixedpars=TRUE)
stopifnot(isTRUE(all.equal(4833.0064065267, cav.msm$minus2loglik, tol=1e-06)))

## how big a dataset can 
if (0) { 
  c3.df <- NULL
  for (i in 1:10) {c22.df <- c2.df; c22.df$PTNUM <- c2.df$PTNUM + 120000*(i-1); c3.df <- rbind(c3.df, c22.df)}
  length(unique(c2.df$PTNUM))
  length(unique(c3.df$PTNUM))
  qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
  c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                qmatrix=qx, death=c(4, 5), method="BFGS", 
                control=list(trace=1, REPORT=1, fnscale=100000))
  c3.msm <- msm(state~years, subject=PTNUM, data=c3.df,
                qmatrix=qx, death=c(4, 5), method="BFGS", 
                control=list(trace=1, REPORT=1, fnscale=100000))
}

## Multiple death states (Jean-Luc's data)
if (developer.local) {
    c2.df <- read.table("~/msm/tests/jeanluc/donneesaveccancerPT.txt", header=TRUE)
    qx <- rbind( c(0, 0.005, 0, 0, 0), c(0, 0, 0.01, 0.02,0), c(0, 0, 0, 0.04, 0.03), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                  qmatrix=qx, death=c(4, 5), method="BFGS", fixedpars = 1:5,
                  control=list(trace=2, REPORT=1, fnscale=100000))
    stopifnot(isTRUE(all.equal(70646.727505836, c2.msm$minus2loglik, tol=1e-06)))
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df,
                  qmatrix=qx, method="BFGS", fixedpars = 1:5,
                  control=list(trace=2, REPORT=1, fnscale=100000))
    stopifnot(isTRUE(all.equal(62915.1638036017, c2.msm$minus2loglik, tol=1e-06)))

    ## Same using an "obstype" vector.
    obstype <- ifelse(c2.df$state %in% c(4,5), 3, 1)
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, qmatrix=qx,
                  obstype=obstype, method="BFGS", fixedpars = 1:5,
                  control=list(trace=2, REPORT=1, fnscale=100000))
    stopifnot(isTRUE(all.equal(70646.727505836, c2.msm$minus2loglik, tol=1e-06)))
    obstype <- rep(1, length(c2.df$state))
    c2.msm <- msm(state~years, subject=PTNUM, data=c2.df, obstype=obstype,
                  qmatrix=qx, method="BFGS", fixedpars = 1:5,
                  control=list(trace=2, REPORT=1, fnscale=100000))
    stopifnot(isTRUE(all.equal(62915.1638036017, c2.msm$minus2loglik, tol=1e-06)))

### G Marshall's diabetic retinopathy data
    marsh.df <- read.table("~/msm/tests/markov/test.dat", col.names=c("subject","eyes","time","duration","hba1"))
    marsh.df$hba1 <- marsh.df$hba1 - mean(marsh.df$hba1)
    marsh.msm <-
      msm(eyes ~ time, subject=subject, qmatrix = rbind(c(0,0.02039,0,0), c(0.007874,0,0.01012,0), c(0,0.01393,0,0.01045), c(0,0,0,0)),
          covariates = ~ hba1, data = marsh.df, fixedpars=TRUE)
    stopifnot(isTRUE(all.equal(335.897217906310, marsh.msm$minus2loglik, tol=1e-06)))
    system.time(marsh.msm <-
      msm(eyes ~ time, subject=subject, qmatrix = rbind(c(0,0.02039,0,0), c(0.007874,0,0.01012,0), c(0,0.01393,0,0.01045), c(0,0,0,0)),
          covariates = ~ hba1, data = marsh.df))
    stopifnot(isTRUE(all.equal(310.989863258621, marsh.msm$minus2loglik, tol=1e-06)))
    stopifnot(isTRUE(all.equal(-0.0496235211442196, qmatrix.msm(marsh.msm, covariates=0)$estimates[1,1], tol=1e-06)))
}


### Aneurysm dataset different in 0.5 (not fromto, includes imputed initial states)


##########    OUTPUT FUNCTIONS    ###################

qmatrix.msm(psor.msm)
stopifnot(isTRUE(all.equal(c(-0.0953882330391683, 0, 0, 0, 0.0953882330391683, -0.163370011525553, 0, 0, 0, 0.163370011525553, -0.255229343798597, 0, 0, 0, 0.255229343798597, 0), as.numeric(qmatrix.msm(psor.msm)$estimates), tol=1e-03)))
stopifnot(isTRUE(all.equal(c(0.0115507014188511, 0, 0, 0, 0.0115507014188511, 0.0195265275850904, 0, 0, 0, 0.0195265275850904, 0.0378507662232158, 0, 0, 0, 0.0378507662232158, 0), as.numeric(qmatrix.msm(psor.msm)$SE), tol=1e-03)))
qmat <- qmatrix.msm(psor.msm, covariates=list(ollwsdrt=0.1, hieffusn=0.4))
stopifnot(isTRUE(all.equal(c(-0.121430585652200, 0, 0, 0, 0.121430585652200, -0.207972362475868, 0, 0, 0, 0.207972362475868, -0.257535341208494, 0, 0, 0, 0.257535341208494, 0), as.numeric(qmat$estimates), tol=1e-03)))
stopifnot(isTRUE(all.equal(c(0.0162156605802465, 0, 0, 0, 0.0162156605802465, 0.0266727053124233, 0, 0, 0, 0.0266727053124233, 0.0364321127089265, 0, 0, 0, 0.0364321127089265, 0), as.numeric(qmat$SE), tol=1e-04)))
try(qmatrix.msm(psor.msm, covariates=list(hieffusn=0.1, foo=0.4))) # deliberate error
qmat <- qmatrix.msm(psor.msm, covariates=list(ollwsdrt=0.1, hieffusn=0.4), cl=0.99)
stopifnot(isTRUE(all.equal(c(-0.171282667596986, 0, 0, 0, 0.0860880282792585, -0.289385121267802, 0, 0, 0, 0.149463467106753, -0.370756460718086, 0, 0, 0, 0.178889538008097, 0), as.numeric(qmat$L), tol=1e-04)))
soj <- qmatrix.msm(psor.msm, covariates=list(ollwsdrt=0.1, hieffusn=0.4), sojourn=TRUE)$sojourn
stopifnot(isTRUE(all.equal(c(8.23515751512713, 4.80833120370037, 3.88296221911705, Inf), as.numeric(soj), tol=1e-03)))
qmatrix.msm(psor.msm, ci="normal", B=2)
qmatrix.msm(psor.msm, ci="boot", B=2)

soj <- sojourn.msm(psor.msm, covariates=list(ollwsdrt=0.1, hieffusn=0.4))
stopifnot(isTRUE(all.equal(c(8.23515751512713, 4.80833120370037, 3.88296221911705, 1.09971073904434, 0.616674252838334, 0.549301375677405, 6.33875136203292, 3.73961380505919, 2.94271599303942, 10.6989240349703, 6.18246967994404, 5.12363260020806), as.numeric(unlist(soj)), tol=1e-04)))
soj <- sojourn.msm(psor.msm, covariates=list(ollwsdrt=0.1, hieffusn=0.4), cl=0.99)
stopifnot(isTRUE(all.equal(5.83830234564607, soj[1,"L"], tol=1e-04)))

stopifnot(isTRUE(all.equal(0.148036812842411, pmatrix.msm(psor.msm, ci="none", t=10)[1,3], tol=1e-04)))
try(pmatrix.msm(psor.msm, t=10, covariates=list(hieffusn=0.1))) # deliberate error
p <- pmatrix.msm(psor.msm, t=10, covariates=list(ollwsdrt=0.1, hieffusn=0.2))
stopifnot(isTRUE(all.equal(0.18196160265907, p[1,3], tol=1e-04)))

set.seed(22061976); stopifnot(isTRUE(all.equal(0.12, pmatrix.msm(psor.msm, ci="normal", B=3)$L[2,3], tol=1e-01, scale=1)))

q <- qratio.msm(psor.msm, c(1,2), c(2,3))
stopifnot(isTRUE(all.equal(c(0.583878474075081, 0.0996029045389022, 0.417943274168735, 0.815694601537263), as.numeric(q), tol=1e-04)))
q <- qratio.msm(psor.msm, c(1,2), c(2,3), cl=0.99)
stopifnot(isTRUE(all.equal(0.376262194364283, as.numeric(q["L"]), tol=1e-04)))
q <- qratio.msm(psor.msm, c(1,1), c(2,3))
stopifnot(isTRUE(all.equal(c(-0.583878474075081, 0.0996029045389022, -0.815694601537263, -0.417943274168735), as.numeric(q), tol=1e-04)))
q <- qratio.msm(psor.msm, c(2,2), c(2,3))
stopifnot(isTRUE(all.equal(c(-1,0,-1,-1), as.numeric(q), tol=1e-06)))

qratio.msm(psor.msm, c(1,2), c(2,3), ci="norm", B=2)
qratio.msm(psor.msm, c(1,2), c(2,3), ci="boot", B=2)

print(interactive())
if (interactive())
  {
      plot.msm(psor.msm)
      plot.msm(psor.msm, from=c(1,3), to=4, range=c(10,30))
      plot.msm(psor.msm, from=c(1,2), to=4, range=c(10,80), legend.pos=c(70,0.1))

      surface.msm(psor.msm)
      surface.msm(psor.msm, type="filled")
      surface.msm(psor.msm, c(3,5), type="filled")
      surface.msm(psor.msm, c(3,4))
      x <- psor.msm$paramdata$params.uniq
      x[6] <- 0
      surface.msm(psor.msm, c(3,4), point=x)
      surface.msm(psor.msm, c(3,4), point=x, xrange=c(-2, -0.6))
      surface.msm(psor.msm, c(3,4), point=x, yrange=c(-1.2, 0.1))
      surface.msm(psor.msm, c(3,4), point=x, np = 5)
      contour(psor.msm)
      persp(psor.msm)
      persp(psor.msm, np=5)
      image(psor.msm)
      plot.prevalence.msm(psor.msm)
      plot.prevalence.msm(psor.msm, col.obs="green", col.exp="brown", lty.obs=2, lty.exp=3, lwd.obs=2, lwd.exp=4, cex.axis=2)

      ##
      plotprog.msm(state ~ months, subject=ptnum, data=psor, legend.pos=c(20,0.99), lwd=3, xlab="Months")
      plotprog.msm(state ~ months, subject=ptnum, data=psor, legend.pos=c(20,0.99), lwd=1, mark.time=FALSE, xlab="Months")
      plot.survfit.msm(psor.msm, lwd=3, xlab="Months")
      plot.survfit.msm(psor.msm, lwd=3, lwd.surv=3, mark.time=FALSE, col.surv="green", lty.surv=4, xlab="Months")
      plot.survfit.msm(psor.msm, lwd=3, lwd.surv=3, ci="normal", B=4, col.surv="green", lty.surv=4, xlab="Months",
                       lty.ci=1, lwd.ci=2, col.ci="purple")
    }


co <- coef.msm(psor.msm)
stopifnot(isTRUE(all.equal(0.498319866154661, co$hieffusn[1,2], tol=1e-04)))

haz <- hazard.msm(psor.msm)
stopifnot(isTRUE(all.equal(0.385347226135311, haz$ollwsdrt[1,2], tol=1e-04)))
stopifnot(isTRUE(all.equal(0.385347226135311, haz$ollwsdrt[2,2], tol=1e-04)))
stopifnot(isTRUE(all.equal(2.35928404626333, haz$hieffusn[1,3], tol=1e-04)))
stopifnot(isTRUE(all.equal(2.35928404626333, haz$hieffusn[3,3], tol=1e-04)))

haz <- hazard.msm(psor.msm, hazard.scale=2)
stopifnot(isTRUE(all.equal(0.148492484690178, haz$ollwsdrt[1,2], tol=1e-04)))
stopifnot(isTRUE(all.equal(0.148492484690178, haz$ollwsdrt[2,2], tol=1e-04)))
#stopifnot(isTRUE(all.equal(5.56622121095267, haz$hieffusn[1,3], tol=1e-04)))
#stopifnot(isTRUE(all.equal(5.56622121095267, haz$hieffusn[3,3], tol=1e-04)))

haz <- hazard.msm(psor.msm, hazard.scale=c(1,2))
stopifnot(isTRUE(all.equal(0.385347226135311, haz$ollwsdrt[1,2], tol=1e-04)))
stopifnot(isTRUE(all.equal(0.385347226135311, haz$ollwsdrt[2,2], tol=1e-04)))
#stopifnot(isTRUE(all.equal(5.56622121095267, haz$hieffusn[1,3], tol=1e-04)))
#stopifnot(isTRUE(all.equal(5.56622121095267, haz$hieffusn[3,3], tol=1e-04)))

stopifnot(isTRUE(all.equal(c(1,2,3), as.numeric(transient.msm(psor.msm)), tol=1e-06)))

stopifnot(isTRUE(all.equal(4, as.numeric(absorbing.msm(psor.msm)), tol=1e-06)))




cat("simple.R: ALL TESTS PASSED\n")

