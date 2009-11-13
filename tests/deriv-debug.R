### Derivatives at MLE for CAV model

## derivsimple sums over fromstate/tostate/timelag/obstype combinations. 1190 of these in cav.msm
if (0){ 
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = twoway4.q, death = TRUE, fixedpars=FALSE, method="BFGS" )
q.mle <- cav.msm$estimates
cav.msm$minus2loglik
likderiv.msm(q.mle, deriv=0, cav.msm$data, cav.msm$qmodel, cav.msm$qcmodel, cav.msm$cmodel, cav.msm$hmodel, cav.msm$paramdata)
likderiv.msm(q.mle, deriv=1, cav.msm$data, cav.msm$qmodel, cav.msm$qcmodel, cav.msm$cmodel, cav.msm$hmodel, cav.msm$paramdata)
colSums(likderiv.msm(q.mle, deriv=3, cav.msm$data, cav.msm$qmodel, cav.msm$qcmodel, cav.msm$cmodel, cav.msm$hmodel, cav.msm$paramdata)) # should be about 0 
}
