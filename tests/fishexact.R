## 30 March 2013: Can't be bothered to do fisher info for exact times.
##  or exact death times.  Does Kay do that?
## Andrew's code doesn't do exact times, so can't test against his
## Theory is different from panel data, because observation times are part of the outcome, not the design.

##  Formula for models with just q, no constraints or covs:

## Info contrib for an observed:
## 1-2 trans
## L = log(exp(-(q12+q13)dt) * q12) = -(q12 + q13) dt + log(q12) = log(pij)
## dL/dq12 = -dt + 1/q12,  dL/dq13 = -dt
## d2L/dq12 dq12 = -1/q12^2,  others = 0
## 1-3 trans
## L = = -(q12 + q13) dt + log(q13)
## dL/dq12 = -dt,  dL/dq13 = -dt + 1/q13
## d2L/dq13^2 = -1/q13^2,  others = 0

# Exact info diag
#1 / (0.01 + 0.004) * c(1/0.01, 1/0.004)
#1 / q[1:2] /(q[1]+q[2])
#1/q[1]^2 * 0.01 / 0.014   # * P(state=2)
#1/q[2]^2 * 0.004 / 0.014  # * P(state=3)

## So:
## For all i-j trans, info diag is 1 / (vec qi) / sum_j(q_ij),   off diags 0
## For i-i transitions, obs of i was by design, and guaranteed to be same as last, so no info contrib.
## TODO: Constraints on q
## For obs trans to j (k is sum over others)
## L = = -(q1k + q1j) dt + log(q1j)

## Constraints within rows

## Constraints between rows

## TODO: Covariates, and constraints on b

## Far too much faff.





## CODE TO USE FOR TESTING IF EVER USE THIS:
if (0) {
fiveq <- rbind(c(0,0.01,0,0,0.002), c(0,0,0.07,0,0.01), c(0,0,0,0.07,0.02), c(0,0,0,0,0.03), c(0,0,0,0,0))
(msmtest5 <- msm(state ~ time, qmatrix = fiveq,  subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=TRUE, opt.method="fisher", control=list(trace=1))) #
msmtest5$paramdata$hess.init
(inf <- msmtest5$paramdata$info.init)
msmtest5$paramdata$allinits
solve(inf$info)%*%inf$deriv

## 3 states
bos$state3 <- bos$state
bos$state3[bos$state %in% 2:4] <- 2; bos$state3[bos$state %in% 5] <- 3
threeq <- rbind(c(0,0.01,0.002), c(0,0,0.01), c(0,0,0))
(msmtest3 <- msm(state3 ~ time, qmatrix = threeq,  subject = ptnum, data = bos, exacttimes=TRUE, fixedpars=TRUE, opt.method="optim",control=list(trace=1))) #
msmtest3$paramdata$hess.init
msmtest3$paramdata$info.init


## First ten rows of BOS data. 2 pts, 8 trans
## trans are 1-2(2), 2-2(4), 2-3(2),
bos10 <- bos[1:10,]
## derivs = sum_s nist dt(i,s)  + nijt / qij

q <- c(0.01, 0.004, 0.01)
from <- c(1,2,2,2,1,2,2,2); to <- c(2,2,2,3,2,2,2,3)
dt <- c(bos$time[2:5] - bos$time[1:4], bos$time[7:10]-bos$time[6:9])
n <- rep(1,8)
qd <- c(-(0.01+0.004), -0.01)

-2*(sum(dt*qd[from]) + sum(log(threeq[cbind(from,to)[from!=to,]]))) # loglik correct

# analytic deriv formula is
# minus sum of all timelags from that state + no of those trans / q
# so deriv wrt log q is
# minus sum of all timelags from that state * q  +  no of those trans

-2*c(q[1]*-sum(dt[c(1,5)]) + 2,  -sum(dt[c(1,5)])*q[2], -sum(dt[c(2,3,4,6,7,8)])*q[3] + 2)
# -2.3, 0.68 -2.01.
# deriv correct

## debug Fisher info obs by obs.  First for first 1-2
mod <- msm(state3~time, subject=ptnum, qmatrix=threeq, data=bos10[1:2,], exacttimes=TRUE, fixedpars=TRUE)
mod$minus2loglik
mod$paramdata$hess.init
msm(state3~time, subject=ptnum, qmatrix=threeq,data=bos10, exacttimes=TRUE, fix0dpars=FALSE, deriv.test=TRUE)



}
