source("local.R")
library(msm)

### Simulation of a multi-state model, and recovering the true model
### Just run these tests locally. 

### Simple three-state progression and death model

if (developer.local)  { 

    nsubj <- 50; nobspt <- 6
    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt), 
                         x = rnorm(nsubj*nobspt), y = rnorm(nsubj*nobspt)* 5 + 2 )
    (three.q <- rbind(c(0, exp(-3), exp(-6)), c(0, 0, exp(-3)), c(0, 0, 0)))

    sim2.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q)
    sim2.df <- simmulti.msm(sim.df, qmatrix=three.q)  # Only retains covariates if there are covariates in the simulated model

    print(crudeinits.msm(state ~ time, subject, qmatrix = rbind(c(0,1,1),c(0,0,1),c(0,0,0)), data=sim2.df))
    sim.mod <- msm(state ~ time, subject=subject, data=sim2.df,
                   qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)))
    print(sim.mod)

### Time-constant covariates 

    nsubj <- 2000; nobspt <- 15
    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt), 
                         x = rep(sample(c(0,1), nsubj, replace=TRUE), each=nobspt),
                         y = rep(sample(c(0,1), nsubj, replace=TRUE), each=nobspt)
                         )
    (three.q <- rbind(c(0, exp(-3), exp(-6)), c(0, 0, exp(-3)), c(0, 0, 0)))
    sim3.df <- simmulti.msm(sim.df, qmatrix=three.q, covariates=list(x = c(-1, 1, 0), y = c(2, 0, -2)))

    sim.mod <- msm(state ~ time, subject=subject, data=sim3.df,
                   qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)),
                   covariates = ~ x+y , covinits = list(x = c(0,0,0)), method="BFGS",
                   control=list(trace=1, REPORT=1, fnscale=10000))
    print(sim.mod)

### Time-varying covariates - much more difficult to estimate
### accurately when the covariates are constantly changing. .

    sim.df <- data.frame(subject = rep(1:nsubj, each=nobspt), time = seq(0, 20, length=nobspt), 
                         x = sample(c(0,1), nsubj*nobspt, replace=TRUE), 
                         y = sample(c(0,1), nsubj*nobspt, replace=TRUE)
                         )
    sim3.df <- simmulti.msm(sim.df, qmatrix=three.q, covariates=list(x = c(-1, 1, 0), y = c(2, 0, -2)))
    sim.mod <- msm(state ~ time, subject=subject, data=sim3.df,
                   qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)),
                   covariates = ~ x+y , covinits = list(x = c(0,0,0)), method="BFGS",
                   control=list(trace=1, REPORT=1, fnscale=10000))
    print(sim.mod)

### Exact observation times

    sim4.df <- NULL
    for (i in 1:100)  {
        s <- sim.msm(three.q, 20)
        sim4.df <- rbind(sim4.df, data.frame(subject=rep(i, length(s$states)), time=s$times, state=s$states))
    }
    print(sim4.df[1:30,])

    try(sim.mod <- msm(state ~ time, subject=subject, data=sim4.df, exacttimes=TRUE, 
                       qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)) )) # bad inits? 
    sim.mod <- msm(state ~ time, subject=subject, data=sim4.df, exacttimes=TRUE, 
                   qmatrix = rbind(c(0, exp(-2), exp(-4)), c(0, 0, exp(-2)), c(0, 0, 0)))
    print(sim.mod)
    sim.mod <- msm(state ~ time, subject=subject, data=sim4.df, exacttimes=FALSE, 
                   qmatrix = rbind(c(0, exp(-2), exp(-4)), c(0, 0, exp(-1)), c(0, 0, 0))) # worse estimates
    print(sim.mod)

### Exact death times

    simdeath.df <- simmulti.msm(sim.df[,1:2], qmatrix=three.q, death=TRUE)

    sim.mod <- msm(state ~ time, subject=subject, data=simdeath.df, death=TRUE, 
                   qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0)))
    print(sim.mod)
    sim.mod <- msm(state ~ time, subject=subject, data=simdeath.df, death=FALSE, 
                   qmatrix = rbind(c(0, exp(-1), exp(-2)), c(0, 0, exp(-1)), c(0, 0, 0))) # very slightly worse estimates of 1-3 transition
    print(sim.mod)

}
