### PACKAGE GLOBAL CONSTANTS ###

### List of allowed hidden Markov model distributions
### and names of parameters for each distribution
### MUST BE KEPT IN THE SAME ORDER as the C variable HMODELS in lik.c
.msm.HMODELPARS <- list(
                        categorical=c("ncats","basecat","prob"),  # vector 
                        identity = NULL,
                        uniform = c("lower", "upper"),
                        normal = c("mean", "sd"),
                        lognormal = c("meanlog", "sdlog"),
                        exponential = c("rate"),
                        gamma = c("shape","rate"),
                        weibull = c("shape","scale"),
                        poisson = c("rate"),
                        binomial = c("size","prob"),
                        truncnorm = c("mean", "sd", "lower", "upper"),
                        metruncnorm = c("mean", "sd", "lower", "upper", "sderr", "meanerr"),
                        meuniform = c("lower", "upper", "sderr", "meanerr"),
                        nbinom = c("disp","prob")
                        )

## TODO - beta, non-central beta, cauchy, chisq, noncentral chisq, F,
## non-central F, geometric, hypergeometric, logistic, negative
## binomial, t, noncentral t.

.msm.HMODELS <- names(.msm.HMODELPARS)

### Parameter in each distribution that can have covariates on it 
.msm.LOCPARS <- c(categorical="p", identity=NA, uniform=NA, normal="mean", lognormal="meanlog",
                  exponential="rate", gamma="rate", weibull="scale", 
                  poisson="rate", binomial="prob", truncnorm="mean",
                  metruncnorm="meanerr", meuniform="meanerr", nbinom="prob")

### Link functions for generalised regressions.
### MUST BE KEPT IN SAME ORDER as LINKFNS in lik.c
.msm.LINKFNS <- c("identity", "log", "qlogis")                    
.msm.INVLINK <- c(identity="identity", log="exp", qlogis="plogis")
                  
### Parameters which are always fixed, never estimated
.msm.AUXPARS <- c("lower", "upper", "which", "size", "meanerr", "ncats", "basecat", "p0", "pbase")

### Parameters which should be defined as integer
.msm.INTEGERPARS <- c("size")

### Defined ranges for parameters
.msm.PARRANGES <- list(qbase=c(0, Inf), p=c(0, 1), lower=c(-Inf,Inf), upper=c(-Inf, Inf),
                       mean=c(-Inf, Inf), sd=c(0, Inf), 
                       meanlog=c(-Inf,Inf), sdlog=c(0, Inf), rate=c(0, Inf), shape=c(0, Inf),
                       prob=c(0, 1), meanerr=c(0, Inf), sderr=c(0, Inf), disp=c(0, Inf))

### Transforms to optimise some parameters on a different scale
.msm.TRANSFORMS <-
  do.call("rbind", 
          lapply(.msm.PARRANGES,
                 function(x) {
                     if (identical(x, c(0, Inf))) c(fn="log",inv="exp")
                     else if (identical(x, c(0, 1))) c(fn="qlogis",inv="plogis")
                     else NULL
                 }
                 ) )

