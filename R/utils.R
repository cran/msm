### msm PACKAGE
### USEFUL FUNCTIONS NOT SPECIFIC TO MULTI-STATE MODELS

### Delta method for approximating the covariance matrix of f(X) given cov(X)

deltamethod <- function(g,       # a formula or list of formulae (functions) giving the transformation g(x) in terms of x1, x2,  etc
                        mean,    # mean, or maximum likelihood estimate, of x
                        cov,     # covariance matrix of x
                        ses=TRUE # return standard errors, else return covariance matrix
                        )
  {
      ## Var (G(x))  =  G'(mu) Var(X) G'(mu)^T
      cov <- as.matrix(cov)
      n <- length(mean)
      if (!is.list(g))
        g <- list(g)
      if ( (dim(cov)[1] != n) || (dim(cov)[2] != n) )
        stop(paste("Covariances should be a ", n, " by ", n, " matrix"))
      syms <- paste("x",1:n,sep="")
      for (i in 1:n)
        assign(syms[i], mean[i])
      gdashmu <- t(sapply(g,
                          function( form ) {
                              as.numeric(attr(eval(
                                                   ## Differentiate each formula in the list
                                                   deriv(form, syms)
                                                   ## evaluate the results at the mean
                                                   ), "gradient"))
                              ## and build the results row by row into a Jacobian matrix
                          }))
      new.covar <- gdashmu %*% cov %*% t(gdashmu)
      if (ses){
          new.se <- sqrt(diag(new.covar))
          new.se
      }
      else
        new.covar
  }

### Matrix exponential

MatrixExp <- function(mat, t = 1, n = 20, k = 3, method="pade")
  {
      if (!is.matrix(mat) || (nrow(mat)!= ncol(mat))) stop("\"mat\" must be a square matrix")
      n <- nrow(mat)
      ev <- eigen(mat)
      if (any ( duplicated(ev$values) ) ) {
          if (method=="series") {
              ## series approximation
              ## adapted from mexp in Jim Lindsey's rmutil library
              mat <- mat*t / 2^k
              sum <- power <- diag(dim(mat)[2])
              for (r in 1:n) {
                  power <- mat %*% power / r
                  sum <- sum + power
              }
              for (i in 1:k)
                sum <- sum %*% sum
              res <- sum
          }
          else if (method == "pade") {
              ## C function adapted from JAGS by Martyn Plummer 
              res <- .C("MatrixExpPadeR", res=double(length(mat)), as.double(mat),
                        as.integer(n), as.double(t), PACKAGE="msm")$res
              res <- matrix(res, nrow=nrow(mat))
          }
          else stop("Method should be \"pade\" or \"series\"")
      }
      else
        ## spectral decomposition
        res <- ev$vectors %*% diag(exp(ev$values * t)) %*% solve(ev$vectors)
      res
  }

identity <- function(x)x

### Truncated normal distribution

dtnorm <- function(x, mean=0, sd=1, lower=-Inf, upper=Inf, log=FALSE)
  {
      ret <- numeric(length(x))
      ret[x < lower | x > upper] <- 0
      ind <- x >=lower & x <=upper
      if (any(ind)) {
          denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, sd)
          xtmp <- dnorm(x, mean, sd, log)
          if (log) xtmp <- xtmp - log(denom) else xtmp <- xtmp/denom
          ret[x >=lower & x <=upper] <- xtmp[ind]
      }
      ret
  }

ptnorm <- function(q, mean=0, sd=1, lower=-Inf, upper=Inf, lower.tail=TRUE, log.p=FALSE)
  {
      ret <- numeric(length(q))
      ret[q < lower] <- 0
      ret[q > upper] <- 1
      ind <- q >=lower & q <=upper
      if (any(ind)) {
          denom <- pnorm(upper, mean, sd) - pnorm(lower, mean, sd)
          if (lower.tail) qtmp <- pnorm(q, mean, sd) - pnorm(lower, mean, sd)
          else qtmp <- pnorm(upper, mean, sd) - pnorm(q, mean, sd)
          if (log.p) qtmp <- log(qtmp) - log(denom) else qtmp <- qtmp/denom
          ret[q >=lower & q <=upper] <- qtmp[ind]
      }
      ret
  }

qtnorm <- function(p, mean=0, sd=1, lower=-Inf, upper=Inf, lower.tail=TRUE, log.p=FALSE)
  {
      ret <- numeric(length(p))
      ret[p == 1] <- Inf
      ret[p == 0] <- -Inf
      ret[p < 0 | p > 1] <- NaN
      ind <- (p > 0 & p < 1)
      if (any(ind)) {
          hind <- seq(along=p)[ind]
          h <- function(y) {
              (ptnorm(y, mean, sd, lower, upper) - p)[hind[i]]
          }
          ptmp <- numeric(length(p[ind]))
          for (i in 1:length(p[ind])) {
              interval <- c(-1, 1)
              while (h(interval[1])*h(interval[2]) >= 0)
                interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
              ptmp[i] <- uniroot(h, interval)$root
          }
          ret[ind] <- ptmp
      }
      if (any(is.nan(ret))) warning("NaNs produced")
      ret
  }

rtnorm <- function(n, mean=0, sd=1, lower=-Inf, upper=Inf)
  {
      ret <- numeric()
      if (length(n) > 1) n <- length(n)
      while (length(ret) < n) {
          y <- rnorm(n - length(ret), mean, sd)
          y <- y[y>=lower & y<=upper]
          ret <- c(ret, y)
      }
      stopifnot(length(ret)==n)
      ret
  }

### Normal distribution with measurement error and optional truncation 

dmenorm <- function(x, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, log = FALSE)
{
    sumsq <- sd*sd + sderr*sderr
    sigtmp <- sd*sderr / sqrt(sumsq)
    mutmp <- ((x - meanerr)*sd*sd + mean*sderr*sderr) / sumsq
    nc <- 1/(pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
    nctmp <- pnorm(upper, mutmp, sigtmp) - pnorm(lower, mutmp, sigtmp)
    nc * nctmp * dnorm(x, meanerr + mean, sqrt(sumsq), 0) 
}

pmenorm <- function(q, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, lower.tail = TRUE, log.p = FALSE)
{
    ret <- numeric(length(q))
    dmenorm2 <- function(x)dmenorm(x, mean=mean, sd=sd, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr)
    for (i in 1:length(q)) { 
        ret[i] <- integrate(dmenorm2, -Inf, q[i])$value
    }
    ret
}

qmenorm <- function(p, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0, lower.tail = TRUE, log.p = FALSE)
{
      ret <- numeric(length(p))
      ret[p == 1] <- Inf
      ret[p == 0] <- -Inf
      ret[p < 0 | p > 1] <- NaN
      ind <- (p > 0 & p < 1)
      if (any(ind)) {
          hind <- seq(along=p)[ind]
          h <- function(y) {
              (pmenorm(y, mean, sd, lower, upper, sderr, meanerr) - p)[hind[i]]
          }
          ptmp <- numeric(length(p[ind]))
          for (i in 1:length(p[ind])) {
              interval <- c(-1, 1)
              while (h(interval[1])*h(interval[2]) >= 0)
                interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
              ptmp[i] <- uniroot(h, interval)$root
          }
          ret[ind] <- ptmp
      }
      if (any(is.nan(ret))) warning("NaNs produced")
      ret
}

rmenorm <- function(n, mean=0, sd=1, lower=-Inf, upper=Inf, sderr=0, meanerr=0)
{
    rnorm(n, meanerr + rtnorm(n, mean, sd, lower, upper), sderr)
}

### Uniform distribution with measurement error

dmeunif <- function(x, lower=0, upper=1, sderr=0, meanerr=0, log = FALSE)
{
    ( pnorm(x, meanerr + lower, sderr) - pnorm(x, meanerr + upper, sderr) ) / (upper - lower)
}

pmeunif <- function(q, lower=0, upper=1, sderr=0, meanerr=0, lower.tail = TRUE, log.p = FALSE)
{
    ret <- numeric(length(q))
    dmeunif2 <- function(x)dmeunif(x, lower=lower, upper=upper, sderr=sderr, meanerr=meanerr)
    for (i in 1:length(q)) { 
        ret[i] <- integrate(dmeunif2, -Inf, q[i])$value
    }
    ret
}

qmeunif <- function(p, lower=0, upper=1, sderr=0, meanerr=0, lower.tail = TRUE, log.p = FALSE)
{
      ret <- numeric(length(p))
      ret[p == 1] <- Inf
      ret[p == 0] <- -Inf
      ret[p < 0 | p > 1] <- NaN
      ind <- (p > 0 & p < 1)
      if (any(ind)) {
          hind <- seq(along=p)[ind]
          h <- function(y) {
              (pmeunif(y, lower, upper, sderr, meanerr) - p)[hind[i]]
          }
          ptmp <- numeric(length(p[ind]))
          for (i in 1:length(p[ind])) {
              interval <- c(-1, 1)
              while (h(interval[1])*h(interval[2]) >= 0)
                interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
              ptmp[i] <- uniroot(h, interval)$root
          }
          ret[ind] <- ptmp
      }
      if (any(is.nan(ret))) warning("NaNs produced")
      ret
}

rmeunif <- function(n, lower=0, upper=1, sderr=0, meanerr=0)
{
    rnorm(n, meanerr + runif(n, lower, upper), sderr)
}
