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
### adapted from mexp in Jim Lindsey's rmutil library

MatrixExp <- function(mat, t = 1, n = 20, k = 3)
  {
      ev <- eigen(mat)
      if (any ( duplicated(ev$values) ) ) {
          ## series approximation
          mat <- mat*t / 2^k
          sum <- power <- diag(dim(mat)[2])
          for (r in 1:n) {
              power <- mat %*% power / r
              sum <- sum + power
          }
          for (i in 1:k)
            sum <- sum %*% sum
          sum
      }
      else
        ## spectral decomposition
        ev$vectors %*% diag(exp(ev$values * t)) %*% solve(ev$vectors)
  }
