\name{deltamethod}
\title{The delta method}
\alias{deltamethod}
\description{
  Delta method for approximating the standard error of a transformation
  \eqn{g(X)} of a random variable \eqn{X = (x1, x2, \ldots)}, given estimates of the
  mean and covariance matrix of \eqn{X}.
}
\usage{
deltamethod(g, mean, cov, ses=TRUE)
}
\arguments{
  \item{g}{A formula representing the transformation. It must have
    arguments labelled \code{x1, x2,\ldots} For example,
    
    \code{~ 1 / (x1 + x2)}

    If the
    transformation returns a vector, then a list of formulae g1, g2,
    \ldots can be provided, for example

    \code{list( ~ x1 + x2,    ~ x1 / (x1 + x2) )}}

    \item{mean}{The estimated mean of \eqn{X}}

    \item{cov}{The estimated covariance matrix of \eqn{X}}

    \item{ses}{If \code{TRUE}, then the standard errors of \eqn{g1(X), g2(X),\ldots} are
      returned. Otherwise the covariance matrix of \eqn{g(X)} is returned.}
  }
}
\value{
  A vector containing the standard errors of \eqn{g1(X), g2(X),\ldots}
  or a matrix containing the covariance of \eqn{g(X)}.
}
\details{
  The delta method expands a differentiable function of a random
  variable about its mean, usually with a first-order Taylor approximation,
  and then takes the variance. For example, an approximation to
  the covariance matrix of \eqn{g(X)} is given by
  
  \deqn{ Cov(g(X)) = g'(mu) Cov(X) [g'(mu)]^T }

  where \eqn{mu} is an estimate of the mean of \eqn{X}.
}
\references{Oehlert, G. W. \emph{A note on the delta method}. American
  Statistician 46(1), 1992}}
}
\examples{

## Simple linear regression, E(y) = alpha + beta x 
x <- 1:100
y <- rnorm(100, 4*x, 5)
toy.lm <- lm(y ~ x)
estmean <- coef(toy.lm)
estvar <- summary(toy.lm)$cov.unscaled

## Approximate standard error of (1 / (alphahat + betahat))
deltamethod (~ 1 / (x1 + x2), estmean, estvar) 
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{math}
