\name{sojourn.msm}
\title{Mean sojourn times from a multi-state model}
\alias{sojourn.msm}
\description{
  Estimate the mean sojourn times in the transient states of a
  multi-state model and their standard errors.
}
\usage{
sojourn.msm(x, covariates="mean")
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}}.

  \item{covariates}{
    The covariate values at which to estimate the mean sojourn times. This can either be:\cr

    the string \code{"mean"}, denoting the means of the covariates in
    the data (this is the default),\cr

    the number \code{0}, indicating that all the covariates should be
    set to zero,\cr

    a list of values, with optional names. For example,

    \code{list(60, 1)}, 
    where the order of the list follows the order of the covariates
    originally given in the model formula, or a named list, e.g.
    
    \code{list (age = 60, sex = 1)}
   
  }
  
}
\value{
    A list with components:

    \item{\code{estimate}}{Estimated mean sojourn times in the transient states.}

    \item{\code{se}}{Corresponding standard errors.}
}
\details{
  The mean sojourn time in a transient state \eqn{r} is estimated by
  \eqn{- 1 / q_{rr}}, where \eqn{q_{rr}} is the \eqn{r}th entry on the diagonal
  of the estimated transition intensity matrix.  Calls
  \code{\link{deltamethod}} to find approximate standard errors.
}
\seealso{
  \code{\link{msm}}, \code{\link{qmatrix.msm}}, \code{\link{deltamethod}}
}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{models}