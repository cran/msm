\name{sojourn.msm}
\title{Mean sojourn times from a multi-state model}
\alias{sojourn.msm}
\description{
  Estimate the mean sojourn times in the transient states of a
  multi-state model and their confidence limits. 
}
\usage{
sojourn.msm(x, covariates="mean", cl=0.95)
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
  
  \item{cl}{Width of the symmetric confidence interval to present. 
    Defaults to 0.95.}
  
}
\value{
    A data frame with components:

    \item{estimates}{Estimated mean sojourn times in the transient states.}
    \item{SE}{Corresponding standard errors.}
    \item{L}{Lower confidence limits.}
    \item{U}{Upper confidence limits.}
}
\details{
  The mean sojourn time in a transient state \eqn{r} is estimated by
  \eqn{- 1 / q_{rr}}, where \eqn{q_{rr}} is the \eqn{r}th entry on the diagonal
  of the estimated transition intensity matrix.  Calls
  \code{\link{deltamethod}} to find approximate standard errors.
  Confidence limits are estimated by assuming normality on the log
  scale. 
}
\seealso{
  \code{\link{msm}}, \code{\link{qmatrix.msm}}, \code{\link{deltamethod}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
