\name{ematrix.msm}
\alias{ematrix.msm}
\title{Misclassification probability matrix}
\description{
  Extract the estimated misclassification probability matrix, and the
  corresponding standard errors,  from a fitted multi-state
  model at a given set of covariate values. 
}
\usage{
ematrix.msm(x, covariates="mean")
}
\arguments{
  
  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}}.

  \item{covariates}{

    The covariate values for which to estimate the misclassification
    probability matrix.  This can either be:\cr

    the string \code{"mean"}, denoting the means of the covariates in
    the data (this is the default),\cr

    the number \code{0}, indicating that all the covariates should be
    set to zero,\cr

    or a list of values, with optional names. For example

    \code{list (60, 1)}

    where the order of the list follows the order of the covariates
    originally given in the model formula, or a named list, 

    \code{list (age = 60, sex = 1)}
  }
 
}
\value{
  A list with components:

  \item{\code{estimate}}{Estimated misclassification probability matrix.}  
  \item{\code{SE}}{Corresponding approximate standard errors.}

}
\details{
  Misclassification probabilities and covariate effects are estimated on the logit
  scale by \code{\link{msm}}. A covariance matrix is estimated from the
  Hessian of the maximised log-likelihood.  The delta method is used to
  obtain from these the standard error of the probabilities on the
  natural scale at arbitrary covariate values. 
}
\seealso{
  \code{\link{qmatrix.msm}}
}
\examples{
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
