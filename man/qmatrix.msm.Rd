\name{qmatrix.msm}
\alias{qmatrix.msm}
\title{Transition intensity matrix}
\description{
  Extract the estimated transition intensity matrix, and the
  corresponding standard errors, from a fitted
  multi-state model at a given set of covariate values. 
}
\usage{
qmatrix.msm(x, covariates="mean", sojourn=FALSE)
}
\arguments{
 
  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}}.

  \item{covariates}{
    The covariate values at which to estimate the intensity matrix. 
    This can either be:\cr

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
  \item{sojourn}{Set to TRUE if the estimated sojourn times and their
    standard errors should also be returned.}
}

\value{
  A list with components:
  
  \item{estimate}{Estimated transition intensity matrix.}
  
  \item{SE}{Corresponding approximate standard errors.}  
  If \code{sojourn} is \code{TRUE}, extra components called
  \code{sojourn} and \code{sojournSE} are included, containing the
  estimate and standard errors, respectively, of the mean sojourn times in
  each transient state.

}

\details{
  Transition intensities and covariate effects are estimated on the log
  scale by \code{\link{msm}}. A covariance matrix is estimated from the
  Hessian of the maximised log-likelihood.  The delta method is used to
  obtain from these the standard error of the intensities on the natural
  scale at arbitrary covariate values. 
}

\seealso{
  \code{\link{pmatrix.msm}}, \code{\link{sojourn.msm}},
  \code{\link{deltamethod}}, \code{\link{ematrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
