\name{ematrix.msm}
\alias{ematrix.msm}
\title{Misclassification probability matrix}
\description{
  Extract the estimated misclassification probability matrix, and 
  corresponding confidence intervals,  from a fitted multi-state
  model at a given set of covariate values. 
}
\usage{
ematrix.msm(x, covariates="mean", cl=0.95)
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
 
  \item{cl}{Width of the symmetric confidence interval to present. 
    Defaults to 0.95.}
}
\value{
  A list with components:

  \item{estimate}{Estimated misclassification probability matrix.}  
  \item{SE}{Corresponding approximate standard errors.}
  \item{L}{Lower confidence limits.}
  \item{U}{Upper confidence limits.}

  The default print method for objects returned by
  \code{\link{ematrix.msm}} presents estimates and confidence limits. To
  present estimates and standard errors, do something like 

  \code{ematrix.msm(x)[c("estimates","SE")]}

}
\details{
  Misclassification probabilities and covariate effects are estimated on the logit
  scale by \code{\link{msm}}. A covariance matrix is estimated from the
  Hessian of the maximised log-likelihood.  From these, the delta method
  is used to obtain standard errors of the probabilities on the natural
  scale at arbitrary covariate values.  Confidence intervals are
  estimated by assuming normality on the logit scale. 
}
\seealso{
  \code{\link{qmatrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
