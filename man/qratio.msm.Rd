\name{qratio.msm}
\alias{qratio.msm}
\title{Estimated ratio of transition intensities}
\description{
  Compute the estimate and approximate standard error of the ratio of
  two estimated transition intensities from a fitted multi-state model
  at a given set of covariate values. 
}
\usage{
qratio.msm(x, ind1, ind2, covariates = "mean")
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}}.

  \item{ind1}{Pair of numbers giving the indices in the intensity matrix
    of the numerator of the ratio, for example, \code{c(1,2)}. }

  \item{ind2}{Pair of numbers giving the indices in the intensity matrix
    of the denominator of the ratio, for example, \code{c(2,1)}. }

  \item{covariates}{
    The covariate values at which to estimate the intensities. 
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
}
\value{
  A list with components \code{estimate} and \code{se}, containing the
  estimate and standard error respectively of the ratio of intensities. 
}
\details{

  For example, we might want to compute the ratio of the progression
  rate and recovery rate for a fitted model \code{disease.msm} with a
  health state (state 1) and a disease state
  (state 2).  In this case, the progression rate is the (1,2) entry of
  the intensity matrix, and the recovery rate is the (2,1) entry.
  Thus to compute this ratio with covariates set to their means, we
  call

  \code{qratio.msm(disease.msm, c(1,2), c(2,1))} .

  Standard errors are estimated by the delta method. 
  
}
\seealso{
  \code{\link{qmatrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
