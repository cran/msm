\name{printnew.msm}
\title{Print a fitted msm model object}
\alias{printnew.msm}
\description{
  Print a fitted msm model object
}
\usage{
printnew.msm(x, covariates=NULL, digits=4, \dots)
}
\arguments{

  \item{x}{Output from \code{\link{msm}}, representing a fitted
    multi-state model object.}

  \item{covariates}{Covariates for which to print ``baseline''
  transition intensities or misclassification probabilities. See
  \code{\link{qmatrix.msm}} for more details. }

  \item{digits}{Minimum number of significant digits, passed to
    \code{\link{format}}. Defaults to 4.}

  \item{...}{Other arguments to be passed to \code{\link{format}}.}
}
\details{
  This is the new method of formatting msm objects for printing.
  The old method was based on printing lists of matrices,
This produced a lot of wasted space for parameters which were zero,
and it was difficult to match corresponding numbers between matrices.
The new method presents all the transition intensities and covariate
effects as a single compact table, and likewise for misclassification
matrices.

Also in the old method, covariate effects were presented as log hazard ratios or
log odds ratios.  The log scale is more convenient mathematically but
unnatural for a human to interpret.  The new method presents hazard
ratios for covariates on transition intensities and odds ratios for
misclassification probabilities.

The old method \code{\link{print.msm}}, is still the default, but the
  new method is intended to be the default in future versions.
  }
\seealso{
  \code{\link{msm}}
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
