\name{print.msm}
\title{Print a fitted multi-state model}
\alias{print.msm}
\description{
  The fitted transition intensity matrix, matrices containing
  covariate effects on intensities, and mean sojourn times from a fitted
  \code{\link{msm}} model are pretty printed, with their corresponding standard errors.
}
\usage{
\method{print}{msm}(x,\dots)
}
\arguments{
  \item{x}{output from \code{\link{msm}} representing a fitted
    multi-state model}
  \item{...}{further arguments passed to or from other methods}
}
\seealso{
  \code{\link{msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
