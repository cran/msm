\name{plot.msm}
\title{Plots of multi-state models}
\alias{plot.msm}
\description{
  This produces a plot of the expected probability of
  survival against time, from each transient state. Survival is defined
  as not entering the final absorbing state.
}
\usage{
\method{plot}{msm}(x, from, to, legend.pos, ...)
}
\arguments{

  \item{x}{Output from \code{\link{msm}}, representing a fitted
    multi-state model object}

  \item{from}{States from which to consider survival from. Defaults to
    the complete set of transient states. }

  \item{to}{Absorbing state to consider. Defaults to the
    highest-labelled absorbing state. }

  \item{legend.pos}{Vector of the \eqn{x} and \eqn{y} position,
    respectively, of the legend.}

    \item{...}{Other arguments to the generic \code{plot} function}
}
\seealso{
    \code{\link{msm}}
}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{models}