\name{plot.msm}
\title{Plots of multi-state models}
\alias{plot.msm}
\description{
  This produces a plot of the expected probability of
  survival against time, from each transient state. Survival is defined
  as not entering the final absorbing state.
}
\usage{
\method{plot}{msm}(x, from, to, covariates, legend.pos, ...)
}
\arguments{

  \item{x}{Output from \code{\link{msm}}, representing a fitted
    multi-state model object}

  \item{from}{States from which to consider survival from. Defaults to
    the complete set of transient states. }

  \item{to}{Absorbing state to consider. Defaults to the
    highest-labelled absorbing state. }

  \item{covariates}{Covariate values for which to evaluate the expected
    probabilities.  This can either be:\cr

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

  \item{legend.pos}{Vector of the \eqn{x} and \eqn{y} position,
    respectively, of the legend.}
  
  \item{...}{Other arguments to the generic \code{plot} function}
}
\seealso{
  \code{\link{msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
