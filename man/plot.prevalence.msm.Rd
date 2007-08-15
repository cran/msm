\name{plot.prevalence.msm}
\alias{plot.prevalence.msm}
\title{Plot of observed and expected prevalences}
\description{
  Provides a rough indication of goodness of fit of a multi-state model,
  by estimating the observed numbers of individuals occupying a state
  at a series of times, and plotting these against forecasts from the
  fitted model, for each state.  Observed prevalences are indicated as
  solid lines, expected prevalences as dashed lines. 
}
\usage{
plot.prevalence.msm(x, mintime, maxtime, interp, ...)
}
\arguments{
    \item{x}{A fitted multi-state model produced by \code{\link{msm}}.}
    \item{mintime}{Minimum time at which to compute the observed and
      expected prevalences of states.}
    \item{maxtime}{Maximum time at which to compute the observed and
      expected prevalences of states.}
    \item{interp}{Interpolation method for observed states, see \code{\link{prevalence.msm}}}
    \item{...}{Further arguments to be passed to the generic \code{\link{plot}}
      function.}
}
\value{
}
\details{
  See \code{\link{prevalence.msm}} for details of the assumptions
  underlying this method.
}
\references{
  Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
  Markov models for analysing incomplete disease history data with
  illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
  805--821.
}
\seealso{
  \code{\link{prevalence.msm}}
}
\keyword{models}
