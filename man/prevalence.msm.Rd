\name{prevalence.msm}
\title{Calculate tables of observed and expected prevalences at a series
of times}
\alias{prevalence.msm}
\description{
  This function is called by \code{\link{summary.msm}},  and can provide a
  rough indication of the goodness of fit of a multi-state model without
  misclassification. 
}
\usage{
prevalence.msm(msm, times, timezero=min(time), initstates)
}
\arguments{
  \item{msm}{A fitted multi-state model without misclassification, produced
    by \code{\link{msm}}. The data for the fitted model must originally
    have been provided as a series of states and observation times, not in
    'from-state, to-state, time-difference' format.}
  \item{times}{Series of times at which to compute the observed and
    expected prevalences of states.}
  \item{timezero}{Initial time of the Markov process. Expected values
    are forecasted from here. Defaults to the minimum of the observation
    times given in the data. }
  \item{initstates}{Optional vector of the same length as the number of
    states. Gives the numbers of individuals occupying each state. The
    default is that all individuals are in state 1 at
    \code{timezero}. }
}
\value{
  A list with components:
  
  \code{Observed} = table of observed numbers of individuals in each state at
  each time

  \code{Observed percentages} = corresponding percentage of the
  individuals at risk at each time  

  \code{Expected} = table of corresponding expected numbers

  \code{Expected percentages} = corresponding percentage of the
  individuals at risk at each time.
  
}

\details{
  For models without misclassification only. To compute `observed'
  prevalences at a time \eqn{t}, individuals are assumed to be in the
  same state as at their last observation time preceding
  \eqn{t}.

  The fitted transition probability matrix is used to forecast expected prevalences from the
  state occupancy at the initial time.  To
  produce the expected number in state \eqn{j} at time \eqn{t} after the
  start, the number of individuals under observation at time \eqn{t}
  (including those who have died, but not those lost to follow-up) is
  multiplied by the probability of transition between the initial state
  and state \eqn{j} in a time interval \eqn{t}.

  This approach only makes sense for processes where all individuals
  start at a common time.  
  
  For an example of this approach, see Gentleman \emph{et
    al.} (1994).

}

\references{
  Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
  Markov models for analysing incomplete disease history data with
  illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
  805--821.
}

\seealso{
  \code{\link{msm}}, \code{\link{summary.msm}}, \code{\link{prevalencemisc.msm}} for the equivalent for models with
    misclassification}
}

\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{models}