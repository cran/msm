\name{prevalence.msm}
\title{Tables of observed and expected prevalences}
\alias{prevalence.msm}
\description{
  This provides a rough indication of the goodness of fit of a
  multi-state model, by estimating the observed numbers of individuals
  occupying each state at a series of times, and comparing these with
  forecasts from the fitted model.   
}
\usage{
prevalence.msm(x, times, timezero=min(time), initstates, covariates="mean", misccovariates="mean")
}
\arguments{
  \item{x}{A fitted multi-state model, either with or without misclassification, produced
    by \code{\link{msm}}. The data for the fitted model must originally
    have been provided as a series of states and observation times, not in
    'from-state, to-state, time-difference' format.}
  \item{times}{Series of times at which to compute the observed and
    expected prevalences of states.}
  \item{timezero}{Initial time of the Markov process. Expected values
    are forecasted from here. Defaults to the minimum of the observation
    times given in the data. }
  \item{initstates}{Optional vector of the same length as the number of
    states. Gives the numbers of individuals occupying each state at the
    initial time. The default is those observed in the data. }
  \item{covariates}{Covariate values for which to forecast expected
    state occupancy.  See \code{\link{qmatrix.msm}}.  Defaults to the
    mean values of the covariates in the data set.}
  \item{misccovariates}{(Misclassification models only) Values of covariates on the misclassification
    probability matrix for which to forecast expected state occupancy.
    Defaults to the mean values of the covariates in the data set.}
}
\value{
  A list with components:
  
  \item{\code{Observed}}{Table of observed numbers of individuals in each state at
    each time}
  
  \item{\code{Observed percentages}}{Corresponding percentage of the
    individuals at risk at each time.}
  
  \item{\code{Expected}}{Table of corresponding expected numbers.}
  
  \item{\code{Expected percentages}}{Corresponding percentage of the
    individuals at risk at each time.}
  
}

\details{
  To compute `observed'  prevalences at a time \eqn{t}, individuals are
  assumed to be in the same state as at their last observation time preceding
  \eqn{t}.  This calculation is rather slow when there are a large
  number of individuals in the data.  

  The fitted transition probability matrix is used to forecast expected prevalences from the
  state occupancy at the initial time.  To
  produce the expected number in state \eqn{j} at time \eqn{t} after the
  start, the number of individuals under observation at time \eqn{t}
  (including those who have died, but not those lost to follow-up) is
  multiplied by the probability of transition between the initial state
  and state \eqn{j} in a time interval \eqn{t}.

  For misclassification models, this aims to assess the fit of the model for the \emph{observed}
  states to the data,  that is the combined Markov progression model for the true states and
  the misclassification model. Thus, expected prevalences of \emph{true}
  states are estimated from the assumed proportion
  occupying each state at the initial time using the fitted transition
  probabiliy matrix. The vector of expected prevalences of true states
  is then multiplied by the fitted misclassification probability matrix
  to obtain the expected prevalences of \emph{observed} states. 

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
  \code{\link{msm}}, \code{\link{summary.msm}}}
}

\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
