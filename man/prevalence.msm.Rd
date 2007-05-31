\name{prevalence.msm}
\title{Tables of observed and expected prevalences}
\alias{prevalence.msm}
\concept{Goodness of fit}
\description{
  This provides a rough indication of the goodness of fit of a
  multi-state model, by estimating the observed numbers of individuals
  occupying each state at a series of times, and comparing these with
  forecasts from the fitted model.   
}
\usage{
prevalence.msm(x, times, timezero=NULL, initstates, covariates="mean",
               misccovariates="mean", piecewise.times, piecewise.covariates,
	       ci.boot=FALSE, cl=0.95, B=500)
}
\arguments{
  \item{x}{A fitted multi-state model produced by \code{\link{msm}}.}
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
  \item{piecewise.times}{Times at which piecewise-constant intensities
    change. See \code{\link{pmatrix.piecewise.msm}} for how
    to specify this.}
  \item{piecewise.covariates}{Covariates on which the piecewise-constant
    intensities depend. See \code{\link{pmatrix.piecewise.msm}} for how
    to specify this.}
  \item{ci.boot}{Calculate a bootstrap confidence interval for the expected prevalences. This is usually time-consuming, and disabled by default. See \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.}
  \item{cl}{Width of the symmetric confidence interval}
  \item{B}{Number of bootstrap replicates}
}
\value{
  A list of matrices, with components:
  
  \item{Observed}{Table of observed numbers of individuals in each state at
    each time}
  
  \item{Observed percentages}{Corresponding percentage of the
    individuals at risk at each time.}
  
  \item{Expected}{Table of corresponding expected numbers.}
  
  \item{Expected percentages}{Corresponding percentage of the
    individuals at risk at each time.}

  Or if \code{ci.boot = TRUE}, the component \code{Expected} is a list
  with components \code{estimates} and \code{ci}. \code{estimates} is
  a matrix of the expected prevalences, and \code{ci} is a list of two
  matrices, containing the confidence limits.  The component
  \item{Expected percentages} has a similar format.

}

\details{
  To compute `observed'  prevalences at a time \eqn{t}, individuals are
  assumed to be in the same state as at their last observation time preceding
  \eqn{t}.  

  The fitted transition probability matrix is used to forecast expected prevalences from the
  state occupancy at the initial time.  To
  produce the expected number in state \eqn{j} at time \eqn{t} after the
  start, the number of individuals under observation at time \eqn{t}
  (including those who have died, but not those lost to follow-up) is
  multiplied by the product of the proportion of individuals in each
  state at the initial time and the transition probability matrix in the
  time interval \eqn{t}.   The proportion of individuals in each state
  at the "initial" time is estimated, if necessary, in the same way as
  the observed prevalences. 

  For misclassification models (fitted using an \code{ematrix}), this
  aims to assess the fit of the full model for the \emph{observed}
  states.  That is, the combined Markov progression model for the true
  states and the misclassification model. Thus, expected prevalences of \emph{true}
  states are estimated from the assumed proportion
  occupying each state at the initial time using the fitted transition
  probabiliy matrix. The vector of expected prevalences of true states
  is then multiplied by the fitted misclassification probability matrix
  to obtain the expected prevalences of \emph{observed} states. 

  For general hidden Markov models, the observed state is taken to be the
  predicted underlying state from the Viterbi algorithm
  (\code{\link{viterbi.msm}}).  The goodness of fit of
  these states to the underlying Markov model is tested. 

  Note that this function assumes intensities are the same for all
  individuals. By default they are taken from the mean values of the
  covariates in the population.  Piecewise-constant intensities may be
  assumed, through the arguments \code{piecewise.times} and
  \code{piecewise.covariates}.
  
  For an example of this approach, see Gentleman \emph{et al.} (1994).
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

\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{models}
