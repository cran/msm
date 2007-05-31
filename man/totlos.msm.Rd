\name{totlos.msm}
\alias{totlos.msm}
\title{Total length of stay}
\description{
  Estimate the expected total length of stay in each transient state,
  for a given period of evolution of a multi-state model. This assumes that the
  transition rates do not change with time. 
}
\usage{
totlos.msm(x, start=1, fromt=0, tot=Inf, covariates="mean", ci.boot=FALSE, cl=0.95, B=500, ...)
}
\arguments{
  \item{x}{A fitted multi-state model, as returned by
    \code{\link{msm}}. }
  
  \item{start}{State at the beginning of the period.}

  \item{fromt}{Time from which to estimate total length of
    stay. Defaults to 0, the beginning of the process.}
  
  \item{tot}{Time up to which total length of stay is estimated. Defaults
    to infinity, giving the expected time spent in the state until
    absorption. For models without an absorbing state, \code{t}
    must be specified.}

  \item{covariates}{
    The covariate values to estimate for.  This can either be:\cr

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
  \item{ci.boot}{Calculate a bootstrap confidence interval. This is usually time-consuming, and disabled by default. See \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.}
  \item{cl}{Width of the symmetric confidence interval}
  \item{B}{Number of bootstrap replicates}

  \item{...}{Further arguments to be passed to the
    \code{\link{integrate}} function to control the numerical
    integration.}
}
\value{
  A vector of expected total lengths of stay for each transient state.
}
\details{
  The expected total length of stay in state \eqn{j} between times
  \eqn{t_1} and \eqn{t_2}, from the point of view of an individual in
  state \eqn{i} at time 0,
  is defined by the integral from \eqn{t_1} to
  \eqn{t_2} of the \eqn{i,j} entry of the transition probability matrix \eqn{P(t)}.
  As the individual entries of \eqn{P(t) = \exp(tQ)}{P(t) = exp(tQ)}
  are not available explicitly in terms of
  \eqn{t} for a general Markov model, this integral is calculated
  numerically, using the \code{\link{integrate}} function.  This may
  take a long time for models with many states where \eqn{P(t)} is
  expensive to calculate. 

  For a model where the individual has only one place to go from each
  state, and each state is visited only once, for example a progressive
  disease model with no recovery or death, these are equal to the mean
  sojourn time in each state.
  However, consider a three-state health-disease-death model with
  transitions from health to disease, health to death, and disease to
  death, where everybody starts healthy.  In this case the mean sojourn
  time in the disease state will be greater than the expected length of
  stay in the disease state.  This is because the mean sojourn time in a
  state is conditional on entering the state, whereas the
  expected total time diseased is a forecast for a healthy
  individual,  who may die before getting the disease. 

}
\seealso{
  \code{\link{sojourn.msm}}, \code{\link{pmatrix.msm}}, \code{\link{integrate}}, \code{\link{boot.msm}}.
}
\author{C. H. Jackson \email{chris.jackson@mrc-bsu.cam.ac.uk}}
\keyword{}
