\name{pmatrix.msm}
\alias{pmatrix.msm}
\title{Transition probability matrix}
\description{
  Extract the estimated transition probability matrix from a fitted multi-state
  model for a given time interval, at a given set of covariate values. 
}
\usage{
pmatrix.msm(x, t, covariates="mean")
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}

  \item{t}{The time interval to estimate the transition probabilities for.}

  \item{covariates}{
    The covariate values at which to estimate the transition
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
}
\value{
  The matrix of estimated transition probabilities \eqn{P(t)} in the given time. 
  Rows correspond to "from-state" and columns to "to-state". 
}
\details{
  For a continuous-time homogeneous Markov process with transition
  intensity matrix
  \eqn{Q},  the probability of occupying state \eqn{s} at time \eqn{u + t}
  conditionally on occupying state \eqn{r} at time \eqn{u} is given by the
  \eqn{(r,s)} entry of the matrix \eqn{P(t) = \exp(tQ)}{P(t) = exp(tQ)}.

  For non-homogeneous processes, where covariates and hence the
  transition intensity matrix are time-dependent, but are
  piecewise-constant within  the time interval \code{[u,
  u+t]}, the function \code{\link{pmatrix.piecewise.msm}} can be used.
  
}
\seealso{
  \code{\link{qmatrix.msm}},  \code{\link{pmatrix.piecewise.msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{}
