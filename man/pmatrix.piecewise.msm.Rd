\name{pmatrix.piecewise.msm}
\alias{pmatrix.piecewise.msm}
\title{Transition probability matrix for processes with piecewise-constant intensities}
\description{
  Extract the estimated transition probability matrix from a fitted multi-state
  model for a given time interval.  This is a generalisation of
  \code{\link{pmatrix.msm}} to non-time-homogeneous Markov models with
  piecewise-constant transition intensity matrices. 
}
\usage{
pmatrix.piecewise.msm(x, t1, t2, times, covariates)
}
\arguments{

  \item{x}{A fitted multi-state model, as returned by \code{\link{msm}}.}

  \item{t1}{The start of the time interval to estimate the transition probabilities for.}

  \item{t2}{The end of the time interval to estimate the transition probabilities
    for.}

  \item{times}{Cut points at which the transition intensity matrix changes.}

  \item{covariates}{
    A list with number of components one greater than the length of
    \code{times}.  Each component of the list is specified in the same
    way as the \code{covariates} argument to \code{\link{pmatrix.msm}}.
    The components correspond to the covariate values in the intervals
    
    \code{(t1, times[1]],  (times[1], times[2]], ...,
      (\times[length(times)], t2]}

    (assuming that all elements of \code{times} are in the interval
    \code{(t1, t2)}).
  }
}
\value{
  The matrix of estimated transition probabilities \eqn{P(t)} for the
  time interval \code{[t1, t2]}.  That is, the probabilities of
  occupying state \eqn{s} at time \eqn{t2}
  conditionally on occupying state \eqn{r} at time \eqn{t1}.
  Rows correspond to "from-state" and columns to "to-state". 
}
\seealso{
  \code{\link{pmatrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{}
