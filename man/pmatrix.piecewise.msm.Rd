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
  time interval \code{[t1, tn]}.  That is, the probabilities of
  occupying state \eqn{s} at time \eqn{tn}
  conditionally on occupying state \eqn{r} at time \eqn{t1}.
  Rows correspond to "from-state" and columns to "to-state". 
}
\details {
  Suppose a multi-state model, \code{x.msm}, has been fitted, in which
  the transition intensity matrix \eqn{Q(x(t))} is modelled in terms of
  time-dependent covariates \eqn{x(t)}.  The transition probability
  matrix \eqn{P(t_1, t_n)}{P(t1, tn)} for the time interval \eqn{(t_1, t_n)}{(t1, tn)} cannot be
  calculated from the estimated intensity matrix as
  \eqn{exp((t_n - t_1) Q}{exp((tn - t1) Q}, because \eqn{Q} varies within the interval
  \eqn{t_1, t_n}{t1, tn}.  However, if the covariates are piecewise-constant, 
  or can be approximated as piecewise-constant, then we can calculate
  \eqn{P(t_1, t_n)}{P(t1, tn)} by multiplying together individual matrices
  \eqn{P(t_i,
  t_{i+1}) = exp((t_{i+1} - t_i) Q}, calculated over intervals where Q is
  constant:

  \deqn {P(t_1, t_n) =  P(t_1, t_2) P(t_2, t_3) \dots P(t_{n-1}, t_n)} {P(t1, tn) =  P(t1, t2) P(t2, t3) \dots P(tn-1, tn)}.
  
}
\seealso{
  \code{\link{pmatrix.msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{}
