\name{crudeinits.msm}
\alias{crudeinits.msm}
\title{Calculate crude initial values for transition intensities}
\description{
  Calculates crude initial values for transition intensities by assuming
  that the data represent the exact transition times of the Markov
  process. 
}
\usage{
crudeinits.msm(state, time, subject, qmatrix, data=NULL, fromto=FALSE,
               fromstate=NULL, tostate=NULL, timelag=NULL, check=FALSE)
}
\arguments{
  \item{state}{Observed states, assumed to be ordered by time
  within each subject.}
  \item{time}{Times at which the elements of \code{state} were observed.}
  \item{subject}{Subject identification numbers corresponding to
    \code{state}. If not given, all observations are assumed to be on
    the same subject.}

  \item{data}{An optional data frame in which the variables represented
    by \code{subject} and \code{state} can be found.}

  \item{qmatrix}{Matrix of indicators for the allowed
    transitions.  Specified as described in the help page for
    \code{\link{msm}}.  An initial value will be estimated for each
    value of qmatrix that is 1.}

  \item{fromto}{If \code{TRUE}, then the data are given as three vectors, 
    \emph{from-state, to-state, time-difference},
    representing the set of observed transitions between states, and the
    time taken by each one. Otherwise, the data are given by \code{formula},
    containing observation times, corresponding observed states, and a
    vector of corresponding subject identification numbers \code{subject}.
  }
  \item{fromstate}{ Starting states for the observed transitions
    (required if \code{fromto == TRUE} ). }

  \item{tostate}{ Finishing states for the observed transitions
    (required if \code{fromto == TRUE} ). }

  \item{timelag}{ Time difference between observing \code{fromstate}
    and \code{tostate}  (required if \code{fromto == TRUE} ). }

  \item{check}{ Should the supplied state, time and subject data
    be checked for inconsistency with the supplied \code{qmatrix}, for
    example, backward transitions in the data with a forward progression
    model.  This check is rather slow, and not performed by default. }
}
\value{
  A named vector of estimates for the initial values of each transition
  intensity.  For simple models with no covariates, this can be used as
  the \code{inits} argument to \code{\link{msm}}.
}
\details{

  If we observe \eqn{n_{rs}}{n_rs} transitions from state \eqn{r} to state \eqn{s}, and a
  total of \eqn{n_r} transitions from state \eqn{r}, then \eqn{q_{rs} /
    q_{rr}}{q_rs / q_rr} can be estimated by \eqn{n_{rs} / n_r}{n_rs / n_r}. Then, given a total of \eqn{T_r} years
  spent in state \eqn{r},  the mean sojourn time \eqn{1 / q_{rr}}{1 / q_rr} can be estimated as \eqn{T_r / n_r}.
  Thus, \eqn{n_{rs} / T_r}{n_rs / T_r} is a crude estimate of \eqn{q_{rs}}{q_rs}.

}
\seealso{
  \code{\link{statetable.msm}}
}
\examples{
## Aortic aneurysm growth data
data(aneur)

## Three-state one-way progression model 
oneway.q <- rbind( c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1), c(0, 0, 0, 0))

## Crude initial values for 1-2, 2-3 and 3-4 transition intensities 
crudeinits.msm(fromstate=from, tostate=to, timelag=dt, fromto=TRUE,
               data=aneur, qmatrix=oneway.q)
}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{}