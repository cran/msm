\name{simmulti.msm}
\title{Simulate multiple trajectories from a multi-state Markov model with arbitrary observation times}
\alias{simmulti.msm}
\description{Simulate a number of individual realisations from a
  multi-state Markov process. Observations of the process are made at
  specified arbitrary times for each individual. 
}
\usage{
simmulti.msm(data, qmatrix, beta, death = FALSE,  tunit = 1)
}

\arguments{

  \item{data}{A data frame with mandatory columns named \code{subject},
    corresponding to subject identification numbers, and \code{time},
    representing observation times. Other named columns of the data frame
    represent covariates.
  }
  
  \item{qmatrix}{
    The transition intensity matrix of the
    Markov process. The diagonal of \code{qmatrix} is ignored,
    and computed as appropriate so that the rows sum to zero. For
    example, a possible \code{qmatrix} for a three state illness-death
    model with recovery is:

    \code{rbind(
    c( 0, 0.1, 0.02 ),
    c( 0.1, 0, 0.01 ),
    c( 0, 0, 0 )
    )}

  }

  \item{beta}{Matrix of linear covariate effects on log transition
    intensities. The rows correspond to different covariates, and the
    columns to the transition intensities. The intensities are ordered
    by reading across rows of the intensity matrix, starting with the
    first, counting the positive off-diagonal elements of the matrix.}

  \item{death}{Vector of indices of the death states.  A death state is
    an absorbing state whose time of entry is known exactly, but the
    individual is assumed to be in an unknown transient state ("alive")
    at the previous instant.  This is the usual situation for times of death in
    chronic disease monitoring data.  For example, if you specify
    \code{death = c(4, 5)} then states 4 and 5 are assumed to be death
    states.

    \code{death = TRUE} indicates that the
    final state is a death state, and \code{death = FALSE} (the default)
    indicates that there is no death state. }

  \item{tunit}{No longer used, from msm version 0.3.2. Death times are
    now assumed to be exact, rather than accurate within one day. 
  }

}

\value{A data frame with columns,
  \item{subject}{Subject identification indicators}
  \item{time}{Observation times}
  \item{state}{Simulated state at the corresponding time}
  plus any supplied covariates. 
}
\details{  
  \code{\link{sim.msm}} is called repeatedly to produce a simulated
  trajectory for each individual. The state at each specified
  observation time is then taken to produce a new column \code{state}.
  The effect of time-dependent covariates on the transition intensity
  matrix for an individual is determined by assuming that the covariate is a step function
  which remains constant in between the individual's observation times.
  If the subject enters an absorbing state, then only the first
  observation in that state is kept in the data frame. Rows corresponding to future
  observations are deleted.  The entry times into states given in
  \code{death} are assumed to be known exactly}.
}
\seealso{ \code{\link{sim.msm}}
}
\examples{
### Simulate 100 individuals with common observation times
sim.df <- data.frame(subject = rep(1:100, rep(13,100)), time = rep(seq(0, 24, 2), 100))
qmatrix <- rbind(c(-0.11,   0.1,  0.01 ),
                 c(0.05,   -0.15,  0.1 ),
                 c(0.02,   0.07, -0.09))
simmulti.msm(sim.df, qmatrix)
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
