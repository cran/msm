\name{prevalencemisc.msm}
\title{Calculate tables of observed and expected prevalences at a series
of times}
\alias{prevalencemisc.msm}
\description{
  This function is called by \code{\link{summary.msm}}, and provides a
  rough indication of the goodness of fit of a multi-state model with
  misclassification.
}
\usage{
prevalencemisc.msm(msm, times, observed=TRUE, expected=TRUE)
}
\arguments{
  \item{msm}{A fitted multi-state model with misclassification, produced
    by \code{\link{msm}}.}
  \item{times}{Series of times. The observed and
    expected prevalences of states will be estimated for the intervals
    between these times.}
  \item{observed}{Should the table of observed counts be produced - 
    defaults to \code{TRUE}.}
  \item{expected}{Should the table of expected counts be produced - 
    defaults to \code{TRUE}.}
}
\value{
  A list with components:
  
  \code{Observed} = table of observed numbers of individuals in each state in
  each time interval

  \code{Observed percentages} = corresponding percentage of the
  individuals at risk in each time interval  

  \code{Expected} = table of corresponding expected numbers

  \code{Expected percentages} = corresponding percentage of the
  individuals at risk in each time interval
  
}

\details{
  For models with misclassification only. This technique is adapted from
  that given by Satten and Longini (1996).

  The observed numbers of individuals in each underlying state in each
  time interval are estimated from the observed states as follows. If an
  individual has one observation in state \eqn{i} in interval \eqn{dt},
  then he/she contributes 1 to the total in underlying state \eqn{i} at
  that time. If an individual has observations in more than one state, for
  example three in state 1 and two in state 2, then the contributions are
  weighted. In this case the individual would contribute 0.6 to the state
  1 total and 0.4 to the state 2 total.

  The expected numbers of individuals in each state in each time
  interval are estimated using one-step-ahead prediction. For each
  observation time, the probability that this lies in each state,
  conditionally on the history of observations preceding that time, is
  computed as a ratio of the partial likelihoods. For each subject, the
  probabilities are averaged within each interval to obtain their
  contribution to the expected number in that state in that time
  interval. For example, if an individual has two observations in an
  interval, and corresponding forecast probabilities of 0.2, and 0.4 for
  state 1, then that individual would contribute 0.3 to the number
  of individuals expected in state 1 in that interval.

  If death times are known exactly (\code{death==TRUE} in the call to
  \code{\link{msm}}), then the observed number of deaths in an interval
  is known. The expected number of deaths in an interval \eqn{t_1, t_2}
  is calculated conditionally on the history before the interval. Let
  \eqn{p_i} be the probability that the \eqn{i}th individual dies in
  \eqn{t_1, t_2}, conditionally on all their observations before
  \eqn{t_1}. Then the expected number of deaths is the sum over
  individuals \eqn{i} of \eqn{p_i}. 

}

\references{
  Satten, G.A. and Longini, I.M.  Markov chains with measurement error:
  estimating the 'true' course of a marker of the progression of human
  immunodeficiency virus disease (with discussion) \emph{Applied
    Statistics} 45(3): 275-309 (1996)
}

\seealso{
  \code{\link{msm}}, \code{\link{summary.msm}}
}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{models}