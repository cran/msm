\name{statetable.msm}
\alias{statetable.msm}
\title{Table of transitions}
\description{
  Calculates a frequency table counting the number of times each pair of
  states were observed in successive observation times.  
}
\usage{
statetable.msm(state, subject, data=NULL)
}
\arguments{
  \item{state}{Observed states, assumed to be ordered by time
  within each subject.}
  \item{subject}{Subject identification numbers corresponding to
    \code{state}. If not given, all observations are assumed to be on
    the same subject.}
  \item{data}{An optional data frame in which the variables represented
    by \code{subject} and \code{state} can be found.}
}
\value{
  A frequency table with starting states as rows and finishing states as
  columns.
}
\details{
  If the data are in "from state, to state, time lag" format, then this table
  can simply be produced by tabulating "from state" against "to state":

  \code{table(fromstate, tostate)}
}
\seealso{
  \code{\link{crudeinits.msm}}
}
\examples{
## Heart transplant data
data(heart)

## 148 deaths from state 1, 48 from state 2 and 55 from state 3.
statetable.msm(state, PTNUM, data=heart)

}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{}