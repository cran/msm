\name{viterbi.msm}
\title{Calculate the most likely path through underlying stages}
\alias{viterbi.msm}

\description{
  For a fitted multi state model with misclassification, the Viterbi
  algorithm recursively constructs the path through the underlying stages with the
  highest probability.
}

\usage{
viterbi.msm(x)
}
\arguments{
  \item{x}{A fitted multi-state model with misclassification, as
  produced by \code{\link{msm}}}
}
\value{
    A data frame with columns:

    \code{subject} = subject identification numbers

    \code{time} = times of observations

    \code{observed} = corresponding observed states

    \code{fitted} = corresponding fitted states found by Viterbi
    recursion.
}
\seealso{
  \code{\link{msm}}
}
\references{
  Durbin, R., Eddy, S., Krogh, A. and Mitchison, G. \emph{Biological
    sequence analysis}, Cambridge University Press, 1998.
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
