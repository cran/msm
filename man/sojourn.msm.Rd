\name{sojourn.msm}
\title{Estimate mean sojourn times from a multi-state model}
\alias{sojourn.msm}
\description{
  Auxiliary function called by \code{\link{msm}}. 
}
\usage{
sojourn.msm(Qmatrices, covmat, foundstderr)
}
\arguments{
  \item{Qmatrices}{List of estimated intensity and covariate matrices as
    given by \code{\link{msm}}.}
  \item{covmat}{Covariance matrix of the full set of estimates from the
    model.}
  \item{foundstderr}{(\code{TRUE} or \code{FALSE}) Was \code{\link{optim}} able to find standard errors? }
}
\value{
    A list with components:

    \code{mean} = estimated mean sojourn times in the transient states

    \code{se} = corresponding standard errors.
}
\details{
  The mean sojourn time in stage \eqn{r} is estimated by \eqn {1 /
    q_{rr}},  where \eqn{q_{r}} is the \eqn{r}th entry on the diagonal
  of the estimated transition intensity matrix.  Calls
  \code{\link{deltamethod}} to find standard errors.
}
\seealso{
  \code{\link{msm}}, \code{\link{deltamethod}}
}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{models}