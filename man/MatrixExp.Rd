\name{MatrixExp}
\title{Matrix exponential}
\alias{MatrixExp}
\description{
  Calculates the exponential of a square matrix.
}
\usage{
MatrixExp(mat, t = 1, n = 20, k = 3)
}
\arguments{
  \item{mat}{A square matrix}
  \item{t}{An optional scaling factor for the eigenvalues of \code{mat}}
  \item{n}{Number of terms in the series approximation to the
    exponential}
  \item{k}{Underflow correction factor}
}
\value{
  The exponentiated matrix \eqn{exp(mat)}.
}
\details{
  The exponential \eqn{E} of a square matrix \eqn{M} is calculated as 

  \deqn{E = U exp(D) U^{-1}}

  where \eqn{D} is a diagonal matrix with the eigenvalues
  of \eqn{M} on the diagonal, \eqn{exp(D)} is a diagonal matrix with the
  exponentiated eigenvalues
  of \eqn{M} on the diagonal, and \eqn{U} is a matrix whose columns are the
  eigenvectors of \eqn{M}.

  However, if \eqn{M} has repeated eigenvalues, then its eigenvector
  matrix is non-invertible. In this case, the matrix exponential is
  calculated using a power series approximation,

  \deqn{exp(M) = I + M + M^2/2 + M^3 / 3! + M^4 / 4! + ...}

  For a continuous-time homogeneous Markov process with transition
intensity matrix
  \eqn{Q},  the probability of occupying state \eqn{s} at time \eqn{u + t}
  conditional on occupying state \eqn{r} at time \eqn{u} is given by the
  \eqn{(r,s)} entry of the matrix \eqn{exp(tQ)}. 
  
  The series approximation method was adapted from the corresponding
  function in Jim Lindsey's R package \code{rmutil}.
}
\references{
  Cox, D. R. and Miller, H. D. \emph{The theory of stochastic
    processes}, Chapman and Hall, London (1965)
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{math}
