\name{boot.msm}
\alias{boot.msm}
\title{Bootstrap resampling for multi-state models}
\description{
  Draw a number of bootstrap resamples, refit a msm model to the
  resamples, and calculate statistics on the refitted models.
}
\usage{
boot.msm(x, stat=pmatrix.msm, B=500, file=NULL)
}
\arguments{
  \item{x}{A fitted msm model, as output by \code{\link{msm}}.}
  \item{stat}{A function to call on each refitted msm model. By default
    this is \code{\link{pmatrix.msm}}, returning the transition
    probability matrix in one time unit. If
    \code{NULL} then no function is computed.}
  \item{B}{Number of bootstrap resamples.}
  \item{file}{Name of a file in which to save partial results after each
    replicate. This is saved using \code{\link{save}} and can be
    restored using \code{\link{load}}, producing an object called
    \code{boot.list} containing the partial results.}  
}
\details{
  The bootstrap datasets are computed by resampling independent
  transitions between pairs of states (for non-hidden models without
  censoring), or independent patient series (for hidden models or
  models with censoring).

  Confidence intervals or standard errors for the corresponding
  statistic can be calculated by summarising the returned list of
  \code{B} replicated outputs.  This is currently implemented for the
  transition probability matrix and total length of stay, see
  \code{\link{pmatrix.msm}}, \code{\link{totlos.msm}}.  At the moment,
  for other outputs, users will have to write their own code to
  summarise the output of \code{\link{boot.msm}}.

  Most of \pkg{msm}'s output functions present confidence intervals
  based on asymptotic standard errors calculated from the
  Hessian. These are expected to be underestimates of the true standard
  errors (Cramer-Rao lower bound).  Bootstrapping may give a more
  accurate estimate of the uncertainty.   
  
  All objects used in the original call to \code{\link{msm}} which
  produced \code{x}, such as the \code{qmatrix}, should be in the
  working environment, or else \code{boot.msm} will produce an
  \dQuote{object not found} error.  This enables \code{boot.msm} to
  refit the original model to the replicate datasets. 
  
  If \code{stat} is \code{NULL}, then \code{B} different \code{msm} model
  objects will be stored in memory. This is unadvisable, as \code{msm}
  objects tend to be large, as they contain the original data used for
  the \code{msm} fit, so this will be wasteful of memory.

  To specify more than one statistic, write a function consisting of a
  list of different function calls, for example,

  \code{stat = function(x) list (pmatrix.msm(x, t=1), pmatrix.msm(x,
    t=2))}
}
\value{
  A list with \code{B} components, containing the result of calling
  function \code{stat} on each of the refitted models.  If \code{stat}
  is \code{NULL}, then each component just contains the refitted
  model.  If one of the \code{B} model fits was unsuccessful and
  resulted in an error, then the corresponding list component will
  contain the error message. 
}
\references{
  Efron, B. and Tibshirani, R.J. (1993) \emph{An Introduction to the
    Bootstrap}, Chapman and Hall.
}
\seealso{
  \code{\link{pmatrix.msm}}, \code{\link{totlos.msm}}
}
\examples{
\dontrun{
  ## Psoriatic arthritis example 
  data(psor)
  psor.q <- rbind(c(0,0.1,0,0),c(0,0,0.1,0),c(0,0,0,0.1),c(0,0,0,0))
  psor.msm <- msm(state ~ months, subject=ptnum, data=psor, qmatrix = psor.q, covariates = ~ollwsdrt+hieffusn, constraint = list(hieffusn=c(1,1,1),ollwsdrt=c(1,1,2)), control = list(REPORT=1,trace=2), method="BFGS")
  ## Bootstrap the baseline transition intensity matrix.  This will take a long time. 
  q.list <- boot.msm(psor.msm, function(x)x$Qmatrices$baseline)
  ## Manipulate the resulting list of matrices to calculate bootstrap standard errors. 
  apply(array(unlist(q.list), dim=c(4,4,5)), c(1,2), sd)
  ## Similarly calculate a bootstrap 95\% confidence interval
  apply(array(unlist(q.list), dim=c(4,4,5)), c(1,2), function(x)quantile(x, c(0.025, 0.975)))
  ## Bootstrap standard errors are larger than the asymptotic standard errors calculated from the Hessian
  psor.msm$QmatricesSE$baseline
}
}
\author{C.H.Jackson <chris.jackson@imperial.ac.uk>}
\keyword{}