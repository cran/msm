## Estimate total length of stay in a given state.



#' Total length of stay, or expected number of visits
#' 
#' Estimate the expected total length of stay, or the expected number of
#' visits, in each state, for an individual in a given period of evolution of a
#' multi-state model.
#' 
#' The expected total length of stay in state \eqn{j} between times \eqn{t_1}
#' and \eqn{t_2}, from the point of view of an individual in state \eqn{i} at
#' time 0, is defined by the integral from \eqn{t_1} to \eqn{t_2} of the
#' \eqn{i,j} entry of the transition probability matrix \eqn{P(t) = Exp(tQ)},
#' where \eqn{Q} is the transition intensity matrix.
#' 
#' The corresponding expected number of visits to state \eqn{j} (excluding the
#' stay in the current state at time 0) is \eqn{\sum_{i!=j} T_i Q_{i,j}}, where
#' \eqn{T_i} is the expected amount of time spent in state \eqn{i}.
#' 
#' More generally, suppose that \eqn{\pi_0}{pi_0} is the vector of
#' probabilities of being in each state at time 0, supplied in \code{start},
#' and we want the vector \eqn{\mathbf{x}}{x} giving the expected lengths of
#' stay in each state.  The corresponding integral has the following solution
#' (van Loan 1978; van Rosmalen et al. 2013)
#' 
#' \deqn{\mathbf{x} =
#' \left[
#' \begin{array}{ll}
#' 1  &  \mathbf{0}_K
#' \end{array}
#' \right]
#' Exp(t Q')
#' \left[
#' \begin{array}{l} \mathbf{0}_K\\I_K
#' \end{array}
#' \right]
#' }{x = [1, 0_K] Exp(t Q') [0_K, I_K]'}
#'
#' where \deqn{Q' = \left[
#' \begin{array}{ll} 0 & \mathbf{\pi}_0\\
#' \mathbf{0}_K  &  Q - rI_K
#' \end{array}
#' \right]
#' }{Q' = rbind(c(0, pi_0), cbind(0_K, Q - r I_K))}
#' 
#' \eqn{\pi_0}{pi_0} is the row vector of initial state probabilities supplied
#' in \code{start}, \eqn{\mathbf{0}_K}{0_K} is the row vector of K zeros,
#' \eqn{r} is the discount rate, \eqn{I_K}{I_K} is the K x K identity matrix,
#' and \eqn{Exp} is the matrix exponential.
#' 
#' Alternatively, the integrals can be calculated numerically, using the
#' \code{\link{integrate}} function.  This may take a long time for models with
#' many states where \eqn{P(t)} is expensive to calculate.  This is required
#' where \code{tot = Inf}, since the package author is not aware of any
#' analytic expression for the limit of the above formula as \eqn{t} goes to
#' infinity.
#' 
#' With the argument \code{num.integ=TRUE}, numerical integration is used even
#' where the analytic solution is available. This facility is just provided for
#' checking results against versions 1.2.4 and earlier, and will be removed
#' eventually. Please let the package maintainer know if any results are
#' different.
#' 
#' For a model where the individual has only one place to go from each state,
#' and each state is visited only once, for example a progressive disease model
#' with no recovery or death, these are equal to the mean sojourn time in each
#' state.  However, consider a three-state health-disease-death model with
#' transitions from health to disease, health to death, and disease to death,
#' where everybody starts healthy.  In this case the mean sojourn time in the
#' disease state will be greater than the expected length of stay in the
#' disease state.  This is because the mean sojourn time in a state is
#' conditional on entering the state, whereas the expected total time diseased
#' is a forecast for a healthy individual, who may die before getting the
#' disease.
#' 
#' In the above formulae, \eqn{Q} is assumed to be constant over time, but the
#' results generalise easily to piecewise-constant intensities.  This function
#' automatically handles models fitted using the \code{pci} option to
#' \code{\link{msm}}. For any other inhomogeneous models, the user must specify
#' \code{piecewise.times} and \code{piecewise.covariates} arguments to
#' \code{\link{totlos.msm}}.
#' 
#' @aliases totlos.msm envisits.msm
#' @param x A fitted multi-state model, as returned by \code{\link{msm}}.
#' @param start Either a single number giving the state at the beginning of the
#' period, or a vector of probabilities of being in each state at this time.
#' @param end States to estimate the total length of stay (or number of visits)
#' in. Defaults to all states.  This is deprecated, since with the analytic
#' solution (see "Details") it doesn't save any computation to only estimate
#' for a subset of states.
#' @param fromt Time from which to estimate.  Defaults to 0, the beginning of
#' the process.
#' @param tot Time up to which the estimate is made.  Defaults to infinity,
#' giving the expected time spent in or number of visits to the state until
#' absorption. However, the calculation will be much more efficient if a finite
#' (potentially large) time is specified: see the "Details" section.  For
#' models without an absorbing state, \code{t} must be specified.
#' @param covariates The covariate values to estimate for.  This can either
#' be:\cr
#' 
#' the string \code{"mean"}, denoting the means of the covariates in the data
#' (this is the default),\cr
#' 
#' the number \code{0}, indicating that all the covariates should be set to
#' zero,\cr
#' 
#' or a list of values, with optional names. For example
#' 
#' \code{list (60, 1)}
#' 
#' where the order of the list follows the order of the covariates originally
#' given in the model formula, or a named list,
#' 
#' \code{list (age = 60, sex = 1)}
#' 
#' @param piecewise.times Times at which piecewise-constant intensities change.
#' See \code{\link{pmatrix.piecewise.msm}} for how to specify this. This is
#' only required for time-inhomogeneous models specified using explicit
#' time-dependent covariates, and should not be used for models specified using
#' "pci".
#' @param piecewise.covariates Covariates on which the piecewise-constant
#' intensities depend. See \code{\link{pmatrix.piecewise.msm}} for how to
#' specify this.
#' @param num.integ Use numerical integration instead of analytic solution (see
#' below).
#' @param discount Discount rate in continuous time.
#' @param env Supplied to \code{\link{totlos.msm}}.  If \code{TRUE}, return the
#' expected number of visits to each state. If \code{FALSE}, return the total
#' length of stay in each state. \code{\link{envisits.msm}} simply calls
#' \code{\link{totlos.msm}} with \code{env=TRUE}.
#' @param ci If \code{"normal"}, then calculate a confidence interval by
#' simulating \code{B} random vectors from the asymptotic multivariate normal
#' distribution implied by the maximum likelihood estimates (and covariance
#' matrix) of the log transition intensities and covariate effects, then
#' calculating the total length of stay for each replicate.
#' 
#' If \code{"bootstrap"} then calculate a confidence interval by non-parametric
#' bootstrap refitting.  This is 1-2 orders of magnitude slower than the
#' \code{"normal"} method, but is expected to be more accurate. See
#' \code{\link{boot.msm}} for more details of bootstrapping in \pkg{msm}.
#' 
#' If \code{"none"} (the default) then no confidence interval is calculated.
#' @param cl Width of the symmetric confidence interval, relative to 1
#' @param B Number of bootstrap replicates
#' @param cores Number of cores to use for bootstrapping using parallel
#' processing. See \code{\link{boot.msm}} for more details.
#' @param ... Further arguments to be passed to the \code{\link{integrate}}
#' function to control the numerical integration.
#' @return A vector of expected total lengths of stay
#' (\code{\link{totlos.msm}}), or expected number of visits
#' (\code{\link{envisits.msm}}), for each transient state.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{sojourn.msm}}, \code{\link{pmatrix.msm}},
#' \code{\link{integrate}}, \code{\link{boot.msm}}.
#' @references C. van Loan (1978). Computing integrals involving the matrix
#' exponential. IEEE Transactions on Automatic Control 23(3)395-404.
#' 
#' J. van Rosmalen, M. Toy and J.F. O'Mahony (2013). A mathematical approach
#' for evaluating Markov models in continuous time without discrete-event
#' simulation.  Medical Decision Making 33:767-779.
#' @keywords models
#' @export totlos.msm
totlos.msm <- function(x, start=1, end=NULL, fromt=0, tot=Inf, covariates="mean",
                       piecewise.times=NULL,
                       piecewise.covariates=NULL,
                       num.integ=FALSE, discount=0,
                       env=FALSE,
                       ci=c("none","normal","bootstrap"), # calculate a confidence interval
                       cl = 0.95, # width of symmetric confidence interval
                       B = 1000, # number of bootstrap replicates
                       cores=NULL,
                       ...)
{
  if (!inherits(x, "msm")) stop("expected x to be a msm model")
  nst <- x$qmodel$nstates
  if (!is.numeric(start) ||
      ((length(start)==1) && (! start %in% 1 : nst)))
    stop("start should be a state in 1, ..., ", nst, " or a vector of length ",nst)
  else if (length(start) == 1) {p0 <- rep(0, nst); p0[start] <- 1; start <- p0}
  else if (length(start) > 1) {
    if (length(start) != nst)
      stop("start should be a state in 1, ..., ", nst, " or a vector of length ",nst)
  }
  if (is.null(end)) end <- 1 : nst
  if (! all(end %in% 1 : nst)) stop("end should be a set of states in 1, ..., ", nst)
  if (!is.numeric(fromt) || !is.numeric(tot) || length(fromt) != 1 || length(tot) != 1 || fromt < 0 || tot < 0)
    stop("fromt and tot must be single non-negative numbers")
  if (fromt > tot) stop("tot must be greater than fromt")
  if (length(absorbing.msm(x)) == 0)
    if (tot==Inf) stop("Must specify a finite end time for a model with no absorbing state")

  if (!is.null(x$pci)){
    piecewise.times <- x$pci
    piecewise.covariates <- msm.fill.pci.covs(x, covariates)
  }
  ncuts <- length(piecewise.times)
  npieces <- length(piecewise.covariates)
  if (!is.null(piecewise.times) && (!is.numeric(piecewise.times) || is.unsorted(piecewise.times)))
    stop("piecewise.times should be a vector of numbers in increasing order")
  if (!is.null(piecewise.covariates) && (npieces != ncuts + 1))
    stop("Number of piecewise.covariate lists must be one greater than the number of cut points")
  if (is.null(piecewise.covariates)) {
    ## define homogeneous model as piecewise with one piece
    npieces <- 1
    covs <- list(covariates)
    ptimes <- c(fromt, tot)
  } else {
    ## ignore all cut points outside [fromt,tot]
    keep <- which((piecewise.times > fromt) & (piecewise.times < tot))
    ## cov value between fromt and min(first cut, tot)
    cov1 <- piecewise.covariates[findInterval(fromt, piecewise.times) + 1]
    covs <- c(cov1, piecewise.covariates[keep+1])
    npieces <- length(covs)
    ptimes <- c(fromt, piecewise.times[keep], tot)
  }

  tmat <- envmat <- matrix(nrow=npieces, ncol=nst)
  if (tot==Inf) {
    tmat[,absorbing.msm(x)] <- Inf # set by hand or else integrate() will fail
    envmat[,absorbing.msm(x)] <- 1
    rem <- setdiff(seq_len(nst), absorbing.msm(x))
  }
  else rem <- seq_len(nst)
  for (i in 1:npieces) {
    from.t <- ptimes[i]
    to.t <- ptimes[i+1]
    Q <- qmatrix.msm(x, covariates=covs[[i]], ci="none")
    if (num.integ || to.t==Inf){
      for (j in rem){
        f <- function(time) {
          y <- numeric(length(time))
          for (k in seq_along(y))
            y[k] <- (start %*% pmatrix.msm(x, time[k], t1=0, covariates=covs[[i]], ci="none")) [j]
          y
        }
        tmat[i,j] <- integrate(f, from.t, to.t, ...)$value
      }
    } else {
      QQ <- rbind(c(0, start),
                  cbind(rep(0,nst), Q - discount*diag(nst)))
      tmat[i,] <- as.vector(c(1, rep(0, nst)) %*%
                            (MatrixExp(to.t*QQ) - MatrixExp(from.t*QQ)) %*%
                            rbind(rep(0, nst), diag(nst)))
    }
    Q0 <- Q; diag(Q0) <- 0
    envmat[i,rem] <- tmat[i,rem] %*% Q0[rem,rem]
  }
  res <- if (env) colSums(envmat) else colSums(tmat)
  names(res) <- rownames(x$qmodel$qmatrix)
  ci <- match.arg(ci)
  t.ci <- switch(ci,
                 bootstrap = totlos.ci.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates,
                                           piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates,
                                           discount=discount, env=env, cl=cl, B=B, cores=cores, ...),
                 normal = totlos.normci.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates,
                                            piecewise.times=piecewise.times, piecewise.covariates=piecewise.covariates,
                                            discount=discount, env=env, cl=cl, B=B, ...),
                 none = NULL)
  res <- if (ci=="none") res[end] else rbind(res, t.ci)[,end]
  class(res) <- c("msm.estbystate", class(res))
  res
}

#' @export
print.msm.estbystate <- function(x, ...){
  print(unclass(x))
}

#' @export
as.data.frame.msm.estbystate <- function(x,...){
  as.data.frame(unclass(x))
}

## Expected number of visits

#' @rdname totlos.msm
#' @export
envisits.msm <- function(x=NULL, start=1, end=NULL, fromt=0, tot=Inf, covariates="mean",
                         piecewise.times=NULL,  piecewise.covariates=NULL,
                         num.integ=FALSE, discount=0,
                         ci=c("none","normal","bootstrap"), # calculate a confidence interval
                         cl = 0.95, # width of symmetric confidence interval
                         B = 1000, # number of bootstrap replicates
                         cores=NULL,
                         ...)
{
  totlos.msm(x=x, start=start, end=end, fromt=fromt, tot=tot, covariates=covariates,
             piecewise.times=piecewise.times,
             piecewise.covariates=piecewise.covariates, num.integ=num.integ,
             discount=discount, env=TRUE, ci=ci, cl=cl, B=B, cores=cores, ...)
}
