### FUNCTIONS FOR SIMULATING FROM MULTI-STATE MODELS

### from help(sample) in base R
resample <- function(x, size, ...)
    if(length(x) <= 1) {
        if(!missing(size) && size == 0) x[FALSE] else x
    } else sample(x, size, ...)

### General function to simulate one individual's realisation from a continuous-time Markov model
### Produces the exact times of transition



#' Simulate one individual trajectory from a continuous-time Markov model
#' 
#' Simulate one realisation from a continuous-time Markov process up to a given
#' time.
#' 
#' The effect of time-dependent covariates on the transition intensity matrix
#' for an individual is determined by assuming that the covariate is a step
#' function which remains constant in between the individual's observation
#' times.
#' 
#' @param qmatrix The transition intensity matrix of the Markov process. The
#' diagonal of \code{qmatrix} is ignored, and computed as appropriate so that
#' the rows sum to zero. For example, a possible \code{qmatrix} for a three
#' state illness-death model with recovery is:
#' 
#' \code{rbind( c( 0, 0.1, 0.02 ), c( 0.1, 0, 0.01 ), c( 0, 0, 0 ) )}
#' 
#' @param maxtime Maximum time for the simulated process.
#' @param covs Matrix of time-dependent covariates, with one row for each
#' observation time and one column for each covariate.
#' @param beta Matrix of linear covariate effects on log transition
#' intensities. The rows correspond to different covariates, and the columns to
#' the transition intensities. The intensities are ordered by reading across
#' rows of the intensity matrix, starting with the first, counting the positive
#' off-diagonal elements of the matrix.
#' @param obstimes Vector of times at which the covariates are observed.
#' @param start Starting state of the process. Defaults to 1.
#' @param mintime Starting time of the process. Defaults to 0.
#' @return A list with components,
#' 
#' \item{states}{Simulated states through which the process moves.  This ends
#' with either an absorption before \code{obstime}, or a transient state at
#' \code{obstime}. }
#' 
#' \item{times}{Exact times at which the process changes to the corresponding
#' states}
#' 
#' \item{qmatrix}{The given transition intensity matrix}
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{simmulti.msm}}
#' @keywords datagen
#' @examples
#' 
#' 
#' qmatrix <- rbind(
#'                  c(-0.2,   0.1,  0.1 ),
#'                  c(0.5,   -0.6,  0.1 ),
#'                  c(0,  0,  0)
#'                  )
#' sim.msm(qmatrix, 30)
#' 
#' @export sim.msm
sim.msm <- function(qmatrix,   # intensity matrix
                    maxtime,   # maximum time for realisations
                    covs=NULL,     # covariate matrix, nobs rows, ncovs cols
                    beta=NULL,     # matrix of cov effects on qmatrix. ncovs rows, nintens cols.
                    obstimes=0, # times at which time-dependent covariates change
                    start = 1,     # starting state
                    mintime = 0    # time to start from
                    )
{
    ## Keep only times where time-dependent covariates change
    if (!is.null(covs)) {
        covs2 <- collapse.covs(covs)
        covs <- covs2$covs
        obstimes <- obstimes[covs2$ind]
    }
    else {obstimes <- mintime; covs <- beta <- 0}
    if (is.vector(beta)) beta <- matrix(beta, ncol=length(beta))
    nct <- length(obstimes)
    nstates <- nrow(qmatrix)
    ## Form an array of qmatrices, one for each covariate change-time
    qmatrices <- array(rep(t(qmatrix), nct), dim=c(dim(qmatrix), nct))
    qmatrices[rep(t(qmatrix)>0, nct)] <- qmatrices[rep(t(qmatrix)>0, nct)] * as.numeric(exp(t(beta)%*%t(covs))) # nintens*nobs
    for (i in 1:nct)
        qmatrices[,,i] <- msm.fixdiag.qmatrix(t(qmatrices[,,i]))
    cur.t <- mintime; cur.st <- next.st <- start; rem.times <- obstimes; t.ind <- 1
    nsim <- 0; max.nsim <- 10
    simstates <- simtimes <- numeric(max.nsim) ## allocate memory up-front for simulated outcome
    absorb <- absorbing.msm(qmatrix=qmatrix)
    ## Simulate up to maxtime or absorption
    while (cur.t < maxtime) {
        nsim <- nsim + 1
        cur.st <- next.st
        simstates[nsim] <- cur.st; simtimes[nsim] <- cur.t
        if (cur.st %in% absorb) break;
        rate <- -qmatrices[cur.st,cur.st, t.ind:length(obstimes)]
        nextlag <- rpexp(1, rate, rem.times-rem.times[1])
        cur.t <- cur.t + nextlag
        t.ind <- which.min((cur.t - obstimes)[cur.t - obstimes > 0])
        rem.times <- cur.t
        if (any(obstimes > cur.t))
            rem.times <- c(rem.times, obstimes[(t.ind+1): length(obstimes)])
        cur.q <- qmatrices[,, t.ind]
        next.st <- resample((1:nstates)[-cur.st], size=1, prob = cur.q[cur.st, -cur.st])
        if (nsim > max.nsim) { ## need more memory for simulated outcome, allocate twice as much
            simstates <- c(simstates, numeric(max.nsim))
            simtimes <- c(simtimes, numeric(max.nsim))
            max.nsim <- max.nsim*2
        }
    }
    ## If process hasn't absorbed by the end, then include a censoring time
    if (cur.t >= maxtime) {
        nsim <- nsim+1
        simstates[nsim] <- cur.st
        simtimes[nsim] <- maxtime
    }
    list(states = simstates[1:nsim], times = simtimes[1:nsim], qmatrix = qmatrix)
}


## Drop rows of a covariate matrix which are identical to previous row
## Similar method to R's unique.data.frame

collapse.covs <- function(covs)
{
    if (nrow(covs)==1) list(covs=covs, ind=1)
    else {
        pcovs <- apply(covs, 1, function(x) paste(x, collapse="\r"))
        lpcovs <- c("\r", pcovs[1:(length(pcovs)-1)])
        ind <- pcovs!=lpcovs
        list(covs=covs[ind,,drop=FALSE], ind=which(ind))
    }
}

### Given a simulated Markov model, get the current state at various observation times
### By default, only keep one observation in the absorbing state

getobs.msm <- function(sim, obstimes, death=FALSE, drop.absorb=TRUE)
{
    absorb <- absorbing.msm(qmatrix=sim$qmatrix)
                                        # Only keep one observation in the absorbing state
    if (drop.absorb && any(sim$states %in% absorb)) {
        if (any(sim$states %in% death))
            keep <- which(obstimes < max(sim$times))
        else {
            lo <- c(-Inf, obstimes[1:(length(obstimes)-1)])
            keep <- which(lo <= max(sim$times))
        }
    }
    else keep <- 1 : length(obstimes)
    obstimes <- obstimes[keep]
    state <- sim$states[rowSums(outer(obstimes, sim$times, ">="))]
    time <- obstimes
    if (any(sim$states %in% death)) { # Keep the exact death time if required
        state <- c(state, sim$states[sim$states %in% death])
        time <- c(time, sim$times[sim$states %in% death])
        state <- state[order(time)]
        time <- time[order(time)]
        keep <- c(keep, max(keep)+1)
    }
    list(state = state, time = time, keep=keep)
}

### Simulate a multi-state Markov or hidden Markov model dataset using fixed observation times

### Would it be better to make specification of covariate model consistent with model fitting function?
### e.g. separate hcovariates and  covariates formulae,
### plus covinits and hcovinits?



#' Simulate multiple trajectories from a multi-state Markov model with
#' arbitrary observation times
#' 
#' Simulate a number of individual realisations from a continuous-time Markov
#' process. Observations of the process are made at specified arbitrary times
#' for each individual, giving panel-observed data.
#' 
#' \code{\link{sim.msm}} is called repeatedly to produce a simulated trajectory
#' for each individual. The state at each specified observation time is then
#' taken to produce a new column \code{state}.  The effect of time-dependent
#' covariates on the transition intensity matrix for an individual is
#' determined by assuming that the covariate is a step function which remains
#' constant in between the individual's observation times.  If the subject
#' enters an absorbing state, then only the first observation in that state is
#' kept in the data frame. Rows corresponding to future observations are
#' deleted.  The entry times into states given in \code{death} are assumed to
#' be known exactly.
#' 
#' @param data A data frame with a mandatory column named \code{time},
#' representing observation times.  The optional column named \code{subject},
#' corresponds to subject identification numbers. If not given, all
#' observations are assumed to be on the same individual. Observation times
#' should be sorted within individuals.  The optional column named \code{cens}
#' indicates the times at which simulated states should be censored. If
#' \code{cens==0} then the state is not censored, and if \code{cens==k}, say,
#' then all simulated states at that time which are in the set
#' \code{censor.states} are replaced by \code{k}.  Other named columns of the
#' data frame represent any covariates, which may be time-constant or
#' time-dependent.  Time-dependent covariates are assumed to be constant
#' between the observation times.
#' @param qmatrix The transition intensity matrix of the Markov process, with
#' any covariates set to zero.  The diagonal of \code{qmatrix} is ignored, and
#' computed as appropriate so that the rows sum to zero. For example, a
#' possible \code{qmatrix} for a three state illness-death model with recovery
#' is:
#' 
#' \code{rbind( c( 0, 0.1, 0.02 ), c( 0.1, 0, 0.01 ), c( 0, 0, 0 ) )}
#' @param covariates List of linear covariate effects on log transition
#' intensities. Each element is a vector of the effects of one covariate on all
#' the transition intensities. The intensities are ordered by reading across
#' rows of the intensity matrix, starting with the first, counting the positive
#' off-diagonal elements of the matrix.
#' 
#' For example, for a multi-state model with three transition intensities, and
#' two covariates \code{x} and \code{y} on each intensity,
#' 
#' \code{covariates=list(x = c(-0.3,-0.3,-0.3), y=c(0.1, 0.1, 0.1))}
#' @param death Vector of indices of the death states.  A death state is an
#' absorbing state whose time of entry is known exactly, but the individual is
#' assumed to be in an unknown transient state ("alive") at the previous
#' instant.  This is the usual situation for times of death in chronic disease
#' monitoring data.  For example, if you specify \code{death = c(4, 5)} then
#' states 4 and 5 are assumed to be death states.
#' 
#' \code{death = TRUE} indicates that the final state is a death state, and
#' \code{death = FALSE} (the default) indicates that there is no death state.
#' @param start A vector with the same number of elements as there are distinct
#' subjects in the data, giving the states in which each corresponding
#' individual begins.  Or a single number, if all of these are the same.
#' Defaults to state 1 for each subject.
#' @param ematrix An optional misclassification matrix for generating observed
#' states conditionally on the simulated true states. As defined in
#' \code{\link{msm}}.
#' @param misccovariates Covariate effects on misclassification probabilities
#' via multinomial logistic regression. Linear effects operate on the log of
#' each probability relative to the probability of classification in the
#' correct state. In same format as \code{covariates}.
#' @param hmodel An optional hidden Markov model for generating observed
#' outcomes conditionally on the simulated true states. As defined in
#' \code{\link{msm}}.  Multivariate outcomes (\code{hmmMV}) are not supported.
#' @param hcovariates List of the same length as \code{hmodel}, defining any
#' covariates governing the hidden Markov outcome models.  Unlike in the
#' \code{msm} function, this should also define the values of the covariate
#' effects. Each element of the list is a named vector of the initial values
#' for each set of covariates for that state.  For example, for a three-state
#' hidden Markov model with two, one and no covariates on the state 1, 2 and 3
#' outcome models respectively,
#' 
#' \code{ hcovariates = list (c(acute=-8, age=0), c(acute=-8), NULL) }
#' @param censor.states Set of simulated states which should be replaced by a
#' censoring indicator at censoring times.  By default this is all transient
#' states (representing alive, with unknown state).
#' @param drop.absorb Drop repeated observations in the absorbing state,
#' retaining only one.
#' @return A data frame with columns, \item{subject}{Subject identification
#' indicators} \item{time}{Observation times} \item{state}{Simulated (true)
#' state at the corresponding time} \item{obs}{Observed outcome at the
#' corresponding time, if \code{ematrix} or \code{hmodel} was supplied}
#' \item{keep}{Row numbers of the original data.  Useful when
#' \code{drop.absorb=TRUE}, to show which rows were not dropped} plus any
#' supplied covariates.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{sim.msm}}
#' @keywords datagen
#' @examples
#' 
#' ### Simulate 100 individuals with common observation times
#' sim.df <- data.frame(subject = rep(1:100, rep(13,100)), time = rep(seq(0, 24, 2), 100))
#' qmatrix <- rbind(c(-0.11,   0.1,  0.01 ),
#'                  c(0.05,   -0.15,  0.1 ),
#'                  c(0.02,   0.07, -0.09))
#' simmulti.msm(sim.df, qmatrix)
#' 
#' @export simmulti.msm
simmulti.msm <- function(data,           # data frame with subject, times, covariates...
                         qmatrix,        # intensity matrix
                         covariates=NULL,  # initial values
                         death = FALSE,  # vector of indicators for "death" states, ie absorbing states whose entry time is known exactly,
                                        # but with unknown transient state at previous instant
                         start,         # starting states of the process, defaults to all 1.
                         ematrix = NULL,# misclassification matrix
                         misccovariates = NULL, # covariates on misclassification probabilities
                         hmodel = NULL,  # hidden Markov model formula
                         hcovariates = NULL,   # covariate effects on hidden Markov model response distribution
                         censor.states = NULL,
                         drop.absorb = TRUE
                         )
{

### Check consistency of qmatrix and covariate inits
    nstates <- nrow(qmatrix)
    msm.check.qmatrix(qmatrix)
    qmatrix <- msm.fixdiag.qmatrix(qmatrix)

### Subject, time and state
    if (!("subject" %in% names(data)))
        data$subject <- rep(1, nrow(data))
    if (!("time" %in% names(data)))
        stop("\"time\" column missing from data")
    data <- as.data.frame(data)
    subject <- data[,"subject"]
    time <- data[,"time"]
    if (is.unsorted(subject)){
        warning("Data are not ordered by subject ID and time - output will be ordered.")
        data <- data[order(subject, time),]
    }
    if (any(duplicated(data[,c("subject", "time")]))){
        warning("Data contain duplicated observation times for a subject - removing duplicates.")
        data <- data[!duplicated(data[,c("subject", "time")]), ]
    }
    subject <- data[,"subject"]; time <- data[,"time"];
    cens <- if (any(colnames(data)=="cens")) data[,"cens"] else rep(0, length(subject))
    msm.check.times(time, subject)
    times <- split(time, subject)
    cens <- split(cens, subject)
    n <- length(unique(subject))

### Covariates on intensities
    covnames <- names(covariates)
    ncovs <- length(covnames)
    misscovs <- setdiff(covnames, names(data))
    if (length(misscovs) > 0)
        stop("Covariates ", paste(misscovs, collapse=", "), " not found in data")
    covs <- if (ncovs > 0) lapply(split(data[,covnames], subject), as.matrix) else NULL
    allcovs <- covnames

### Covariates on misclassification
    misccovnames <- names(misccovariates)
    nmisccovs <- length(misccovnames)
    misscovs <- setdiff(misccovnames, names(data))
    if (length(misscovs) > 0)
        stop("Misclassification covariates ", paste(misscovs, collapse=", "), " not found in data")
    misccovs <- if (nmisccovs > 0) lapply(split(data[,setdiff(misccovnames,covnames)], subject), as.matrix) else NULL
    allmisccovs <- misccovnames

### Covariates on HMM
    if (!is.null(hcovariates)) {
        if (is.null(hmodel)) stop("hcovariates specified, but no hmodel")
        hcovnames <- unique(names(unlist(hcovariates)))
        if (length(hcovariates) != nstates)
            stop("hcovariates of length ", length(hcovariates), ", expected ", nstates)
        msm.check.hmodel(hmodel, nstates)
        misscovs <- setdiff(hcovnames, names(data))
        if (length(misscovs) > 0)
            stop("Covariates ", paste(misscovs, collapse=", "), " not found in data")
        hcovs <- lapply(split(data[,setdiff(hcovnames, setdiff(misccovnames,covnames))], subject), as.matrix)
    }
    else hcovs <- hcovnames <- NULL

### Extra variables to return
    extravars <- setdiff(names(data), unique(c("subject","time", "cens", covnames, misccovnames, hcovnames)))
    extradat <- if(length(extravars)>0) split(data[,extravars,drop=FALSE], subject) else NULL

### Starting states
    if (missing(start)) start <- rep(1, n)
    else if (length(start) == 1) start <- rep(start, n)
    else if (length(start) != n)
        stop("Supplied ", length(start), " starting states, expected 1 or ", n)

    nq <- length(qmatrix[qmatrix > 0])
    misspeccovs <- covnames[sapply(covariates, length) != nq]
    if (length(misspeccovs) > 0)
        stop("Initial values for covariates ", paste(misspeccovs, collapse=", "), " should be of length ", nq)
    beta <- do.call("rbind", as.list(covariates))
    ne <- length(ematrix[ematrix > 0])
    misspeccovs <- covnames[sapply(misccovariates, length) != ne]
    if (length(misspeccovs) > 0)
        stop("Initial values for misclassification covariates ", paste(misspeccovs, collapse=", "), " should be of length ", ne)
    beta.misc <- do.call("rbind", as.list(misccovariates))

### Check death argument. Logical values allowed for backwards compatibility
### (TRUE means final state is death, FALSE means no death state)
    statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
    if (is.logical(death) && death==TRUE)  {death <- nstates}
    else if (is.logical(death) && death==FALSE) {death <- 0}
    else if (length(setdiff(unique(death), 1:nstates)) > 0)
        stop(paste("Death states indicator contains states not in",statelist))

### Simulate a realisation for each person
    state <- numeric()
    keep.data <- numeric()
    subj <- split(subject, subject)
    subj.num <- match(subject, unique(subject))
    for (pt in 1:n)
    {
        sim.mod <- sim.msm(qmatrix, max(times[[pt]]), covs[[pt]], beta, times[[pt]], start[pt], min(times[[pt]]))
        obsd <- getobs.msm(sim.mod, times[[pt]], death, drop.absorb)
        pt.data <- cbind(subj[[pt]][obsd$keep], obsd$time, obsd$state, cens[[pt]][obsd$keep])
        if (!is.null(covnames)) pt.data <- cbind(pt.data, covs[[pt]][obsd$keep,,drop=FALSE])
        if (!is.null(misccovnames)) pt.data <- cbind(pt.data, misccovs[[pt]][obsd$keep,,drop=FALSE])
        if (!is.null(hcovnames)) pt.data <- cbind(pt.data, hcovs[[pt]][obsd$keep,,drop=FALSE])
        if (!is.null(extravars)) pt.data <- cbind(pt.data, extradat[[pt]][obsd$keep,,drop=FALSE])
        inds <- which(subj.num==pt)
        pt.data <- cbind(pt.data, keep=inds[obsd$keep])
        keep.data <- rbind(keep.data, pt.data)
      }
    colnames(keep.data) <- c("subject","time","state","cens",union(union(covnames,misccovnames),hcovnames),extravars,"keep")
    keep.data <- as.data.frame(keep.data)

### Simulate some misclassification or a HMM conditionally on the underlying state
    if (!is.null(ematrix)) {
        if (!all(dim(ematrix) == dim(qmatrix)))
            stop("Dimensions of qmatrix and ematrix should be equal")
        keep.data <- cbind(keep.data, obs=simmisc.msm(keep.data$state, ematrix, beta.misc, keep.data[,misccovnames,drop=FALSE]))
    }
    else if (!is.null(hmodel))
        keep.data <- cbind(keep.data, obs=simhidden.msm(keep.data$state, hmodel, nstates, hcovariates, keep.data[,hcovnames,drop=FALSE]))

### Replace state at censor times by censoring indicators
    censor <- unique(keep.data$cens[keep.data$cens != 0])
    if (is.null(censor.states)) censor.states <- 1:(nstates-1)
    if (!is.null(keep.data$obs))
      keep.data$obs <- ifelse(keep.data$cens > 0 & keep.data$state %in% censor.states, keep.data$cens, keep.data$obs)
    else
      keep.data$state <- ifelse(keep.data$cens > 0 & keep.data$state %in% censor.states, keep.data$cens, keep.data$state)
    keep.data$cens <- NULL

#    attr(keep.data, "keep") <- obsd$keep
    keep.data
}


## Simulate misclassification conditionally on an underlying state

simmisc.msm <- function(state, ematrix, beta, misccovs)
{
    ostate <- state
    if (is.null(ematrix))
        warning("No misclassification matrix given, assuming no misclassification")
    else {
        if (any(ematrix < 0)) stop("Not all elements of ematrix are > 0")
        if (any(ematrix > 1)) stop("Not all elements of ematrix are < 1")
        if (nrow(ematrix) != ncol(ematrix)) stop("Number of rows and columns of ematrix are not equal")
        nstates <- nrow(ematrix)
        ematrix <- msm.fixdiag.ematrix(ematrix)
        ostate <- state
        beta.states <- t(row(ematrix))[t(ematrix) > 0 & !(row(ematrix)==col(ematrix))]
        for (i in 1:nstates)
          if (any(state==i)) {
                n <- length(state[state==i])
                if (!is.null(beta)) { # covariates on misclassification probabilities
                    X <- as.matrix(misccovs[state==i,])
                    b <- beta[,beta.states == i,drop=FALSE]
                    p <- matrix(rep(ematrix[i,], n), nrow=n, byrow=TRUE)
                    mu <- log(p / p[,i])
                    emu <- array(0, dim=dim(p))
                    miscstates <- setdiff(which(ematrix[i,] > 0), i)
                    for (j in seq_along(miscstates))
                        emu[,miscstates[j]] <- exp(mu[,miscstates[j]] + X %*% b[,j])
                    emu[,i] <- 1
                    emu[,ematrix[i,]==0] <- 0
                    p <- emu / rowSums(emu)
                    for (j in 1:n)
                        ostate[state==i][j] <- resample(1:nstates, size=1, prob=p[j,], replace=TRUE)
                }
                else ostate[state==i] <- resample(1:nstates, size=n, prob=ematrix[i,], replace=TRUE)
            }
    }
    ostate
}

## Simulate HMM outcome conditionally on an underlying state

simhidden.msm <- function(state, hmodel, nstates, beta=NULL, x=NULL)
{
    y <- state
    msm.check.hmodel(hmodel, nstates)
    for (i in 1:nstates)
        if (any(state==i)) {
            if (inherits(hmodel[[i]], "hmmMVdist")) 
              stop("multivariate HMM outcomes not supported for simulation")
            ## don't change the underlying state if the HMM is the null (identity) model
            if (!(hmodel[[i]]$label=="identity" && (length(hmodel[[i]]$pars) == 0)))  {
                ## simulate from the sampling function "r" in the HMM object
                ## transform the location parameter by covariates if necessary
                rcall <- list(n=length(state[state==i]))
                if (!is.null(beta[[i]])) {
                    link <- get(hmodel[[i]]$link)
                    invlink <- get(.msm.INVLINK[hmodel[[i]]$link])
                    locpar <- .msm.LOCPARS[hmodel[[i]]$label]
                    loc <- hmodel[[i]]$pars[locpar]
                    loc <- invlink(link(loc) + as.matrix(x[state==i,names(beta[[i]])]) %*% beta[[i]])
                    rcall[[paste("r",locpar,sep="")]] <- loc
                }
                rfn <- hmodel[[i]]$r
                y[state==i] <- do.call("rfn", rcall)
            }
        }
    y
}


### Simulate data from fitted model with same observation scheme
### Used for parametric bootstrap in pearson.msm



#' Simulate from a Markov model fitted using msm
#' 
#' Simulate a dataset from a Markov model fitted using \code{\link{msm}}, using
#' the maximum likelihood estimates as parameters, and the same observation
#' times as in the original data.
#' 
#' This function is a wrapper around \code{\link{simmulti.msm}}, and only
#' simulates panel-observed data.  To generate datasets with the exact times of
#' transition, use the lower-level \code{\link{sim.msm}}.
#' 
#' Markov models with misclassified states fitted through the \code{ematrix}
#' option to \code{\link{msm}} are supported, but not general hidden Markov
#' models with \code{hmodel}.  For misclassification models, this function
#' includes misclassification in the simulated states.
#' 
#' This function is used for parametric bootstrapping to estimate the null
#' distribution of the test statistic in \code{\link{pearson.msm}}.
#' 
#' @param x A fitted multi-state model object as returned by \code{\link{msm}}.
#' @param drop.absorb Should repeated observations in an absorbing state be
#' omitted.  Use the default of \code{TRUE} to avoid warnings when using the
#' simulated dataset for further \code{\link{msm}} fits.  Or set to
#' \code{FALSE} if exactly the same number of observations as the original data
#' are needed.
#' @param drop.pci.imp In time-inhomogeneous models fitted using the \code{pci}
#' option to \code{\link{msm}}, censored observations are inserted into the
#' data by \code{\link{msm}} at the times where the intensity changes, but
#' dropped by default when simulating from the fitted model using this
#' function. Set this argument to \code{FALSE} to keep these observations and
#' the corresponding indicator variable.
#' @return A dataset with variables as described in \code{\link{simmulti.msm}}.
#' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
#' @seealso \code{\link{simmulti.msm}}, \code{\link{sim.msm}},
#' \code{\link{pearson.msm}}, \code{\link{msm}}.
#' @keywords models
#' @export simfitted.msm
simfitted.msm <- function(x, drop.absorb=TRUE, drop.pci.imp=TRUE){
    sim.df <- x$data$mf
    x$data <- expand.data(x)
    sim.df$"cens" <- ifelse(sim.df$"(state)" %in% 1:x$qmodel$nstates, 0, sim.df$"(state)") # 0 if not censored, cens indicator if censored, so that censoring is retained in simulated data.  TODO used in pearson?
    if (x$qcmodel$ncovs > 0) {
        sim.df <- cbind(sim.df, x$data$mm.cov)
        cov.effs <- lapply(x$Qmatrices, function(y)t(y)[t(x$qmodel$imatrix)==1])[x$qcmodel$covlabels]
    } else cov.effs <- NULL
    if (x$ecmodel$ncovs > 0) {
        sim.df <- cbind(sim.df, x$data$mm.mcov)
        misccov.effs <- lapply(x$Ematrices, function(y)t(y)[t(x$emodel$imatrix)==1])[x$ecmodel$covlabels]
    } else misccov.effs <- NULL
    names(sim.df) <- replace(names(sim.df), match(c("(time)","(subject)"), names(sim.df)),
                             c("time","subject"))
    if (!is.null(sim.df$"(subject.weights)")) names(sim.df)[names(sim.df)=="(subject.weights)"] = "subject.weights"
    sim.df$state <- NULL # replace observed with simulated state
    if (any(union(names(cov.effs), names(misccov.effs)) %in% c("state","time","subject","cens")))
        stop("Not supported with covariates named \"state\", \"time\", \"subject\" or \"cens\"") # TODO?
    boot.df <- simmulti.msm(data=sim.df,
                            qmatrix=qmatrix.msm(x, covariates=0, ci="none"),
                            covariates=cov.effs,
                            death=FALSE,
                            ematrix=ematrix.msm(x, covariates=0, ci="none"),
                            misccovariates=misccov.effs,
                            drop.absorb=drop.absorb
                            )
    if (drop.pci.imp & !is.null(boot.df$"(pci.imp)")) {
        boot.df <- boot.df[!boot.df$"(pci.imp)",]
        boot.df$"(pci.imp)" <- NULL
    }
    boot.df
}
