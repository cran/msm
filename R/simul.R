### FUNCTIONS FOR SIMULATING FROM MULTI-STATE MODELS

### General function to simulate one individual's realisation from a continuous-time Markov model
### Produces the exact times of transition

sim.msm <- function(qmatrix,   # intensity matrix
                    maxtime,   # maximum time for realisations
                    covs=NULL,     # covariate matrix
                    beta=NULL,     # matrix of cov effects on qmatrix
                    obstimes=NULL, # times at which time-dependent covariates change
                    start = 1,     # starting state
                    mintime = 0    # time to start from 
                    )
{
    nstates <- dim(qmatrix)[1]
    simstates <- start
    simtimes <- mintime
    if (!is.null(covs)) covs <- as.matrix(covs)
    i <- 1
    absorb <- (sum ( qmatrix[1, -1] ) == 0)
    qmatrix <- msm.fixdiag.qmatrix(qmatrix)

### Assume that time-dependent covariates are constant in between observation times.
### Gets the nearest covariate value prior to t to use to compute the intensity matrix at t 
    getQcov <- function(cur.time, # Current time of Markov process
                        obstimes, # Observation times of Markov process
                        covs      # Covariate matrix at the observation times
                        )
      {
          ret <- NULL
          covs <- as.matrix(covs)
          for (i in 1 : (length(obstimes)-1) )
            if ( (cur.time >= obstimes[i]) && (cur.time < obstimes[i+1]) )
              ret <- covs[i,]
          if (is.null(ret)) stop ("Observation times inconsistent with mintime")
          ret                     # Returns the covariate value at the previous observation time
      }
    
    while ( (simtimes[i] < maxtime) & (!absorb) ){
        cur.st <- simstates[i]
        cur.t <- simtimes[i]        
        if (!is.null(covs)){
            cur.qmatrix <- t(qmatrix)
            cur.cov <- as.matrix(getQcov(cur.t, obstimes, covs))
            intens <- cur.qmatrix[cur.qmatrix > 0]
            intens <- intens * exp( t(beta) %*% cur.cov )
            cur.qmatrix[cur.qmatrix > 0] <- intens
            cur.qmatrix <- msm.fixdiag.qmatrix(t(cur.qmatrix))
        }
        else cur.qmatrix <- qmatrix
        absorb <- (sum ( cur.qmatrix[cur.st, -cur.st] ) == 0)
        if (!absorb) {
            nextprobs <- cur.qmatrix[cur.st, ] / sum ( cur.qmatrix[cur.st, -cur.st] )
            nextprobs[cur.st] <- 0
            cumprobs <- cumsum(nextprobs)
            nextstate <- min( (1:nstates) [runif(1) < cumprobs])
            nextlag <- rexp (1,  - cur.qmatrix[cur.st, cur.st] )
            simstates <- c(simstates, nextstate)
            simtimes <- c(simtimes, cur.t + nextlag)
            i <- i+1
        }
    }
    list(states = simstates, times = simtimes, qmatrix = qmatrix)
}


### Given a simulated Markov model, get the current state at various observation times
### Only keep one observation in the absorbing state

getobs.msm <- function(sim,          # output from simMSM
                       obstimes,     # fixed observation times
                       death = FALSE, # indicators for death states
                       tunit = 1.0   # observation time unit in days (e.g. time in months, tunit = 30)
                       ) 
  {
      censtime <- max(obstimes)
      nsim <- length(sim$states)
      nobs <- length(obstimes)
      obsstate <- numeric()
      keep.time <- numeric()
      nstates <- dim(sim$qmatrix)[1]
      absorbing <- absorbing.msm(qmatrix=sim$qmatrix)
      absorbed <- FALSE
      cur.j <- 1
      for (i in 1:nobs) {
          if (absorbed) break
          j <- cur.j
          found <- FALSE
          while (!found & (j < nsim)){
              if ((obstimes[i] >= sim$times[j]) & (obstimes[i] < sim$times[j+1]) ){
                  keep.time[i] <- i
                  obsstate[i] <- sim$states[j]
                  cur.j <- j
                  found <- TRUE
              }
              else if ((obstimes[i] >= max(sim$times)) && (sim$states[j+1] %in% absorbing)){
                  obsstate[i] <- nstates
                  if (sim$states[j+1] %in% death)
                    obstimes[i] <- sim$times[j+1] # death times are observed exactly. 
                  keep.time[i] <- i
                  absorbed <- TRUE
                  found <- TRUE
              }
              j <- j+1
          }
      }
      list(state = obsstate, time = obstimes[keep.time], keep = keep.time)
  }


### Simulate a multi-state Markov or hidden Markov model dataset using fixed observation times

### Would it be better to make specification of covariate model consistent with model fitting function? 
### e.g. separate hcovariates and  covariates formulae,
### plus covinits and hcovinits? 

simmulti.msm <- function(data,           # data frame with subject, times, covariates... 
                         qmatrix,        # intensity matrix
                         covariates=NULL,  # initial values
                         death = FALSE,  # vector of indicators for "death" states, ie absorbing states whose entry time is known exactly,
                                        # but with unknown transient state at previous instant
                         start,         # starting states of the process, defaults to all 1.
                         ematrix = NULL,# misclassification matrix
                         hmodel = NULL,  # hidden Markov model formula
                         hcovariates = NULL   # covariate effects on hidden Markov model response distribution
                         )
  {
            
### Check consistency of qmatrix and covariate inits
      nstates <- dim(qmatrix)[1]
      msm.check.qmatrix(qmatrix)
      qmatrix <- msm.fixdiag.qmatrix(qmatrix)

### Subject, time and state
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
      subject <- data[,"subject"]; time <- data[,"time"]
      msm.check.times(time, subject)
      times <- split(time, subject)
      n <- length(unique(subject))

### Covariates on intensities 
      covnames <- names(covariates)
      ncovs <- length(covnames)
      misscovs <- setdiff(covnames, names(data))
      if (length(misscovs) > 0)
        stop("Covariates ", paste(misscovs, collapse=", "), " not found in data")
      covs <- if (ncovs > 0) lapply(split(data[,covnames], subject), as.matrix) else NULL
      allcovs <- covnames

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
          hcovs <- lapply(split(data[,setdiff(hcovnames, covnames)], subject), as.matrix)
      }
      else hcovs <- hcovnames <- NULL
      
### Starting states
      if (missing(start)) start <- rep(1, n)
      else if (length(start) != n)
        stop("Supplied ", length(start), " starting states, expected ", n)

      nq <- length(qmatrix[qmatrix > 0])
      misspeccovs <- covnames[sapply(covariates, length) != nq]
      if (length(misspeccovs) > 0)
        stop("Initial values for covariates ", paste(misspeccovs, collapse=", "), " should be of length ", nq)
      beta <- do.call("rbind", as.list(covariates))

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
      for (pt in 1:n)
        {
            sim.mod <- sim.msm(qmatrix, max(times[[pt]]), covs[[pt]], beta, times[[pt]], start[pt], min(times[[pt]]))
            obsd <- getobs.msm(sim.mod, times[[pt]], death)
            keep.data <- rbind(keep.data,
                               cbind(subj[[pt]][obsd$keep], obsd$time,
                                     covs[[pt]][obsd$keep,], hcovs[[pt]][obsd$keep,], obsd$state))
        }
      colnames(keep.data) <- c("subject","time",covnames,setdiff(hcovnames, covnames),"state")
      keep.data <- as.data.frame(keep.data)

### Simulate some misclassification or a HMM conditionally on the underlying state
      if (!missing(ematrix))
        keep.data <- cbind(keep.data, obs=simmisc.msm(keep.data$state, ematrix))
      else if (!missing(hmodel))
        keep.data <- cbind(keep.data, obs=simhidden.msm(keep.data$state, hmodel, nstates, hcovariates, keep.data[,hcovnames,drop=FALSE]))
      keep.data 
  }

## Simulate misclassification conditionally on an underlying state

simmisc.msm <- function(state, ematrix)
  {
      ostate <- state
      if (is.null(ematrix))
        warning("No misclassification matrix given, assuming no misclassification")
      else {
          if (any(ematrix) < 0) stop("Not all elements of ematrix are > 0")
          if (any(ematrix) > 1) stop("Not all elements of ematrix are < 1")
          if (nrow(ematrix) != ncol(ematrix)) stop("Number of rows and columns of ematrix are not equal")
          nstates <- nrow(ematrix)
          ematrix <- msm.fixdiag.ematrix(ematrix)
          ostate <- state
          for (i in seq(nstates))
            if (any(state[state==i]))
              ostate[state==i] <- sample(seq(nstates), size=length(state[state==i]), prob=ematrix[i,], replace=TRUE)
      }
      ostate
  }

## Simulate HMM outcome conditionally on an underlying state

simhidden.msm <- function(state, hmodel, nstates, beta=NULL, x=NULL)
  {
      y <- state
      msm.check.hmodel(hmodel, nstates)
      for (i in seq(nstates))
        if (any(state==i)) {
            ## don't change the underlying state if the HMM is the null (identity) model
            if (!(hmodel[[i]]$label=="identity" && (length(hmodel[[i]]$pars) == 0)))  {
                ## simulate from the sampling function "r" in the HMM object
                ## transform the location parameter by covariates if necessary
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
