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
    diag(qmatrix) <- rep(0, nstates); diag(qmatrix) <- - apply(qmatrix, 1, sum)

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
            cur.qmatrix <- t(cur.qmatrix)
            diag(cur.qmatrix) <- rep(0, nstates); diag(cur.qmatrix) <- - apply(cur.qmatrix, 1, sum)
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
      absorbing <- which(diag(sim$qmatrix == 0))
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


### Simulate a multi-state model dataset using fixed observation times

simmulti.msm <- function(data,           # data frame with subject, times, covariates... 
                         qmatrix,        # intensity matrix
                         beta = NULL,    # list of covariate effects on log intensities
                         death = FALSE,  # vector of indicators for "death" states, ie absorbing states whose entry time is known exactly,
                                        # but with unknown transient state at previous instant
                         tunit = 1.0 # no longer used
                         )
  {
      if (!("subject" %in% names(data)))
        stop("\"subject\" column missing from data")
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
      covnames <- setdiff(names(data), c("subject","time"))
      ncovs <- length(covnames)
      times <- split(time, subject)
      covs <- if (ncovs > 0) lapply(split(data[,covnames], subject), as.matrix) else NULL
      n <- length(unique(subject))
            
### Check consistency of qmatrix
      nstates <- dim(qmatrix)[1]
      if (nstates != dim(qmatrix)[2])
        stop("Number of rows and columns of qmatrix should be equal")
      diag(qmatrix) <- rep(0, nstates)
      if (any ( qmatrix < 0 ) )
        stop("qmatrix should not have negative off-diagonal elements")

### Check covariate effects
      if (ncovs > 0) {
          if (is.null(beta)) stop("Covariate effects \"beta\" not provided")
          if (!is.matrix(beta))
            beta <- matrix(beta, ncol=length(beta))
          if (nrow(beta) != ncovs)
            stop(paste("Expected",ncovs,"rows in beta corresponding to different covariates, found", nrow(beta)))
          if (!is.null(ncol(beta)) & ncol(beta) != length(qmatrix[qmatrix > 0]))
            stop(paste("Expected",length(qmatrix[qmatrix > 0]),"columns in covariate matrix, found", ncol(beta)))
      }

### Check death argument. Logical values allowed for backwards compatibility (TRUE means final state is death, FALSE means no death state)       
      statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
      if (is.logical(death) && death==TRUE)  {death <- nstates}
      else if (is.logical(death) && death==FALSE) {death <- 0}
      else if (length(setdiff(unique(death), 1:nstates)) > 0)
        stop(paste("Death states indicator contains states not in",statelist))
      if (!missing(tunit)) warning("tunit argument is no longer used. death times are now assumed exact")
      
### Simulate a realisation for each person
      state <- numeric()
      keep.data <- numeric()
      subj <- split(subject, subject)
      for (pt in 1:n)
        {
            sim.mod <- sim.msm(qmatrix, max(times[[pt]]), covs[[pt]], beta, times[[pt]], 1, min(times[[pt]]))
            obsd <- getobs.msm(sim.mod, times[[pt]], death)
            keep.data <- rbind(keep.data,
                               cbind(subj[[pt]][obsd$keep], obsd$time, covs[[pt]][obsd$keep,], obsd$state))
        }
      colnames(keep.data) <- c("subject","time",covnames,"state")
      keep.data
  }
