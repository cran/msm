### Function to fit Markov multi-state models in continuous time
### with either arbitrary observation times or observed exact transition times
### with or without misclassification between true and underlying states

msm <- function(formula,   # formula with  observed Markov states   ~  observation times (required)
                qmatrix,    # matrix of 1s and 0s with indices of allowed transitions (diagonal is ignored) (required)
                misc = FALSE,
                ematrix = NULL,    # matrix of 1s and 0s with indices of allowed misclassfications (diagonal is ignored) (required)
                inits,      # initial values of optimisation or fixed values (required)
                subject = NULL, # optional, defaults to all the same if not given
                covariates = NULL, # formula specifying covariates on transition rates.
                constraint = NULL, # which intensities have covariates on them (as in Marshall et al.)
                misccovariates = NULL, # formula specifying covariates on misclassification probs
                miscconstraint = NULL, # which misc probs have covariates on them (as in Marshall et al.)
                covmatch = "previous",   # take the covariate value from the previous or next observation
                initprobs = NULL,  # initial state occupancy probabilities
                data=list(),       # optional data frame in which to interpret variable names
                fromto = FALSE,
                fromstate, #
                tostate,   #  data required if fromto is TRUE
                timelag,   #
                death = FALSE,  # if final state is death and death date known within a day
                tunit = 1.0, # unit in days of the given time variable (death date assumed known within a day)
                exacttimes = FALSE,
                fixedpars = NULL, # specify which parameters to fix
                ... # options to optim
                )
  {

### INTENSITY MATRIX
      nstates <- dim(qmatrix)[1]
      qmatrix <- msm.check.qmatrix(qmatrix)
      qvector <- as.numeric(t(qmatrix)) # matrix to vector by filling rows first
      nintens <- sum(qmatrix)

### MISCLASSIFICATION MATRIX
      if (misc) {
          if (missing(ematrix)) stop("Misclassification matrix not given")
          ematrix <- msm.check.ematrix(ematrix, qmatrix)
          evector <- as.numeric(t(ematrix)) # matrix to vector by filling rows first
          nmiscs <- sum(ematrix)
          if (is.null(initprobs))
            initprobs <- c(1, rep(0, nstates-1))
      }
      else {
          nmiscs <- 0; ematrix <- evector <- NULL
      }

### BASIC DATA: BY TRANSITION PAIRS
      if (fromto) {
          if (misc)
            stop("\"fromto\" style of data not allowed for misclassification model")
          state <- if (!missing(fromstate)) eval(substitute(fromstate), data, parent.frame()) else stop("From-state not given")
          tostate <- if (!missing(tostate)) eval(substitute(tostate), data, parent.frame()) else stop("To-state not given")
          time <- if (!missing(timelag)) eval(substitute(timelag), data, parent.frame()) else stop("Time lag not given")
          subject <- NULL
          nobs <- length(state)
          mf <- na.omit(data.frame(state=state, tostate=tostate, time=time, ind=1:nobs))
          state <- mf$state; tostate <- mf$tostate; time <- mf$time
          statetimerows.kept <- mf$ind
          subjrows.kept <- 1:nobs
      }
      
### BASIC DATA: BY OBSERVATION TIME
      else { 
          ## form vectors with subject ID, state and corresponding observation time
          ## Store subject internally as a factor. 
          mf <- model.frame(formula, data=data)
          state <- mf[,1]
          time <- mf[,2]
          droprows <- as.numeric(attr(mf, "na.action"))
          nobs <- length(c(state, droprows))
          statetimerows.kept <- (1:nobs)[! ((1:nobs) %in% droprows)]
          subject <-
            if (missing(subject)) rep(1, dim(mf)[1])
            else eval(substitute(subject), data, parent.frame())
          subjrows.kept <- (1:nobs) [!is.na(subject)]
          subject <- factor(na.omit(subject))
          fromstate <- tostate <- timelag <- NULL
      }
      npts <- length(unique(subject))      
      covrows.kept <- misccovrows.kept <- 1:nobs

### COVARIATES
      if (!is.null(covariates)) {
          pc <- msm.process.covs(covariates, data, constraint, nobs, nintens)
          covvec <- pc$covvec;  constrvec <- pc$constrvec;  covlabels <- pc$covlabels
          ncovs <- pc$ncovs;  ncoveffs <- pc$ncoveffs; covrows.kept <- pc$kept.rows
          covmeans <- pc$covmeans; covfactor <- pc$covfactor
      }
      else {
          ncovs <- ncoveffs <- 0
          covvec <- constrvec <- covlabels <- covmeans <- covfactor <- NULL
      }      

### MISCLASSIFICATION COVARIATES
      if (!is.null(misccovariates) & (misc)) {
          pc <- msm.process.covs(misccovariates, data, miscconstraint, nobs, nmiscs)
          misccovvec <- pc$covvec;  miscconstrvec <- pc$constrvec;  misccovlabels <- pc$covlabels
          nmisccovs <- pc$ncovs;  nmisccoveffs <- pc$ncoveffs; misccovrows.kept <- pc$kept.rows
          misccovmeans <- pc$covmeans; misccovfactor <- pc$covfactor
      }
      else if (is.null (misccovariates) | (!misc)) {
          if (!is.null (misccovariates))
            warning("misccovariates have been specified, but misc is FALSE. Ignoring misccovariates.")
          nmisccovs <- nmisccoveffs <- 0
          misccovvec <- miscconstrvec <- misccovlabels <- misccovmeans <- misccovfactor <- NULL
      }

### DROP MISSING DATA      
      final.rows <- intersect(intersect(statetimerows.kept, subjrows.kept),
                              intersect(covrows.kept, misccovrows.kept))
      if (!fromto) subject <- factor(subject[ subjrows.kept %in% final.rows ])
      time <- time[ statetimerows.kept %in% final.rows ]
      state <- state[ statetimerows.kept %in% final.rows ]
      if (fromto) tostate <- tostate[ statetimerows.kept %in% final.rows ]
      if (ncovs>0) covvec <- as.numeric(t(matrix(covvec, ncol=ncovs)[covrows.kept %in% final.rows]))
      if (nmisccovs>0) misccovvec <- as.numeric(t(matrix(misccovvec, ncol=nmisccovs)[misccovrows.kept %in% final.rows]))
      nmiss <- length(setdiff(1:nobs, final.rows))
      plural <- if (nmiss==1) "" else "s"
      if (nmiss > 0) warning(paste ( nmiss, " record", plural, " dropped due to missing values",sep=""))
      msm.check.consistency(qmatrix, misc, fromto=fromto, state=state, subject=subject, time=time,
                            tostate=tostate, death=death, exacttimes=exacttimes, tunit=tunit)
      nobs <- length(state)

### FORM LIST OF INITIAL PARAMETERS
      npars <- nintens + ncoveffs + nmiscs + nmisccoveffs
      plabs <- c(rep("qbase",nintens), rep("qcov", ncoveffs), rep("ebase",nmiscs), rep("ecov",nmisccoveffs))
      if (npars != length(inits)) {
          err.msg <- paste(length(inits),"initial values supplied, expected",
                           nintens,"intensities,",
                           ncoveffs,"covariate effects",
                           if (misc) paste(", ",nmiscs,"misclassification probabilities and",
                                           nmisccoveffs,"misclassification covariate effects"))
          stop(err.msg)
      }
      inits[plabs=="qbase"] <- log(inits[plabs=="qbase"]) ## optimise transition rates on the log scale
      inits[plabs=="ebase"] <- logit(inits[plabs=="ebase"]) ## optimise misclassification probs on the logit scale
      covmatch <- if (covmatch=="previous") 0 else if (covmatch=="next") 1 else stop ("covmatch should be \"previous\" or \"next\"")

### FIX SOME PARAMETERS IF REQUIRED
      notfixed <- setdiff(1:npars, fixedpars)
      fixdiff <- setdiff(fixedpars, 1:npars)
      if (length(fixdiff) > 0)
        stop ( paste ("Elements of fixedpars should be in 1, ..., ",npars, ", fixedpars contains",paste(fixdiff,collapse=" ") ) )
      allinits <- inits
      if (length(fixedpars)==length(inits)) ## all parameters are fixed
        fixed <- TRUE
      else {
          inits <- inits[notfixed]
          fixed <- FALSE
      }

### CALCULATE LIKELIHOOD AT INITIAL VALUES...
      if (fixed) {
          likval <- lik.msm(inits, allinits, misc, subject, time, state, tostate, fromto,
                            qvector, evector, covvec, constrvec, misccovvec, miscconstrvec, 
                            initprobs, nstates, nintens, nmiscs, nobs, npts,
                            ncovs, ncoveffs, nmisccovs, nmisccoveffs, covmatch,
                            death, tunit, exacttimes, fixedpars, plabs)
          likval <- likval$endlik
          params <- inits
          covmat <- NULL
          foundse <- FALSE
      }

### ... OR DO MAXIMUM LIKELIHOOD ESTIMATION
      else {
          opt <- optim(inits, lik.msm, hessian=TRUE,  ...,# arguments to optim
                       allinits=allinits, misc=misc, subject=subject, time=time,  # arguments to lik
                       state=state, tostate=tostate, fromto=fromto,
                       qvector=qvector, evector=evector, covvec=covvec, constrvec=constrvec,
                       misccovvec=misccovvec, miscconstrvec=miscconstrvec, initprobs=initprobs,
                       nstates=nstates, nintens=nintens, nmiscs=nmiscs, nobs=nobs, npts=npts,
                       ncovs=ncovs, ncoveffs=ncoveffs, nmisccovs=nmisccovs, nmisccoveffs=nmisccoveffs,
                       covmatch=covmatch, death=death, tunit=tunit, exacttimes=exacttimes,
                       fixedpars=fixedpars, plabs=plabs)
          params <- allinits
          params[notfixed] <- opt$par
          if (all(eigen(opt$hessian)$values > 0)) {
              covmat <- matrix(0, nrow=npars, ncol=npars)
              covmat[notfixed,notfixed] <- solve(0.5 * opt$hessian) 
              foundse <- TRUE
          }
          else {
              covmat <- "Non-positive definite approximate variance-covariance"
              foundse <- FALSE
          }
      }
      estimates.t <- params  # transformed estimates
      estimates.t[plabs=="qbase"] <- exp(params[plabs=="qbase"])  
      estimates.t[plabs=="ebase"] <- expit(params[plabs=="ebase"])

### REARRANGE THE VECTOR OF PARAMETER ESTIMATES (LOG-INTENSITIES, MISC PROBS AND
### COVARIATE EFFECTS) INTO LISTS OF MATRICES
      output <- msm.form.output(qmatrix, nstates, nintens,
                                npars, ncovs, constrvec, fixedpars, fixed, covlabels,
                                params, covmat, foundse, 0)
      Qmatrices <- output$Matrices
      QmatricesSE <- if (fixed) NULL else output$MatricesSE 

      if (misc) {
          output <- msm.form.output(ematrix, nstates, nmiscs,
                                    npars, nmisccovs, miscconstrvec, fixedpars, fixed, misccovlabels,
                                    params, covmat, foundse,
                                    nintens+ncoveffs)
          Ematrices <- output$Matrices
          EmatricesSE <- if (fixed) NULL else output$MatricesSE
          names(Ematrices)[1] <- "logitbaseline"
          if (foundse & !fixed) names(EmatricesSE)[1] <- "logitbaseline"
      }
      else {
          Ematrices <- EmatricesSE <- NULL
      }

### FORM A MSM OBJECT FROM THE RESULTS
      
      msmobject <- list (
                         misc = misc,
                         Qmatrices = Qmatrices, 
                         QmatricesSE = QmatricesSE, 
                         minus2loglik = if (fixed) likval else opt$value,
                         estimates = params,
                         estimates.t = estimates.t,
                         covmat = covmat,
                         foundse = foundse,
                         data = list(nobs=nobs, npts=npts, state=state, time=time, subject=subject,
                           fromto=fromto, tostate=tostate,
                           covvec=covvec, covmeans=covmeans, 
                           covlabels=covlabels, covfactor=covfactor, ncovs=ncovs,
                           nmisccovs=nmisccovs, tunit=tunit),
                         model = list(qvector=qvector, evector=evector,
                           nstates=nstates, nintens=nintens, nmiscs=nmiscs, 
                           constrvec=constrvec,  ncoveffs=ncoveffs,
                           covmatch=covmatch, initprobs=initprobs, death=death, 
                           exacttimes=exacttimes)
                         )
      
      attr(msmobject, "fixed") <- fixed
      class(msmobject) <- "msm"

      q <- qmatrix.msm(msmobject) # intensity matrix with centered covariates
      msmobject$Qmatrices$baseline <- q$estimates
      msmobject$QmatricesSE$baseline <- q$SE
      if (misc) {
          msmobject$Ematrices <- Ematrices
          msmobject$EmatricesSE <- EmatricesSE
          e <- ematrix.msm(msmobject) # misc matrix with centered covariates
          msmobject$Ematrices$baseline <- e$estimates
          msmobject$EmatricesSE$baseline <- e$SE
          msmobject$model$nmisccoveffs <- nmisccoveffs
          msmobject$model$initprobs <- initprobs
      }
      if (nmisccovs > 0){
          msmobject$data$misccovvec <- misccovvec;  msmobject$data$misccovmeans <- misccovmeans
          msmobject$data$misccovlabels <- misccovlabels;  msmobject$model$miscconstrvec <- miscconstrvec
          msmobject$data$misccovfactor <- misccovfactor
      }
      
      ## Calculate mean sojourn times with centered covariates
      msmobject$sojourn <- sojourn.msm(msmobject)
      
      msmobject
  }

expit <- function(x){exp(x) / (1 + exp(x))}
logit <- function(x){log (x / (1 - x)) }

### Wrapper for the C code which evaluates the -2*log-likelihood for a Markov multi-state model with misclassification
### This is optimised by nlm

lik.msm <- function(params, allinits, misc, subject, time, state, tostate, fromto, qvector, evector, covvec, constrvec, misccovvec, miscconstrvec, 
                    initprobs, nstates, nintens, nmiscs, nobs, npts, ncovs, ncoveffs, nmisccovs, nmisccoveffs, covmatch,
                    death, tunit, exacttimes, fixedpars, plabs)
  {      
      p <- length(params)
      state <- state - 1  # In R, work with states 1, ... n. In C, work with states 0, ... n-1
      npars <- nintens + ncoveffs + nmiscs + nmisccoveffs
      notfixed <- setdiff(1:npars, fixedpars)
      optplabs <- plabs[notfixed]
      if (!is.null(fixedpars))  fixedpars <- fixedpars - 1
      else fixedpars <- -1
      params[optplabs=="qbase"] <- exp(params[optplabs=="qbase"])   ## optimise transition intensities on log scale
      params[optplabs=="ebase"] <- expit(params[optplabs=="ebase"]) ## optimise misclassification probs on logit scale
      allinits[plabs=="qbase"] <- exp(allinits[plabs=="qbase"])   
      allinits[plabs=="ebase"] <- expit(allinits[plabs=="ebase"])
      if (!misc) {
          evector <- misccovvec <- miscconstrvec <- initprobs <- nms <- NULL
          if (fromto)
            tostate <- tostate - 1
          nms <- 0
      }
      else {
          ematrix <- t(matrix(evector, nrow=nstates))
          esum <- apply(ematrix, 1, sum)
          nms <- max ( seq(along=esum)[esum > 0])
      }
      do.what <- 1
      lik <- .C("msmCEntry",
                as.integer(do.what),
                as.double(params),
                as.double(allinits),
                as.integer(misc),
                as.integer(p),
                as.integer(subject),
                as.double(time),
                as.integer(state),
                as.integer(tostate),
                as.integer(fromto),
                as.integer(qvector),
                as.integer(evector),
                as.double(covvec),
                as.integer(constrvec),
                as.double(misccovvec),
                as.integer(miscconstrvec),
                as.double(initprobs),
                as.integer(nstates),
                as.integer(nms),
                as.integer(nintens),
                as.integer(nmiscs),
                as.integer(nobs),
                as.integer(npts),
                as.integer(ncovs),
                as.integer(ncoveffs),
                as.integer(nmisccovs),
                as.integer(nmisccoveffs),
                as.integer(covmatch),
                as.integer(death),
                as.double(tunit),
                as.integer(exacttimes),
                as.integer(fixedpars),
                as.double(NULL),
                as.integer(NULL),
                endlik = double(1)
                )      
      lik$endlik
  }


### Process covariates and covariates constraints, in preparation for being passed to the likelilhood optimiser
### This function is called for both sets of covariates (transition rates and the misclassification probs)

msm.process.covs <- function(covariates, # formula:  ~ cov1 + cov2 + ...
                             data,
                             constraint,
                             nobs,       # number of observations including missing values
                             nmatrix     # number of transition intensities / misclassification probs
                             )
  {
      ## form covariate matrix to vectorise and pass to optimisation
      mm <- as.data.frame(model.matrix(covariates, data=data))
      covlabels <- names(mm)[-1]
      ncovs <- length(covlabels)
      mm <- mm[-1]
      mf <- model.frame(covariates, data=data)
      droprows <- as.numeric(attr(mf, "na.action"))
      kept.rows <- (1:nobs)[! ((1:nobs) %in% droprows)]
      ## centre the covariates about their means
      covmeans <- apply(mm, 2, mean)
      covstds <- apply(mm, 2, sd)
      covfactor <- sapply(mf, is.factor)
      mm <- sweep(mm, 2, covmeans)
      ##       mm <- sweep(mm, 2, covstds, "/")
      covvec <- unlist(mm)
      if (is.null(constraint))
        constrvec <- 1:(nmatrix*ncovs)
      else {
          ## check and parse the list of constraints on covariates
          for (i in names(constraint))
            if (!(is.element(i, covlabels))){
                factor.warn <- if (is.factor(data[, i]))
                  "\n\tFor factor covariates, specify constraints using covnameCOVVALUE = c(...)"
                else ""                  
                stop(paste("Covariate \"", i, "\" in constraint statement not in model.", factor.warn, sep=""))
            }
          constrvec <- numeric()
          maxc <- 0
          for (i in seq(along=covlabels)){
              ## build complete vectorised list of constraints for covariates in covariates statement
              ## so. e.g. constraints = (x1=c(3,3,4,4,5), x2 = (0.1,0.2,0.3,0.4,0.4))
              ##     turns into constrvec = c(1,1,2,2,3,4,5,6,7,7) with seven distinct covariate effects
              if (is.element(covlabels[i], names(constraint))) {
                  thiscon <- constraint[[match(covlabels[i], names(constraint))]]
                  if (length(thiscon) != nmatrix)
                    stop("\"",names(constraint)[i],"\"","constraint of length",length(constraint[i]),"should be",nmatrix)
                  constrvec <- c(constrvec, maxc + codes(factor(thiscon)))
                  maxc <- max(constrvec)
              }
              else constrvec <- c(constrvec, (i-1)*nmatrix + 1:nmatrix)
          }
      }
      ncoveffs <- max(unique(constrvec))
      list(covvec=covvec, # vector of concatenated covariate data
           constrvec=constrvec, 
           covlabels=covlabels, # covariate names
           covfactor=covfactor, # indicators for whether each covariate is a factor
           ncovs=ncovs,         # number of covariates
           ncoveffs=ncoveffs,   # number of distinct covariate effect parameters
           covmeans=covmeans,
           covstds=covstds,
           kept.rows = kept.rows # rows of original data file not containing missing covariate values
           )
  }


### Build a list of intensity / misclassification matrices from vector of parameter estimates
### One baseline, one corresponding to linear effect of each covariate

msm.form.output <- function(matrix, # matrix of 0/1 indicators for allowed intensities/misclassifications
                            nstates, # number of states
                            nmatrix, # number of allowed intensities/misclassifications
                            npars, ncovs, constrvec, fixedpars, fixed, covlabels, # ... inherited from msm
                            params, covmat, foundse, # ... inherited from msm
                            params.offset # starting index in full vector of parameters
                            )
  {
### REARRANGE THE OPTIMISED INTENSITIES AND COVARIATE EFFECTS
### 'Baseline' log intensities are with covariates equal to their means.
      Matrices <- list()
      MatricesSE <- list()
      fixtf <- rep(FALSE, npars)
      fixtf[fixedpars] <- TRUE
      ## build intensity matrix and matrix of SEs from vector of estimates
      for (i in 0:ncovs) {
          matrixname <- if (i==0) "logbaseline" else covlabels[i] # name of the current output matrix.
          mat <- t(matrix) # I fill matrices by row, while R fills them by column. Is this sensible...?
          ## the relevant elements of the parameter vector
          parinds <-
            if (i==0)
              params.offset + 1:nmatrix
            else
              (params.offset + nmatrix + constrvec[((i-1)*nmatrix + 1)  :  (i*nmatrix)])
          mat[t(matrix)==1] <- params[parinds]
          mat <- t(mat)
          dimnames(mat) <- list(paste("Stage",1:nstates), paste("Stage",1:nstates))
          if (foundse && !fixed){
              intenscov <- covmat[parinds, parinds]
              intensse <- sqrt(diag(intenscov))
              semat <- t(matrix)
              semat[t(matrix)==1] <- intensse 
              semat <- t(semat)
              diag(semat) <- rep(0, nstates)
              dimnames(semat)  <- list(paste("Stage",1:nstates), paste("Stage",1:nstates))
          }
          else if (!fixed){
              semat <- "Unstable Hessian at the estimates"
          }
          Matrices[[matrixname]] <- mat
          if (!fixed) MatricesSE[[matrixname]] <- semat
      }
      list(Matrices=Matrices,     # list of baseline log intensities/logit misc probability matrix
                                        # and linear effects of covariates
           MatricesSE=MatricesSE  # corresponding matrices of standard errors
           )
  }

### Check the consistency of the supplied data with the specified model 

msm.check.consistency <- function(qmatrix, misc, fromto=FALSE, subject=NULL,
                                  state, tostate=NULL, time,
                                  death=FALSE, exacttimes=FALSE, tunit=NULL)
  {
      msm.check.state(nrow(qmatrix), state, fromto, tostate)
      if (!misc)
        msm.check.model(state, subject, qmatrix, fromto, tostate)      
      if (!fromto)
        msm.check.times(time, subject)
      if (death & !exacttimes)
        msm.check.tunit(fromto, time, subject, timelag, tunit)
      invisible()
  }

### Check transition matrix indicators

msm.check.qmatrix <- function(qmatrix)
{
    qmatrix <- as.matrix(qmatrix)
    if (nrow(qmatrix) != ncol(qmatrix))
      stop("Number of rows and columns of qmatrix should be equal")
    diag(qmatrix) <- 0
    if (!all ( is.element(qmatrix, c(0,1)) ) )
      stop("Not all off-diagonal elements of qmatrix are 1 or 0")
    qmatrix
}

## Check misclassification matrix indicators

msm.check.ematrix <- function(ematrix, qmatrix)
{
    ematrix <- as.matrix(ematrix)
    if (!all(dim(qmatrix) == dim(ematrix)))
      stop("Dimensions of qmatrix and ematrix should be the same")
    diag(ematrix) <- 0
    if (!all ( is.element(ematrix, c(0,1)) ))
      stop("Not all off-diagonal elements of ematrix are 1 or 0")
    ematrix
}

### Check elements of state vector

msm.check.state <- function(nstates, state, fromto=FALSE, tostate=NULL)
  {
      statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
      if (fromto) {
          if (length(setdiff(unique(state), 1:nstates)) > 0)
            stop(paste("From-state vector contains elements not in",statelist))
          if (length(setdiff(unique(tostate), 1:nstates)) > 0)
            stop(paste("To-state vector contains elements not in",statelist))          
      }
      else if (length(setdiff(unique(state), 1:nstates)) > 0)
        stop(paste("State vector contains elements not in",statelist))
      invisible()
  }

### CHECK IF TRANSITION PROBABILITIES FOR DATA ARE ALL NON-ZERO
### (e.g. check for backwards transitions when the model is irreversible)

msm.check.model <- function(state, subject, qmatrix, fromto=FALSE, tostate=NULL)
{
    n <- length(state)
    diag(qmatrix) <- 0
    diag(qmatrix) <- - apply(qmatrix, 1, sum)
    Pmat <- MatrixExp(qmatrix)
    Pmat[Pmat < 1e-16] <- 0
    if (fromto) fromstate <- state
    else {
        fromstate <- c(NA, state[1:(n-1)])
        subj.num <- as.numeric(factor(subject))
        fromstate[subj.num != c(NA, subj.num[1:(n-1)])] <- NA
        tostate <- state
    }
    tostate <- tostate[!is.na(fromstate)]
    fromstate <- fromstate[!is.na(fromstate)]
    unitprob <- apply(cbind(fromstate, tostate), 1, function(x) { Pmat[x[1], x[2]] } )
    if (identical(all.equal(min(unitprob, na.rm=TRUE), 0),  TRUE))
      {
          badobs <- min ( (1:n) [unitprob==0], na.rm = TRUE)
          stop (paste ("Data inconsistent with transition matrix for model without misclassification:\n",
                       "individual", if(fromto) "" else subject[badobs], "moves from state", fromstate[badobs],
                       "to state", tostate[badobs], "at observation", badobs, "\n") )
      }
    invisible()
}

msm.check.times <- function(time, subject)
  {
### CHECK IF OBSERVATIONS ARE ORDERED IN TIME WITHIN SUBJECT
      orderedpt <- tapply(time, subject, function(x) { all(!duplicated(x) & order(x)==1:length(x)) })
      if (any (!orderedpt)) {
          badsubjs <- sort(unique(subject))[ !orderedpt ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          stop (paste ("Observations within subject", plural, " ", badlist, " are not in strictly increasing order of time", sep="") )
      }
### CHECK IF ANY INDIVIDUALS HAVE ONLY ONE OBSERVATION
      nobspt <- tapply(subject, subject, length)
      if (any (nobspt == 1)) {
          badsubjs <- sort(unique(subject))[ nobspt == 1 ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          has <-  if (length(badsubjs)==1) "has" else "have"
          stop (paste ("Subject", plural, " ", badlist, " only ", has, " one observation", sep="") )
      }
### CHECK IF OBSERVATIONS WITHIN A SUBJECT ARE ADJACENT
      ind <- tapply(1:length(subject), subject, length)
      imin <- tapply(1:length(subject), subject, min)
      imax <- tapply(1:length(subject), subject, max)
      adjacent <- (ind == imax-imin+1)
      if (any (!adjacent)) {  
          badsubjs <- sort(unique(subject))[ !adjacent ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          stop (paste ("Observations within subject", plural, " ", badlist, " are not adjacent in the data file", sep="") )
      }
      invisible()
  }

msm.check.tunit <- function(fromto, time, subject, timelag, tunit)
{
### CHECK SETTING OF tunit, the unit in days of observed time variable
### IF death time known to within a day, then dt should not be less than 1 / tunit
    n <- length(time)
    if (!fromto) {
        prevtime <- c(NA, time[1:(n-1)])
        subject <- as.numeric(factor(subject))
        prevtime[subject != c(NA, subject[1:(n-1)])] <- NA
        timelag <- time - prevtime
    }
    eps <- 1e-05 ## fuzz factor for approximate decimals
    if (any (timelag - 1/tunit + eps<= 0, na.rm=TRUE)) {
        badobs <- min ( (1:n) [ (timelag - 1/tunit + eps < 0) ], na.rm = TRUE)
        stop (paste ("tunit =", tunit, "is too small: time difference between observations", badobs-1,
                     "and", badobs, "is less than 1 /", tunit) )
    }
    invisible()
}

## Table of 'transitions': previous state versus current state

statetable.msm <- function(state, subject, data=NULL)
{
    if(!is.null(data)) {
        data <- as.data.frame(data)
        state <- eval(substitute(state), data, parent.frame())
    }
    n <- length(state)
    if (!is.null(data))
      subject <-
        if(missing(subject)) rep(1,n) else eval(substitute(subject), data, parent.frame())
    subject <- as.numeric(factor(subject))
    prevsubj <- c(NA, subject[1:(n-1)])
    previous <- c(NA, state[1:(n-1)])
    previous[prevsubj!=subject] <- NA
    ntrans <- table(previous, state)
    names(dimnames(ntrans)) <- c("from", "to")
    ntrans
}

## Calculate crude initial values for transition intensities by assuming observations represent the exact transition times

crudeinits.msm <- function(state, time, subject, qmatrix, data=NULL, fromto=FALSE,
                           fromstate=NULL, tostate=NULL, timelag=NULL, check=FALSE)
  {
      if (fromto) {
          if (missing(fromstate) || missing(tostate) || missing(timelag))
            stop("fromstate, tostate and timelag must be specified if fromto is TRUE")
          if(!is.null(data)) {
              data <- as.data.frame(data)
              state <- eval(substitute(fromstate), data, parent.frame())
              tostate <- eval(substitute(tostate), data, parent.frame())
              timelag <- eval(substitute(timelag), data, parent.frame())
          }
          tottime <- tapply(timelag, state, sum)
          ntrans <- table(state, tostate)
      }
      else {
          if (missing(state) || missing(time) || missing(subject))
            stop ("state, time, and subject must be specified if fromto is FALSE")
          if(!is.null(data)) {
              state <- eval(substitute(state), data, parent.frame())
              time <- eval(substitute(time), data, parent.frame())
          }
          n <- length(state)
          if (!is.null(data)) subject <- if(missing(subject)) rep(1,n) else eval(substitute(subject), data, parent.frame())
          subject <- as.numeric(factor(subject))
          nextsubj <- c(subject[2:n], NA)
          lastsubj <- (subject != nextsubj)
          timecontrib <- ifelse(lastsubj, NA, c(time[2:n], 0) - time)
          tottime <- tapply(timecontrib[!lastsubj], state[!lastsubj], sum) # total time spent in each state
          ntrans <- statetable.msm(state, subject, data=NULL) # table of transitions
      }
      qmatrix <- msm.check.qmatrix(qmatrix)
      msm.check.state(nrow(qmatrix), state, fromto, tostate)
      if (check) 
        msm.check.model(state, subject, qmatrix, fromto, tostate)
      nst <- nrow(qmatrix)
      estmat <- matrix(0, nst, nst)
      rownames(estmat) <- colnames(estmat) <- paste(1:nst)
      tab <- sweep(ntrans, 1, tottime, "/") 
      for (i in 1:nst) # Include zero rows for states for which there were no transitions
        for (j in 1:nst)
          if ((paste(i) %in% rownames(tab)) && (paste(j) %in% colnames(tab)))
            estmat[paste(i), paste(j)] <- tab[paste(i),paste(j)]
      inits <- t(estmat)[t(qmatrix)==1]
      names(inits) <- paste(t(row(qmatrix)), t(col(qmatrix)), sep="-")[t(qmatrix)==1]
      inits
  }

### Force dynamic loading of the C code library (msm.so)
.First.lib <- function(lib, pkg) library.dynam( "msm", pkg, lib )

