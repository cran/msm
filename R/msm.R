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
      qmatrix <- as.matrix(qmatrix)
      nstates <- dim(qmatrix)[1]
      if (nstates != dim(qmatrix)[2])
        stop("Number of rows and columns of qmatrix should be equal")
      diag(qmatrix) <- rep(0, nstates)
      if (!all ( is.element(qmatrix, c(0,1)) ) )
        stop("Not all off-diagonal elements of qmatrix are 1 or 0")    
      qvector <- as.numeric(t(qmatrix)) # matrix to vector by filling rows first
      nintens <- sum(qmatrix)

### MISCLASSIFICATION MATRIX
      if (misc) {
          if (missing(ematrix)) stop("Misclassification matrix not given")
          ematrix <- as.matrix(ematrix)
          if (!all(dim(qmatrix) == dim(ematrix)))
            stop("Dimensions of qmatrix and ematrix should be the same")
          diag(ematrix) <- rep(0, nstates)
          if (!all ( is.element(ematrix, c(0,1)) ))
            stop("Not all off-diagonal elements of ematrix are 1 or 0")
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
          subject <- na.omit(subject)
          fromstate <- tostate <- timelag <- NULL
      }
      npts <- length(unique(subject))
      
      covrows.kept <- misccovrows.kept <- 1:nobs
### COVARIATES
      if (!is.null(covariates)) {
          pc <- msm.process.covs(covariates, data, constraint, nobs, nintens)
          covvec <- pc$covvec;  constrvec <- pc$constrvec;  covlabels <- pc$covlabels
          ncovs <- pc$ncovs;  ncoveffs <- pc$ncoveffs; covrows.kept <- pc$kept.rows
          covmeans <- pc$covmeans
      }
      else {
          ncovs <- ncoveffs <- 0
          covvec <- constrvec <- covlabels <- NULL
      }      

### MISCLASSIFICATION COVARIATES
      if (!is.null(misccovariates) & (misc)) {
          pc <- msm.process.covs(misccovariates, data, miscconstraint, nobs, nmiscs)
          misccovvec <- pc$covvec;  miscconstrvec <- pc$constrvec;  misccovlabels <- pc$covlabels
          nmisccovs <- pc$ncovs;  nmisccoveffs <- pc$ncoveffs; misccovrows.kept <- pc$kept.rows
          misccovmeans <- pc$covmeans
      }
      else if (is.null (misccovariates) | (!misc)) {
          if (!is.null (misccovariates))
            warning("misccovariates have been specified, but misc is FALSE. Ignoring misccovariates.")
          nmisccovs <- nmisccoveffs <- 0
          misccovvec <- miscconstrvec <- misccovlabels <- NULL
      }

### DROP MISSING DATA
      
      final.rows <- intersect(intersect(statetimerows.kept, subjrows.kept),
                              intersect(covrows.kept, misccovrows.kept))
      subject <- subject[ subjrows.kept %in% final.rows ]
      time <- time[ statetimerows.kept %in% final.rows ]
      state <- state[ statetimerows.kept %in% final.rows ]
      tostate <- tostate[ statetimerows.kept %in% final.rows ]
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
      if (length(fixedpars)==length(inits))
        fixed <- TRUE
      else {
          inits <- inits[notfixed]
          fixed <- FALSE
      }

### CALCULATE LIKELIHOOD AT INITIAL VALUES...
      if (fixed) {
          likval <- lik.msm(inits, allinits, misc, subject, time, state, tostate, fromto, qvector, evector, covvec, constrvec, misccovvec, miscconstrvec, 
                            initprobs, nstates, nintens, nmiscs, nobs, npts, ncovs, ncoveffs, nmisccovs, nmisccoveffs, covmatch,
                            death, tunit, exacttimes, fixedpars, plabs)
###          bits <- likval$likbits
###          likval <- likval$endlik
          params <- inits
          covmat <- ses <- NULL
          foundstderr <- FALSE
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
              ses <- rep(0, npars)                                  
              for (i in 1:npars) 
                if (!is.element(i, fixedpars)){
                    if (plabs[i]=="qbase") ## recover SEs of intensities on non-log scale
                      ses[i] <- deltamethod(~ exp(x1), params[i], diag(covmat)[i])
                    else if (plabs[i]=="ebase") ## recover SEs of probabilities on non-logit scale
                      ses[i] <- deltamethod(~ exp(x1)/(1 + exp(x1)), params[i], diag(covmat)[i])
                }
              ses[plabs=="qcov"] <- sqrt(diag(covmat)[plabs=="qcov"])
              ses[plabs=="ecov"] <- sqrt(diag(covmat)[plabs=="ecov"])
              foundstderr <- TRUE
          }
          else {
              covmat <- ses <- "Non-positive definite approximate variance-covariance"
              foundstderr <- FALSE
          }
      }
      params[plabs=="qbase"] <- exp(params[plabs=="qbase"])
      params[plabs=="ebase"] <- expit(params[plabs=="ebase"])

### Rearrange the optimised intensities, misc probs and covariate effects into lists of matrices
      output <- msm.form.output(qmatrix, nstates, nintens,
                                npars, ncovs, constrvec, fixedpars, fixed, covlabels,
                                params, covmat, ses, foundstderr, 0)
      Qmatrices <- output$Matrices
      QmatricesSE <- if (fixed) NULL else output$MatricesSE 
      diag(Qmatrices[["baseline"]]) <- - apply(Qmatrices[["baseline"]], 1, sum)

      if (misc) {
          output <- msm.form.output(ematrix, nstates, nmiscs,
                                    npars, nmisccovs, miscconstrvec, fixedpars, fixed, misccovlabels,
                                    params, covmat, ses, foundstderr,
                                    nintens+ncoveffs)
          Ematrices <- output$Matrices
          EmatricesSE <- if (fixed) NULL else output$MatricesSE
          diag(Ematrices[["baseline"]]) <- 1 - apply(Ematrices[["baseline"]], 1, sum)
      }
      else {
          Ematrices <- EmatricesSE <- "Model without misclassfication"
      }
      ## calculate mean sojourn times with centered covariates
      sojourn <- sojourn.msm(Qmatrices, covmat, foundstderr)      
      minus2loglik <- if (fixed) likval else opt$value
      
      if (ncovs > 0) { 
          ## calculate intensity matrix with covariates set to zero
          qmean <- Qmatrices[[1]]; diag(qmean) <- 0
          logqcenter <- log(qmean)
          for (i in 1:ncovs)
            logqcenter <- logqcenter - Qmatrices[[i+1]] * covmeans[i]
          qcenter <- exp(logqcenter) 
      }
      else qcenter <- NULL
      
      if (nmisccovs > 0) {
          ## calculate misc matrix with covariates set to zero
          logitecenter <- logit(Ematrices[[1]])
          for (i in 1:nmisccovs)
            logitecenter <- logitecenter - Ematrices[[misccovlabels[i]]] * misccovmeans[i]
          ecenter <- expit(logitecenter)
      }
      else ecenter <- NULL
      
      returned <- list (
                        misc = misc,
                        Qmatrices = Qmatrices,
                        QmatricesSE = QmatricesSE,
                        qcenter = qcenter,
                        Ematrices = Ematrices,
                        EmatricesSE = EmatricesSE,
                        ecenter = ecenter, 
                        sojourn = sojourn,
                        minus2loglik = minus2loglik, 
                        Pmatrix = function(t) {MatrixExp(Qmatrices[[1]]*t)},
                        estimates = params,
                        covmat = covmat,
                        data = list(nobs=nobs, npts=npts, state=state, time=time, subject=subject,
                          fromto=fromto, tostate=tostate,
                          covvec=covvec, misccovvec=misccovvec,
                          covlabels=covlabels, misccovlabels=misccovlabels,
                          ncovs=ncovs, nmisccovs=nmisccovs, tunit=tunit),
                        model = list(qvector=qvector, evector=evector,
                          nstates=nstates, nintens=nintens, nmiscs=nmiscs, 
                          constrvec=constrvec, miscconstrvec=miscconstrvec,
                          ncoveffs=ncoveffs, nmisccoveffs=nmisccoveffs, 
                          covmatch=covmatch, initprobs=initprobs, death=death, 
                          exacttimes=exacttimes),
                        )
      
      attr(returned, "fixed") <- fixed
      class(returned) <- "msm"
      returned
  }

# msmmisc <- msm # don't retain for backwards compatibility

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
                endlik = double(1),
                )
      
###      list(endlik=lik$endlik, likbits=lik$likbits)
      lik$endlik
  }

### methods for msm objects

print.msm <- function(x,...)
  {
      cat ("\n Multi-state Markov models in continuous time \n")
      cat (" -------------------------------------------- \n\n")

      printmatlist <- function(M, MSE, matrixtype, transformtype) {
          base <- if (is.list(M)) M[[1]] else M
          nmatrix <- dim(base)[1]
          foundse <- !is.character(MSE$baseline)
          covmessage <- if(length(M) == 1) "" else "with covariates set to their means"
          cat ("  * Matrix of", matrixtype, covmessage," \n\n")
          print (M[["baseline"]], na.print="."); cat("\n")
          cat ("    corresponding standard errors \n\n")
          ses <- MSE[["baseline"]]
          print (ses, na.print="."); cat("\n")
          if (length(M) == 1)
            cat("  * No covariates on", matrixtype, "\n\n")
          else {
              for (cov in names(M)[-1]) {
                  ests <- M[[cov]]
                  cat ("  * Linear effects on", transformtype, matrixtype, "of", cov, "\n")
                  print (ests, na.print="."); cat("\n")
                  cat ("    corresponding standard errors \n\n")
                  ses <- MSE[[cov]]
                  print (ses, na.print="."); cat("\n")
              }
          }
      }

      if (!attr(x, "fixed")){
          cat (" Maximum likelihood estimates: \n\n")
          printmatlist (x$Qmatrices, x$QmatricesSE, "transition intensities", "log")
          if (x$misc)
            printmatlist (x$Ematrices, x$EmatricesSE, "misclassification probabilities", "logit")
          covmessage <- if(length(x$Qmatrices) == 1) "" else "with covariates set to their means"
          cat ("  * Mean sojourn times in transient states", covmessage, "\n\n")
          print( x$sojourn ); cat("\n")
      }

      cat (" -2 * log-likelihood: ", x$minus2loglik, "\n")
  }

summary.msm <- function(object, # fitted model
                        prevtimes = NULL,  # times at which to compare observed and expected prevalences
                        timezero = NULL,
                        initstates = NULL,
                        hazard.scale = 1,
                        ...
                        )
  {
      cat("Calculating tables of observed and expected stage occupancy...\n")
      if (object$misc) {
          prevalences <- NULL
          prevalences <- prevalencemisc.msm(object)
      }
      else if (!object$data$fromto)
        prevalences <- prevalence.msm(object, times=prevtimes, timezero, initstates)
      else prevalences <- NULL
      if (object$data$ncovs > 0) {
          if (missing (hazard.scale))
            hazard.scale <- rep(1, length(object$data$covlabels))
          hazard <- hazard.msm(object)
      }
      else {hazard <- hazard.scale <- NULL}
      ret <- list(prevalences=prevalences,
                  hazard=hazard,
                  hazard.scale=hazard.scale)
      class(ret) <- "summary.msm"
      ret
  }

print.summary.msm <- function(x,...)
  {
      if (!is.null(x$prevalences)) {
          cat("\nObserved numbers of individuals occupying stages at each time\n\n")
          print(x$prevalences$Observed, digits=2)
          cat("\nExpected numbers of individuals occupying stages at each time\n\n")
          print(x$prevalences$Expected, digits=2)
          cat("\nObserved prevalences of stages (percentages of population at risk)\n\n")
          print(x$prevalences$"Observed percentages", digits=2)
          cat("\nExpected prevalences of stages (percentages of population at risk)\n\n")
          print(x$prevalences$"Expected percentages", digits=2)
      }
      i <- 1
      for (cov in names(x$hazard)) {
          cat ("\nTransition hazard ratios corresponding to covariate effects\n\n" )
          cat (cov, " ( unit of",x$hazard.scale[i],")\n")
          print(round(x$hazard[[cov]], 2))
          i <- i+1
      }
      invisible()
  }

plot.msm <- function(x, from=NULL, to=NULL, legend.pos=NULL, ...)
  {
      qmat <- x$Qmatrices[["baseline"]]
      nstates <- dim(qmat)[1]
      transients <- (1:nstates) [ apply(qmat, 1, min) < 0]
      absorbing <-  (1:nstates) [ apply(qmat, 1, min) == 0]
      if (is.null(to))       # find an absorbing state
        to <- max(absorbing)
      if (is.null(from))
        from <- transients
      pijt <- function(Qmatrix,i,j,t) {MatrixExp(Qmatrix,t) [i,j]}
      rg <- range(x$data$time)
      timediff <- (rg[2] - rg[1]) / 50
      times <- seq(rg[1], rg[2], timediff)
      pr <- numeric()
      for (t in times)
        pr <- c(pr, pijt(qmat, from[1], to, t))
      plot(times, 1 - pr, type="l", xlab="Time", ylab="Fitted survival probability",
           ylim=c(0,1), lty = 1, ...)
      lt <- 2
      for (st in from[-1]){ 
          pr <- numeric()
          for (t in times)
            pr <- c(pr, pijt(qmat, st, to, t))
          lines(times, 1 - pr, type="l", lty = lt)
          lt <- lt+1
      }
      if (length(legend.pos) != 2)
        legend.pos <- c(max(times) - 15*timediff, 1)
      legend(legend.pos[1], legend.pos[2], legend=paste("From stage",from), lty = 1:(lt-1))
      invisible()
  }


sojourn.msm <- function(Qmatrices,  # first matrix is intensity matrix with covariates = means. Rest are linear effect matrices for each cov.
                        covmat,     # covariance matrix of all parameters
                                        # covvecmat,  # matrix of additional covariate effects? not implemented.
                        foundstderr # whether standard errors are available
                        )
  {
      nstates <- dim(Qmatrices[[1]])[1]
      qmatrix <- Qmatrices[[1]];  qmatrix[qmatrix > 0] <- 1; diag(qmatrix) <- 0
      nintens <- sum(qmatrix)
      sojstages <- (1:nstates) [apply(qmatrix, 1, sum) > 0] # transient states
      meansoj <-  - 1 / diag(Qmatrices[["baseline"]])
      meansoj <- meansoj[ sojstages ]
      names (meansoj) <- paste("Stage",sojstages)
      if (foundstderr) {
          ## calculate approximate standard errors of the sojourn times using the delta method
          sesoj <- NULL
          indmat <- t(qmatrix);   indmat[t(qmatrix)==1] <- 1:nintens;   indmat <- t(indmat) # matrix of indices of estimate vector 
          for (i in 1:nintens) {
              if ( i %in% sojstages ) {   # calculate SE of sojourn time corresponding to each transient row of the q-matrix
                  nintrow <- sum(qmatrix[i,-i])
                  sojformula <- as.formula( paste("~ -1 / ( ", paste("exp(x", 1:nintrow, ")", sep="", collapse=" + "), " )") )
                  means <-  log ( Qmatrices[[1]][i, qmatrix[i,]==1] )
                  covs <- covmat[indmat[i,-i], indmat[i,-i]]
                  sesoj <- c(sesoj, deltamethod(sojformula,  means, covs) )
              }
          }
          names (sesoj) <- paste("Stage",sojstages)
      }
      else if (is.null(covmat))
        sesoj <- NULL
      else
        sesoj <- "Unstable Hessian at the estimates"
      list(mean=meansoj, SE=sesoj)
  }


### Table of observed and expected prevalences

prevalence.msm <- function(msm,
                           times=NULL,
                           timezero=NULL,
                           initstates=NULL
                           )
  {
      if (msm$misc) stop("Use prevalencemisc.msm for models with misclassification")
      if (msm$data$fromto) stop("Need observation time data for a prevalence table (not \"fromstate-tostate\" style)")
      time <- msm$data$time; state <- msm$data$state; subject <- msm$data$subject
      if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)
      nstates <- msm$model$nstates
      absorbing <-  (1:nstates) [ apply(msm$Qmatrices$baseline, 1, min) == 0]
      if (is.null(timezero))
        timezero <- min(msm$data$time)
      if (is.null(initstates)){
          initstates <- table(state[time==timezero])
          z <- setdiff(paste(1:nstates), names(initstates))
          zero <- rep(0, length(z)); names(zero) <- z; 
          initstates <- c(initstates, zero)
          initstates <- initstates[order(names(initstates))]
      }
      initprobs <- initstates / sum(initstates)
      Qmatrix <- if (is.list(msm$Qmatrices)) msm$Qmatrices[[1]] else msm$Qmatrices

### Estimate observed state occupancies in the data at a series of times
### Assume previous observed state is retained until next observation time
      getcontrib <- function(pt, subject, time, state, times, nstates){
          rows <- subject==pt
          contrib <- matrix(0, nrow=length(times), ncol=nstates)
          risk <- rep(0, length(times))
          for (i in seq(along=times)){
              t <- times[i]
              endtime <- max(time[rows])
              endstate <- state[rows][length(state[rows])]
              if ( t <= endtime ) {
                  risk[i] <-  1
                  for (j in 1:(length(time[rows])-1)){
                      if ( ( (time[rows][j+1] > t) & (time[rows][j] <= t)) |
                          ( (time[rows][j+1] == t) & (t == endtime) ) )
                          {
                              currstate <- state[rows][j]
                              contrib[i, currstate] <- 1
                          }
                  }
              }
              else if (t >  endtime  &  endstate %in% absorbing) {
                  risk[i] <- 1
                  contrib[i, endstate] <- 1
              }
          }
          list(risk, contrib)   
      }
      obstab <- matrix(0, nrow=length(times), ncol=nstates)
      risk <- rep(0, length(times))
      for (pt in unique(subject)){
          cont <- getcontrib(pt, subject, time, state, times, nstates)
          risk <- risk + cont[[1]]
          obstab <- obstab + cont[[2]]
      }
      obstab <- cbind(obstab, apply(obstab, 1, sum))
      dimnames(obstab) <- list(times, c(paste("Stage",1:nstates), "Total"))
      obsperc <- 100*obstab[,1:nstates] / obstab[, nstates+1]

### Work out expected state occupancies from trans probs
      exptab <- initstates
      for (j in 2 : length(times) ) {
          pmat <- MatrixExp(Qmatrix, times[j] - timezero)
          expj <- risk[j] * initprobs %*% pmat
          exptab <- rbind(exptab, expj)
      }
      exptab <- cbind(exptab, apply(exptab, 1, sum))
      dimnames(exptab) <- list(times, c(paste("Stage",1:nstates),"Total"))
      expperc <- 100*exptab[,1:nstates] / exptab[, nstates+1]

      res <- list(observed=obstab, obsperc=obsperc, expected=exptab, expperc=expperc)
      names(res) <- c("Observed", "Observed percentages", "Expected", "Expected percentages")
      res      

  }


### Estimated observed and expected values for models with misclassification 

prevalencemisc.msm <- function(msm, 
                               times = NULL,
                               observed = TRUE,
                               expected = TRUE
                               )
  {
      time <- msm$data$time; state <- msm$data$state; subject <- msm$data$subject
      if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)
      nstates <- msm$model$nstates
      nobs <- msm$data$nobs
      if (!msm$misc) stop("Use prevalence.msm for models without misclassification")
      ematrix <- t(matrix(msm$model$evector, nrow=nstates))
      esum <- apply(ematrix, 1, sum)
      nms <- max ( seq(along=esum)[esum > 0])
      
      ## Observed : For each time-interval (e.g. year) in 'times',
      ## calculate the proportion of observations in each state
      if (msm$model$death) {
          deathstring <- "and numbers of deaths"
      }
      else deathstring <- ""
      if (observed) {
          cat(paste("Calculating approximate observed state prevalences",deathstring,"in intervals ...\n"))
          ptsplit <- split(data.frame(state, time), subject)
          
          getcont <- function(state, time){
              ct <- cut(time, times, include.lowest=TRUE)
              tct <- split(state, ct)
              contribs <- sapply(tct,
                                 function(x) {
                                     keep <- if (msm$model$death) x!=nstates else TRUE
                                     av <- table(x[keep]) / length(x[keep])
                                     av <- av[paste(1:nstates)]
                                     if (msm$model$death & nstates %in% x) # death time known exactly. 
                                       av[nstates] <- 1
                                     av[is.na(av)] <- 0
                                     av
                                 }
                                 )
              othertimes <- setdiff( levels(ct),  dimnames(contribs)[[2]] )
              ## Include zero contributions for time intervals in which the patient was not observed
              zerocontribs <- matrix(0, nstates, length(othertimes))
              dimnames(zerocontribs) <- list(paste(1:nstates),  othertimes)
              alltimes <- cbind(contribs, zerocontribs)
              alltimes <- alltimes[, levels(ct)]
              dimnames(alltimes)[[1]] <- paste(1:nstates)
              alltimes
          }          

          ## Calculate contributions to the overall table for individuals and then add them up
          contlist <- lapply(ptsplit, function(x){getcont(x[,1], x[,2])} )
          obstab <- matrix(0, ncol=length(times) - 1, nrow=nstates)
          for (i in contlist)
            obstab <- obstab + i
          obstab <- t(obstab)
          obstab.tmp <- cbind(obstab[,1:nms], apply(obstab[,1:nms], 1, sum))
          dimnames(obstab.tmp) <- list(levels(cut(time, times, include.lowest=TRUE)), c(paste("Stage",1:nms), "Total"))
          if (msm$model$death){
              ## Include column for observed numbers of deaths
              obstab <- cbind(obstab.tmp, obstab[,nstates])
              dimnames(obstab)[[2]][nstates+1] <- "Number of deaths"
          }
          else obstab <- obstab.tmp
          obsperc <- 100*obstab[,1:nms] / obstab[,nms+1]
      }
      else {obstab <- numeric(); obsperc <- numeric()}
      
      if (expected) {
          ## Expected state occupancy is calculated using one-step-ahead state prediction
          cat(paste("Calculating expected state prevalences as one-step-ahead predictions",deathstring,"in intervals ...\n"))
          do.what <- 3
          ## Call C routine to calculate P( obs i+1 = state r | obs 1...i ) for each observation i and state r
          explist <- .C("msmCEntry",
                        as.integer(do.what),
                        as.double (msm$estimates),
                        as.double (msm$estimates),
                        as.integer(1),
                        as.integer(length(msm$estimates)),
                        as.integer (subject),
                        as.double (time),
                        as.integer (state - 1),
                        as.integer(NULL),
                        as.integer(1),
                        as.integer(msm$model$qvector),
                        as.integer(msm$model$evector),
                        as.double (msm$data$covvec),
                        as.integer (msm$model$constrvec),
                        as.double (msm$data$misccovvec),
                        as.integer(msm$model$miscconstrvec),
                        as.double(msm$model$initprobs),
                        as.integer (nstates),
                        as.integer (nms),
                        as.integer (msm$model$nintens),
                        as.integer (msm$model$nmiscs),
                        as.integer (nobs),
                        as.integer (msm$data$npts),
                        as.integer (msm$data$ncovs),
                        as.integer (msm$model$ncoveffs),
                        as.integer (msm$data$nmisccovs),
                        as.integer (msm$model$nmisccoveffs),
                        as.integer (msm$model$covmatch),
                        as.integer (msm$model$death),
                        as.double (msm$data$tunit),
                        as.integer (msm$model$exacttimes),
                        as.integer(-1),
                        as.double (times),
                        as.integer (length(times)),
                        result = double (nms * nobs + length(times) - 1) )
          ## Matrix of one-step-ahead predicted probabilities of states. Rows are observations, columns are states
          exptab <-  matrix(explist$result[1:(nobs*nms)], nrow=nobs, ncol=nms)
          ptsplit <- split(cbind(time, as.data.frame(exptab)), subject)
          ## Average the probabilities within each prediction interval for each person (as in Satten and Longini)
          getcont <- function(probs, time){
              ct <- cut(time, times, include.lowest=TRUE)
              tct <- split(as.data.frame(probs), ct)
              contribs <- sapply(tct, function(x) apply(x, 2, mean) )
              othertimes <- setdiff( levels(ct),  dimnames(contribs)[[2]] )
              zerocontribs <- matrix(0, nms, length(othertimes))
              dimnames(zerocontribs) <- list(paste(1:nms),  othertimes)
              alltimes <- cbind(contribs, zerocontribs)
              alltimes <- alltimes[, levels(ct)]
              dimnames(alltimes)[[1]] <- paste(1:nms)
              alltimes
          }
          contlist <- lapply(ptsplit, function(x){getcont(x[,-1],x[,1])} )
          exptab <- matrix(0, ncol=length(times) - 1, nrow=nms)
          for (i in contlist)
            exptab <- exptab + i
          exptab <- t(exptab)
          exptab <- cbind(exptab, apply(exptab, 1, sum))
          dimnames(exptab) <- list(levels(cut(time, times, include.lowest=TRUE)), c(paste("Stage",1:nms), "Total"))
          if (msm$model$death) {
              ## Include column for predicted numbers of deaths
              expdeaths <- explist$result[(nobs * nms + 1):(nobs * nms + length(times) - 1)]
              exptab <- cbind(exptab, expdeaths)
              dimnames(exptab)[[2]][nstates+1] <- "Number of deaths"
          }
          expperc <- 100*exptab[,1:nms] / exptab[, nms+1]
      }
      else {exptab <- numeric(); expperc <- numeric()}
      res <- list(observed=obstab, obsperc=obsperc, expected=exptab, expperc=expperc)
      names(res) <- c("Observed", "Observed percentages", "Expected", "Expected percentages")
      res
  }

### Obtain hazard ratios from estimated effects of covariates on log-transition rates

hazard.msm <- function(msm, 
                       hazard.scale = 1)
  {
      keep <- (msm$Qmatrices[["baseline"]] > 0)
      nstates <- dim(msm$Qmatrices[["baseline"]]) [1]
      keepvec <- as.vector(keep)
      fromlabs <- rep(row.names(keep), nstates) [keepvec]
      tolabs <- rep(dimnames(keep)[[2]], rep(nstates, nstates)) [keepvec]
      if (msm$data$ncovs > 0) {
          haz.list <- list()
          if (is.matrix(msm$QmatricesSE$baseline)) {
              for (cov in msm$data$covlabels) {
                  haz.rat <- exp(hazard.scale*msm$Qmatrices[[cov]])[keepvec]
                  l95 <- exp(hazard.scale*(msm$Qmatrices[[cov]] - qnorm(0.975)*msm$QmatricesSE[[cov]]) )[keepvec]
                  u95 <- exp(hazard.scale*(msm$Qmatrices[[cov]] + qnorm(0.975)*msm$QmatricesSE[[cov]]) )[keepvec]
                  haz.tab <- cbind(haz.rat, l95, u95)
                  dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs),
                                            c("HR", "L95", "U95"))
                  haz.list[[cov]] <- haz.tab
              }
          }
          else {
              for (cov in msm$data$covlabels) {
                  haz.tab <- as.matrix(exp(hazard.scale*msm$Qmatrices[[cov]])[keepvec])
                  dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs), "HR")
                  haz.list[[cov]] <- haz.tab
              }
          }
      }
      else haz.list <- "No covariates on transition intensities"
      haz.list
  }


### Delta method for approximating the covariance matrix of f(X) given cov(X)

deltamethod <- function(g,       # a formula or list of formulae (functions) giving the transformation g(x) in terms of x1, x2,  etc
                        mean,    # mean, or maximum likelihood estimate, of x
                        cov,     # covariance matrix of x
                        ses=TRUE # return standard errors, else return covariance matrix
                        )
  {
      ## Var (G(x))  =  G'(mu) Var(X) G'(mu)^T
      cov <- as.matrix(cov)
      n <- length(mean)
      if (!is.list(g))
        g <- list(g)
      if ( (dim(cov)[1] != n) || (dim(cov)[2] != n) )
        stop(paste("Covariances should be a ", n, " by ", n, " matrix"))
      syms <- paste("x",1:n,sep="")
      for (i in 1:n)
        assign(syms[i], mean[i])
      gdashmu <- t(sapply(g,
                          function( form ) {
                              as.numeric(attr(eval(
                                                   ## Differentiate each formula in the list
                                                   deriv(form, syms)
                                                   ## evaluate the results at the mean
                                                   ), "gradient"))
                              ## and build the results row by row into a Jacobian matrix
                          }))
      new.covar <- gdashmu %*% cov %*% t(gdashmu)
      if (ses){
          new.se <- sqrt(diag(new.covar))
          new.se
      }
      else
        new.covar
  }

### Matrix exponential
### adapted from mexp in Jim Lindsey's rmutil library

MatrixExp <- function(mat, t = 1, n = 20, k = 3)
  {
      ev <- eigen(mat)
      if (any ( duplicated(ev$values) ) ) {
          ## series approximation
          mat <- mat*t / 2^k
          sum <- power <- diag(dim(mat)[2])
          for (r in 1:n) {
              power <- mat %*% power / r
              sum <- sum + power
          }
          for (i in 1:k)
            sum <- sum %*% sum
          sum
      }
      else
        ## spectral decomposition
        ev$vectors %*% diag(exp(ev$values * t)) %*% solve(ev$vectors)
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
      mm <- sweep(mm, 2, covmeans)
      ##       mm <- sweep(mm, 2, covstds, "/")
      covvec <- unlist(mm)
      if (is.null(constraint))
        constrvec <- 1:(nmatrix*ncovs)
      else {
          ## check and parse the list of constraints on covariates
          for (i in names(constraint))
            if (!(is.element(i, covlabels))){
                if (is.factor(data[, i]))
                  factor.warn <- "\n\tFor factor covariates, specify constraints using covnameCOVVALUE = c(...)"
                stop(paste("Covariate \"", i, "\" in constraint statement unknown.", factor.warn, sep=""))
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
                            params, covmat, ses, foundstderr, # ... inherited from msm
                            params.offset # starting index in full vector of parameters
                            )
  {
### REARRANGE THE OPTIMISED INTENSITIES AND COVARIATE EFFECTS
### 'Baseline' intensities are with covariates equal to their means.
      Matrices <- list()
      MatricesSE <- list()
      fixtf <- rep(FALSE, npars)
      fixtf[fixedpars] <- TRUE
      ## build intensity matrix and matrix of SEs from vector of estimates
      for (i in 0:ncovs) {
          matrixname <- if (i==0) "baseline" else covlabels[i] # name of the current output matrix.
          mat <- t(matrix) # I fill matrices by row, while R fills them by column. (uh, FIXME?)
          ## the relevant elements of the parameter vector
          parinds <-
            if (i==0)
              params.offset + 1:nmatrix
            else
              (params.offset + nmatrix + constrvec[((i-1)*nmatrix + 1)  :  (i*nmatrix)])
          mat[t(matrix)==1] <- params[parinds]
          mat <- t(mat)
          dimnames(mat) <- list(paste("Stage",1:nstates), paste("Stage",1:nstates))
          if (foundstderr && !fixed){
              intenscov <- covmat[parinds, parinds]
              intensse <- ses[parinds]
              semat <- t(matrix)
              semat[t(matrix)==1] <- intensse 
              semat <- t(semat)
              diag(semat) <- rep(NA, nstates)
              dimnames(semat)  <- list(paste("Stage",1:nstates), paste("Stage",1:nstates))
          }
          else if (!fixed){
              semat <- "Unstable Hessian at the estimates"
          }
          Matrices[[matrixname]] <- mat
          if (!fixed) MatricesSE[[matrixname]] <- semat
      }
      list(Matrices=Matrices,     # list of baseline intensities/misc probability matrix and linear effects of covariates
           MatricesSE=MatricesSE  # corresponding matrices of standard errors
           )
  }

### Force dynamic loading of the C code library (msm.so)
.First.lib <- function(lib, pkg) library.dynam( "msm", pkg, lib )

### Check the consistency of the supplied data with the specified model 

msm.check.consistency <- function(qmatrix, misc, fromto=FALSE, subject=NULL, state, tostate=NULL, time, death=FALSE, exacttimes=FALSE, tunit=NULL)
  {
      retval <- TRUE
      n <- length(state)
      nstates <- dim(qmatrix)[1]
      statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
      if (fromto) {
          if (length(setdiff(unique(state), 1:nstates)) > 0){
              stop(paste("From-state vector contains elements not in",statelist))
              retval <- FALSE
          }
          if (length(setdiff(unique(tostate), 1:nstates)) > 0){
              stop(paste("To-state vector contains elements not in",statelist))
              retval <- FALSE
          }
      }
      else if (length(setdiff(unique(state), 1:nstates)) > 0) {
          stop(paste("State vector contains elements not in",statelist))
          retval <- FALSE
      }

      if (!misc) {
### CHECK IF TRANSITION PROBABILITIES FOR DATA ARE ALL NON-ZERO
### (e.g. check for backwards transitions when the model is irreversible)
          diag(qmatrix) <- - apply(qmatrix, 1, sum)
          Pmat <- MatrixExp(qmatrix)
          Pmat[Pmat < 1e-16] <- 0
          if (fromto) {
              prevstate <- state
              state <- tostate
          }
          else {
              prevstate <- c(NA, state[1:(n-1)])
              prevstate[subject != c(NA, subject[1:(n-1)])] <- NA
          }
          unitprob <- apply(cbind(prevstate, state), 1, function(x) { Pmat[x[1], x[2]] } )
          if (min(unitprob, na.rm=TRUE) == 0) {
              badobs <- min ( (1:n) [unitprob == 0], na.rm = TRUE)
              stop (paste ("Data inconsistent with transition matrix for model without misclassification:\n",
                           "individual", if(fromto) "" else subject[badobs], "moves from state", prevstate[badobs],
                           "to state", state[badobs], "at observation", badobs, "\n") )
              retval <- FALSE
          }
      }
      
      if (!fromto) {
          retval <- msm.check.times(time, subject)
      }

      if (death & !exacttimes) {
          retval <- msm.check.tunit(fromto, time, subject, timelag, tunit)
      }
      retval
  }

msm.check.times <- function(time, subject)
  {
      retval <- TRUE
### CHECK IF OBSERVATIONS ARE ORDERED IN TIME WITHIN SUBJECT
      orderedpt <- tapply(time, subject, function(x) { any(duplicated(x)) | all(order(x)==1:length(x)) })
      if (any (!orderedpt)) {
          badsubjs <- unique(subject)[ !orderedpt ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          stop (paste ("Observations within subject", plural, " ", badlist, " are not strictly increasing", sep="") )
          retval <- FALSE
      }
### CHECK IF ANY INDIVIDUALS HAVE ONLY ONE OBSERVATION
      nobspt <- tapply(subject, subject, length)
      if (any (nobspt == 1)) {
          badsubjs <- unique(subject)[ nobspt == 1 ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          stop (paste ("Subject", plural, " ", badlist, " only have one observation", sep="") )
          retval <- FALSE
      }
      retval
  }

msm.check.tunit <- function(fromto, time, subject, timelag, tunit)
{
### CHECK SETTING OF tunit, the unit in days of observed time variable
### IF death time known to within a day, then dt should not be less than 1 / tunit
    retval <- TRUE
    n <- length(time)
    if (!fromto) {
        prevtime <- c(NA, time[1:(n-1)])
        prevtime[subject != c(NA, subject[1:(n-1)])] <- NA
        timelag <- time - prevtime
    }
    eps <- 1e-05 ## fuzz factor for approximate decimals
    if (any (timelag - 1/tunit + eps<= 0, na.rm=TRUE)) {
        badobs <- min ( (1:n) [ (timelag - 1/tunit + eps < 0) ], na.rm = TRUE)
        stop (paste ("tunit =", tunit, "is too small: time difference between observations", badobs-1,
                     "and", badobs, "is less than 1 /", tunit) )
        retval <- FALSE
    }
    retval
}

### Viterbi algorithm for reconstructing the most likely path through underlying states
### This is all done in C

viterbi.msm <- function(msm)
  {
      if (!msm$misc) stop("Viterbi algorithm is for models with misclassification")
      ematrix <- t(matrix(msm$model$evector, nrow=msm$model$nstates))
      esum <- apply(ematrix, 1, sum)
      nms <- max ( seq(along=esum)[esum > 0])
      do.what <- 2
      vit <- .C("msmCEntry",
                as.integer(do.what),
                as.double (msm$estimates),
                as.double (msm$estimates),
                as.integer(1),
                as.integer(length(msm$estimates)),
                as.integer (msm$data$subject),
                as.double (msm$data$time),
                as.integer (msm$data$state - 1),
                as.integer(NULL),
                as.integer(1),
                as.integer(msm$model$qvector),
                as.integer(msm$model$evector),
                as.double (msm$data$covvec),
                as.integer (msm$data$constrvec),
                as.double (msm$data$misccovvec),
                as.integer(msm$data$miscconstrvec),
                as.double(msm$model$initprobs),
                as.integer (msm$model$nstates),
                as.integer (nms),
                as.integer (msm$model$nintens),
                as.integer (msm$model$nmiscs),
                as.integer (msm$data$nobs),
                as.integer (msm$data$npts),
                as.integer (msm$data$ncovs),
                as.integer (msm$model$ncoveffs),
                as.integer (msm$data$nmisccovs),
                as.integer (msm$model$nmisccoveffs),
                as.integer (msm$model$covmatch),
                as.integer (msm$model$death),
                as.integer (msm$data$tunit),
                as.integer (msm$model$exacttimes),
                as.integer(-1),
                as.double(NULL),
                as.integer(NULL),
                fitted = double (msm$data$nobs)
                )
      
      fitted <- vit$fitted + 1
      data.frame(subject = msm$data$subject,
                 time = msm$data$time,
                 observed = msm$data$state,
                 fitted = fitted)
  }
