### Function to fit Markov multi-state models in continuous time
### with either arbitrary observation times or observed exact transition times
### with or without misclassification between true and underlying states

msm <- function(formula,   # formula with  observed Markov states   ~  observation times (required)
                subject = NULL, # optional, defaults to all the same if not given
                data=list(),       # data frame in which to interpret variable names
                qmatrix,    # matrix of 1s and 0s with indices of allowed transitions (diagonal is ignored) (required)
                gen.inits = FALSE, # generate initial values for transition intensities using crudeinits.msm
                ematrix = NULL,    # matrix of 1s and 0s with indices of allowed misclassfications (diagonal is ignored) (required)
                hmodel = NULL,  # list of constructors for hidden emission distributions 
                obstype = NULL, # optional, defaults to all 1 (snapshots) if not given
                covariates = NULL, # formula specifying covariates on transition rates.
                covinits = NULL,      # initial values for covariate effects
                constraint = NULL, # which intensities have covariates on them (as in Marshall et al.)
                misccovariates = NULL, # formula specifying covariates on misclassification probs
                misccovinits = NULL,      # initial values for misclassification covariate effects
                miscconstraint = NULL, # which misc probs have covariates on them
                hcovariates = NULL, # list of formula specifying covariates model for each hidden state 
                hcovinits = NULL,      # initial values for covariate effects on hidden emission distribution
                hconstraint = NULL, # constraints on hidden Markov model parameters
                qconstraint = NULL, # constraints on equality of baseline intensities
                econstraint = NULL, # constraints on equality of baseline misc probs
                initprobs = NULL,  # initial state occupancy probabilities
                death = FALSE,  # 'death' states, ie, entry time known exactly, but unknown transient state at previous instant
                exacttimes = FALSE, # TRUE is shortcut for all obstype 2. 
                censor = NULL,
                censor.states = NULL,
                cl = 0.95, # width of confidence intervals 
                fixedpars = NULL, # specify which parameters to fix. TRUE for all parameters
                center = TRUE, # center covariates at their means
                opt.method = c("optim","nlm"),
                hessian = TRUE,
                use.deriv = FALSE,
                deriv.test=FALSE,
                ... # options to optim or nlm
                )
  {            
      if (missing(formula)) stop("state ~ time formula not given")
      subject <- if (missing(subject)) NULL else eval(substitute(subject), data, parent.frame())
      obstype <- if (missing(obstype)) NULL else eval(substitute(obstype), data, parent.frame())
      if (missing(data)) data <- environment(formula)

### MODEL FOR TRANSITION INTENSITIES 
      qmodel <- msm.form.qmodel(qmatrix, qconstraint, exacttimes, gen.inits, formula, subject, data, censor, censor.states)
      
### MISCLASSIFICATION MODEL
      if (!missing(ematrix)) {
          emodel <- msm.form.emodel(ematrix, econstraint, initprobs, qmodel)
      }
      else emodel <- list(misc=FALSE, npars=0, ndpars=0)

### GENERAL HIDDEN MARKOV MODEL
      if (!missing(hmodel)) {
          msm.check.hmodel(hmodel, qmodel$nstates)
          if (!missing(hcovariates)) msm.check.hcovariates(hcovariates, qmodel)
          hmodel <- msm.form.hmodel(hmodel, hconstraint, initprobs, qmodel)
      }
      else {
          if (!missing(hcovariates)) stop("hcovariates have been specified, but no hmodel")
          hmodel <- list(hidden=FALSE, models=rep(0, qmodel$nstates), totpars=0, ncoveffs=0) # might change later if misc
      }
### CONVERT OLD STYLE MISCLASSIFICATION MODEL TO NEW GENERAL HIDDEN MARKOV MODEL
      if (emodel$misc) {
          hmodel <- msm.emodel2hmodel(emodel, qmodel)
      }
      else emodel <- list(misc=FALSE, npars=0, ndpars=0)

### DEATH STATES. Logical values allowed for backwards compatibility (TRUE means final state is death, FALSE means no death state)
      dmodel <- msm.form.dmodel(death, qmodel, hmodel)  # returns death, ndeath, 
      if (dmodel$ndeath > 0 && qmodel$exacttimes) warning("Ignoring death argument, as all states have exact entry times")
### CENSORING MODEL
      cmodel <- msm.form.cmodel(censor, censor.states, qmodel$qmatrix)

      msmdata.obs <- msm.form.data(formula, subject, obstype, covariates, data,
                                   hcovariates, misccovariates, qmodel, emodel, hmodel, cmodel, dmodel, exacttimes, center)

      if (hmodel$hidden || (cmodel$ncens > 0)) {
          msmdata <- msm.aggregate.hmmdata(msmdata.obs)
          msmdata$fromstate <- msmdata$tostate <- msmdata$timelag <- numeric(0)
      }
      else { 
          ## To speed calculation of the likelihood for the simple model (no
          ## HMM or censoring) data are aggregated by distinct fromstate,
          ## tostate, timelag, covariates combinations
          msmdata <- msm.obs.to.fromto(msmdata.obs)
          msm.check.model(msmdata$fromstate, msmdata$tostate, msmdata$obs, msmdata$subject, msmdata$obstype, qmodel$qmatrix, cmodel)
          msmdata <- msm.aggregate.data(msmdata)
          msmdata$subject <- msmdata$state <- msmdata$time <- numeric(0)
          for (i in c("subject", "time", "state")) msmdata[[i]] <- msmdata.obs[[i]]
          msmdata$cov <- msmdata.obs$covmat
      }

### MODEL FOR COVARIATES ON INTENSITIES
      qcmodel <-
        if (!is.null(covariates))
          msm.form.covmodel(msmdata$covdata, constraint, qmodel$npars, covinits)
        else {
            if (!is.null(constraint)) warning("constraint specified but no covariates")
            list(npars=0, ncovs=0, ndpars=0)
        }
### MODEL FOR COVARIATES ON MISCLASSIFICATION PROBABILITIES
      if (!emodel$misc || is.null(misccovariates))
        ecmodel <- list(npars=0, ncovs=0)
      if (!is.null(misccovariates)) {
          if (!emodel$misc) {
              warning("misccovariates have been specified, but misc is FALSE. Ignoring misccovariates.")
          }
          else {
              ecmodel <- msm.form.covmodel(msmdata$misccovdata, miscconstraint, emodel$npars, misccovinits)
              hcovariates <- msm.misccov2hcov(misccovariates, emodel)
              hcovinits <- msm.misccovinits2hcovinits(misccovinits, hcovariates, emodel, ecmodel)
          }
      }
### MODEL FOR COVARIATES ON GENERAL HIDDEN PARAMETERS
      if (!is.null(hcovariates)) {
          hmodel <- msm.form.hcmodel(hmodel, msmdata$hcovdata, hcovinits, hconstraint)
          if (emodel$misc) 
            hmodel$covconstr <- msm.form.hcovconstraint(miscconstraint, hmodel)
      }
      else if (hmodel$hidden) {
          hmodel <- c(hmodel, list(ncovs=rep(rep(0, hmodel$nstates), hmodel$npars), ncoveffs=0))
          class(hmodel) <- "hmodel"
      }
      if (hmodel$hidden && !emodel$misc) { 
          hmodel$constr <- msm.form.hconstraint(hconstraint, hmodel)
          hmodel$covconstr <- msm.form.hcovconstraint(hconstraint, hmodel)
      }
      
### FORM LIST OF INITIAL PARAMETERS, MATCHING PROVIDED INITS WITH SPECIFIED MODEL, FIXING SOME PARS IF REQD
      p <- msm.form.params(qmodel, qcmodel, emodel, hmodel, fixedpars)
      
      if (deriv.test) {
          ## Validate the new code for analytic derivatives against numeric derivatives
          likwrap <- function(x, ...){
              pars <- list(unlist(list(...)))
              do.call("lik.msm", c(pars, x))
          }
          myenv <- new.env()
          assign("x", list(msmdata, qmodel, qcmodel, cmodel, hmodel, p), env = myenv)
          for (i in 1:p$nopt)
            assign(paste("p", i, sep=""), p$inits[i], env = myenv)
          pvec <- paste("p",1:p$nopt,sep="")
          foo <- numericDeriv(as.call(lapply(as.list(c("likwrap", "x", pvec)), as.name)), pvec, myenv)
          an.d <- deriv.msm(p$inits, msmdata, qmodel, qcmodel, cmodel, hmodel, p)
          num.d <- attr(foo,"gradient")
          error <- max(abs((an.d - num.d)/num.d))
          return (list(analytic.deriv=an.d, numeric.deriv=num.d, error=error))
      }
      
### CALCULATE LIKELIHOOD AT INITIAL VALUES...
      if (p$fixed) {
          p$lik <- lik.msm(p$inits, msmdata, qmodel, qcmodel, cmodel, hmodel, p)
          p$deriv <- if (use.deriv) deriv.msm(p$inits, msmdata, qmodel, qcmodel, cmodel, hmodel, p) else NULL
          p$params.uniq <- p$allinits[!duplicated(p$constr)]
          p$params <- p$allinits[!duplicated(p$constr)][p$constr]
          p$foundse <- FALSE
          p$covmat <- NULL
      }

### ... OR DO MAXIMUM LIKELIHOOD ESTIMATION
      else {
          p$params <- p$allinits
          gr <- if (!hmodel$hidden && cmodel$ncens==0 && use.deriv) deriv.msm else NULL
          opt.method <- match.arg(opt.method)
          if (opt.method == "optim") {
              opt <- optim(p$inits, lik.msm, hessian=hessian, gr=gr, ...,# arguments to optim
                           msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel, 
                           cmodel=cmodel, hmodel=hmodel, paramdata=p)
              p$lik <- opt$value
              p$params[p$optpars] <- opt$par
          }
          else if (opt.method == "nlm") {
              nlmfn <- function(par) {
                  ret <- lik.msm(par, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel, 
                                 cmodel=cmodel, hmodel=hmodel, paramdata=p)
                  if (!is.null(gr))
                    attr(ret, "gradient") <- deriv.msm(par, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel, 
                                                       cmodel=cmodel, hmodel=hmodel, paramdata=p)
                  ret
              }
              opt <- nlm(nlmfn, p$inits, hessian=hessian, ...)
              p$lik <- opt$minimum
              p$params[p$optpars] <- opt$estimate
          }
          p$opt <- opt
          p$params.uniq <- p$params[!duplicated(p$constr)]
          p$params <- p$params[!duplicated(p$constr)][p$constr]
          if (hessian && all(eigen(opt$hessian)$values > 0)) {
              p$foundse <- TRUE
              p$covmat <- matrix(0, nrow=p$npars, ncol=p$npars)
              p$covmat[p$optpars,p$optpars] <- solve(0.5 * opt$hessian)
              p$covmat.uniq <- p$covmat[!duplicated(p$constr),!duplicated(p$constr), drop=FALSE]
              p$covmat <- p$covmat[!duplicated(p$constr),!duplicated(p$constr), drop=FALSE][p$constr,p$constr, drop=FALSE]
              p$ci <- cbind(p$params - qnorm(1 - 0.5*(1-cl))*sqrt(diag(p$covmat)),
                            p$params + qnorm(1 - 0.5*(1-cl))*sqrt(diag(p$covmat)))
              p$ci[p$fixedpars,] <- NA
          }
          else {
              p$foundse <- FALSE
              p$covmat <- p$ci <- NULL
          }
      }

      p$estimates.t <- p$params  # Calculate estimates and CIs on natural scale
      for (lab in rownames(.msm.TRANSFORMS)) {
          p$estimates.t[p$plabs==lab] <- get(.msm.TRANSFORMS[lab,"inv"])(p$params[p$plabs==lab])
          if (p$foundse)
            p$ci[p$plabs==lab] <- get(.msm.TRANSFORMS[lab,"inv"])(p$ci[p$plabs==lab, ])
      }
      
### REARRANGE THE VECTOR OF PARAMETER ESTIMATES (LOG-INTENSITIES, MISC PROBS AND
### COVARIATE EFFECTS) INTO LISTS OF MATRICES
      output <- msm.form.output("intens", qmodel, qcmodel, p)
      Qmatrices <- output$Matrices
      QmatricesSE <- if (p$fixed) NULL else output$MatricesSE 
      QmatricesL <- if (p$fixed) NULL else output$MatricesL
      QmatricesU <- if (p$fixed) NULL else output$MatricesU

      if (emodel$misc) {
          output <- msm.form.output("misc", emodel, ecmodel, p)
          Ematrices <- output$Matrices
          EmatricesSE <- if (p$fixed) NULL else output$MatricesSE
          EmatricesL <- if (p$fixed) NULL else output$MatricesL
          EmatricesU <- if (p$fixed) NULL else output$MatricesU
          names(Ematrices)[1] <- "logitbaseline"
          if (p$foundse & !p$fixed) names(EmatricesSE)[1] <- names(EmatricesL)[1] <-
            names(EmatricesU)[1] <- "logitbaseline"
      }
      else {
          Ematrices <- EmatricesSE <- EmatricesL <- EmatricesU <- NULL
      }
      if (hmodel$hidden) {
          hmodel <- msm.form.houtput(hmodel, p)
      }
      
### FORM A MSM OBJECT FROM THE RESULTS      
      msmobject <- list (
                         call = match.call(),
                         Qmatrices = Qmatrices, 
                         QmatricesSE = QmatricesSE, 
                         QmatricesL = QmatricesL,
                         QmatricesU = QmatricesU, 
                         minus2loglik = p$lik,
                         deriv = p$deriv,
                         estimates = p$params,
                         estimates.t = p$estimates.t,
                         fixedpars = p$fixedpars,
                         covmat = p$covmat,
                         ci = p$ci,
                         opt = p$opt, 
                         foundse = p$foundse,
                         data = msmdata,
                         qmodel = qmodel,
                         emodel = emodel,
                         qcmodel = qcmodel,
                         ecmodel = ecmodel,
                         hmodel = hmodel, 
                         cmodel = cmodel,
                         paramdata=p
                         )      
      attr(msmobject, "fixed") <- p$fixed
      class(msmobject) <- "msm"
      q <- qmatrix.msm(msmobject) # intensity matrix with centered covariates
      msmobject$Qmatrices$baseline <- q$estimates
      msmobject$QmatricesSE$baseline <- q$SE
      msmobject$QmatricesL$baseline <- q$L
      msmobject$QmatricesU$baseline <- q$U
      if (emodel$misc) {
          msmobject$Ematrices <- Ematrices
          msmobject$EmatricesSE <- EmatricesSE
          msmobject$EmatricesL <- EmatricesL
          msmobject$EmatricesU <- EmatricesU
          e <- ematrix.msm(msmobject) # misc matrix with centered covariates
          msmobject$Ematrices$baseline <- e$estimates
          msmobject$EmatricesSE$baseline <- e$SE
          msmobject$EmatricesL$baseline <- e$L
          msmobject$EmatricesU$baseline <- e$U
      }
      ## Calculate mean sojourn times with centered covariates
      msmobject$sojourn <- sojourn.msm(msmobject)      
      msmobject
 }

msm.check.qmatrix <- function(qmatrix)
  {
      if (!is.numeric(qmatrix) || ! is.matrix(qmatrix))
        stop("qmatrix should be a numeric matrix")
      if (nrow(qmatrix) != ncol(qmatrix))
        stop("Number of rows and columns of qmatrix should be equal")
      q2 <- qmatrix; diag(q2) <- 0
      if (any(q2 < 0))
        stop("off-diagonal entries of qmatrix should not be negative")
      invisible()
  }

msm.fixdiag.qmatrix <- function(qmatrix)
  {
      diag(qmatrix) <- 0
      diag(qmatrix) <- - rowSums(qmatrix)
      qmatrix
  }

msm.fixdiag.ematrix <- function(ematrix)
  {
      diag(ematrix) <- 0
      diag(ematrix) <- 1 - rowSums(ematrix)
      ematrix
  }

msm.form.qmodel <- function(qmatrix, qconstraint=NULL, exacttimes=FALSE, gen.inits=FALSE, formula, subject, data, censor, censor.states)
  {
### INTENSITY MATRIX (INPUT: qmatrix, qconstraint; OUTPUT: nstates, nintens, qmatrix, qvector, baseconstr, nintenseffs)
      if (gen.inits)
        qmatrix <- crudeinits.msm(formula, subject, qmatrix, data, censor, censor.states)
      msm.check.qmatrix(qmatrix)
      nstates <- dim(qmatrix)[1]
      diag(qmatrix) <- 0
      diag(qmatrix) <- - rowSums(qmatrix)
      if (is.null(rownames(qmatrix)))
        rownames(qmatrix) <- colnames(qmatrix) <- paste("State", seq(nstates))
      else if (is.null(colnames(qmatrix))) colnames(qmatrix) <- rownames(qmatrix)
      imatrix <- ifelse(qmatrix > 0, 1, 0)
      inits <- t(qmatrix)[t(imatrix)==1]
      npars <- sum(imatrix)
      if (!is.null(qconstraint)) {
          if (!is.numeric(qconstraint)) stop("qconstraint should be numeric")
          if (length(qconstraint) != npars)
            stop("baseline intensity constraint of length " ,length(qconstraint), ", should be ", npars)
          constr <- match(qconstraint, unique(qconstraint))
      }
      else 
        constr <- 1:npars
      ndpars <- max(constr)
      qmodel <- list(nstates=nstates, npars=npars, imatrix=imatrix, qmatrix=qmatrix, inits=inits,
                        constr=constr, ndpars=ndpars, exacttimes=exacttimes)
      class(qmodel) <- "msmqmodel"
      qmodel
  }

msm.check.ematrix <- function(ematrix, nstates)
  {
      if (!is.numeric(ematrix) || ! is.matrix(ematrix))
        stop("ematrix should be a numeric matrix")
      if (nrow(ematrix) != ncol(ematrix))
        stop("Number of rows and columns of ematrix should be equal")
      if (!all(dim(ematrix) == nstates))
        stop("Dimensions of qmatrix and ematrix should be the same")
      if (!all ( ematrix >= 0 | ematrix <= 1) )
        stop("Not all elements of ematrix are between 0 and 1")
      invisible()
  }

msm.form.emodel <- function(ematrix, econstraint=NULL, initprobs=NULL, qmodel)
  {
      msm.check.ematrix(ematrix, qmodel$nstates)
      diag(ematrix) <- 0
      imatrix <- ifelse(ematrix > 0 & ematrix < 1, 1, 0)
      diag(ematrix) <- 1 - rowSums(ematrix)
      if (is.null(rownames(ematrix)))
        rownames(ematrix) <- colnames(ematrix) <- paste("State", seq(qmodel$nstates))
      else if (is.null(colnames(ematrix))) colnames(ematrix) <- rownames(ematrix)
      dimnames(imatrix) <- dimnames(ematrix)
      npars <- sum(imatrix)
      nstates <- nrow(ematrix)
      inits <- t(ematrix)[t(imatrix)==1]
      if (is.null(initprobs))
        initprobs <- c(1, rep(0, qmodel$nstates-1))
      else {
          if (!is.numeric(initprobs)) stop("initprobs should be numeric")
          if (length(initprobs) != qmodel$nstates) stop("initprobs of length ", length(initprobs), ", should be ", qmodel$nstates)
          initprobs <- initprobs / sum(initprobs)
      }
      if (!is.null(econstraint)) {
          if (!is.numeric(econstraint)) stop("econstraint should be numeric")
          if (length(econstraint) != npars)
            stop("baseline misclassification constraint of length " ,length(econstraint), ", should be ", npars)
          constr <- match(econstraint, unique(econstraint))
      }
      else 
        constr <- 1:npars
      ndpars <- max(constr)
      emodel <- list(misc=TRUE, npars=npars, nstates=nstates, imatrix=imatrix, ematrix=ematrix, inits = inits,
                     constr=constr, ndpars=ndpars, initprobs=initprobs)
      class(emodel) <- "msmemodel"
      emodel
  }

### Extract data from supplied arguments, check consistency, drop missing data.
### Returns dataframe of cleaned data in observation time format
### Covariates returned in covmat

msm.form.data <- function(formula, subject=NULL, obstype=NULL, covariates=NULL, data=NULL,
                          hcovariates=NULL, misccovariates=NULL, qmodel, emodel, hmodel, cmodel, dmodel, exacttimes, center)
  {
      ## Parse the model formula of subject and time, getting missing values
      if (!inherits(formula, "formula")) stop("\"formula\" argument should be a formula")
      mf <- model.frame(formula, data=data)
      state <- mf[,1]
      if (!hmodel$hidden || emodel$misc)
        msm.check.state(qmodel$nstates, state=state, cmodel$censor)  ## replace after splitting form.hmodel
      time <- mf[,2]
      if (is.null(subject)) subject <- rep(1, nrow(mf))
      msm.check.times(time, subject)
      obstype <- msm.form.obstype(obstype, length(state), state, dmodel, exacttimes)
      droprows <- as.numeric(attr(mf, "na.action"))
      n <- length(c(state, droprows))
      statetimerows.kept <- (1:n)[! ((1:n) %in% droprows)]
      subjrows.kept <- (1:n) [!is.na(subject)]
      otrows.kept <- (1:n) [!is.na(obstype)]
      
      ## Parse covariates formula and extract data
      covdata <- misccovdata <- list(ncovs=0, covmat=numeric(0))
      if (!is.null(covariates)) {
          covdata <- msm.form.covdata(covariates, data, center)
      }
      if (!is.null(misccovariates) && emodel$misc) {
          misccovdata <- msm.form.covdata(misccovariates, data, center)
          hcovariates <- lapply(ifelse(rowSums(emodel$imatrix)>0, deparse(misccovariates), deparse(~1)), as.formula)
      }
      hcovdata <- vector(qmodel$nstates, mode="list")
      if (!is.null(hcovariates)) {
          for (i in seq(qmodel$nstates)) {
              if (!is.null(hcovariates) && !is.null(hcovariates[[i]]))
                hcovdata[[i]] <- msm.form.covdata(hcovariates[[i]], data, center)
              else hcovdata[[i]] <- list(ncovs=0)
          }
      }
      ## List of which covariates are in which model
      all.covlabels <- unique(c(covdata$covlabels, unlist(lapply(hcovdata, function(x)x$covlabels))))
      covdata$whichcov <- match(covdata$covlabels, all.covlabels)
      if(!is.null(hcovariates))
        for (i in seq(along=hcovdata))
          hcovdata[[i]]$whichcov <- match(hcovdata[[i]]$covlabels, all.covlabels)
        
      ## Drop missing data 
      final.rows <- intersect(statetimerows.kept, subjrows.kept)
      final.rows <- intersect(final.rows, otrows.kept)
      if (covdata$ncovs > 0)
        final.rows <- intersect(final.rows, covdata$covrows.kept)
      if (!is.null(hcovariates))
        for (i in seq(along=hcovariates))
          if (hcovdata[[i]]$ncovs > 0)
            final.rows <- intersect(final.rows, hcovdata[[i]]$covrows.kept)
      subject <- factor(subset(subject, subjrows.kept %in% final.rows))
      time <- subset(time, statetimerows.kept %in% final.rows)
      state <- subset(state, statetimerows.kept %in% final.rows)
      obstype <- subset(obstype, otrows.kept %in% final.rows)
      covmat <- numeric()
      if (covdata$ncovs > 0) {
          covmat <- subset(covdata$covmat, covdata$covrows.kept %in% final.rows)
          covdata$covmat <- NULL
      }
      for (i in seq(along=hcovariates)) {
          if (hcovdata[[i]]$ncovs > 0) {
              hcovdata[[i]]$covmat <- subset(hcovdata[[i]]$covmat, hcovdata[[i]]$covrows.kept %in% final.rows)
              covmat <- cbind(covmat, as.matrix(hcovdata[[i]]$covmat))
              hcovdata[[i]]$covmat <- NULL
          }
      }
      if (length(all.covlabels) > 0)
        covmat <- as.data.frame(covmat, optional=TRUE)[all.covlabels]
      nobs <- length(final.rows)
      nmiss <- n - nobs
      plural <- if (nmiss==1) "" else "s"
      if (nmiss > 0) warning(nmiss, " record", plural, " dropped due to missing values")
      dat <- list(state=state, time=time, subject=subject, obstype=obstype, nobs=nobs, npts=length(unique(subject)), 
                  ncovs=length(all.covlabels), covlabels=all.covlabels,
                  covdata=covdata, misccovdata=misccovdata, hcovdata=hcovdata, covmat=covmat)
#      dat <- c(dat, covmat)
      class(dat) <- "msmdata"
      dat
  }

### Check elements of state vector. For simple models and misc models specified with ematrix 
### No check is performed for hidden models

msm.check.state <- function(nstates, state=NULL, censor)
  {
      statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
      states <- c(1:nstates, censor)
      if (!is.null(state)) {
          if (length(setdiff(unique(state), states)) > 0)
            stop("State vector contains elements not in ",statelist)
          miss.state <- setdiff(states, unique(state))
          if (length(miss.state) > 0)
            warning("State vector doesn't contain observations of ",paste(miss.state, collapse=","))
      }
      invisible()
  }

msm.check.times <- function(time, subject)
  {
### Check if any individuals have only one observation
      nobspt <- table(subject)
      if (any (nobspt == 1)) {
          badsubjs <- sort(unique(subject))[ nobspt == 1 ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          has <-  if (length(badsubjs)==1) "has" else "have"
          warning ("Subject", plural, " ", badlist, " only ", has, " one observation")
      }
### Check if observations within a subject are adjacent
      ind <- tapply(1:length(subject), subject, length)
      imin <- tapply(1:length(subject), subject, min)
      imax <- tapply(1:length(subject), subject, max)
      adjacent <- (ind == imax-imin+1)
      if (any (!adjacent)) {  
          badsubjs <- sort(unique(subject))[ !adjacent ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          stop ("Observations within subject", plural, " ", badlist, " are not adjacent in the data")
      }
### Check if observations are ordered in time within subject
      orderedpt <- ! tapply(time, subject, is.unsorted)
      if (any (!orderedpt)) {
          badsubjs <- sort(unique(subject))[ !orderedpt ]
          badlist <- paste(badsubjs, collapse=", ")
          plural <- if (length(badsubjs)==1) "" else "s"
          stop ("Observations within subject", plural, " ", badlist, " are not ordered by time")
      }
      invisible()
  }

### Convert observation time data to from-to format

msm.obs.to.fromto <- function(dat)
  {
      n <- length(dat$state)
      subj.num <- as.numeric(dat$subject)
      prevsubj <- c(-Inf, subj.num[1:(n-1)])
      firstsubj <- subj.num != prevsubj
      nextsubj <- c(subj.num[2:n], Inf)
      lastsubj <- subj.num != nextsubj
      fromstate <- c(-Inf, dat$state[1:(n-1)])[!firstsubj]
      tostate <- dat$state[!firstsubj]
      timelag <- diff(dat$time)[!firstsubj[-1]]
      subject <- dat$subject[!firstsubj]
      obstype <- dat$obstype[!firstsubj]
      obs <- seq(n)[!firstsubj]
      datf <- list(fromstate=fromstate, tostate=tostate, timelag=timelag, subject=subject, obstype=obstype,
                   time=dat$time, obs=obs, firstsubj=firstsubj, npts=dat$npts, ncovs=dat$ncovs, covlabels=dat$covlabels,
                   covdata=dat$covdata, hcovdata=dat$hcovdata)
      if (datf$ncovs > 0) {
          ## match time-dependent covariates with the start of the transition 
#          datf <- c(datf, subset(as.data.frame(dat[dat$covlabels], optional=TRUE), !lastsubj))
          datf$covmat <- subset(as.data.frame(dat$covmat, optional=TRUE), !lastsubj)
      } ## n.b. don't need to  use this function for misc models 
      class(datf) <- "msmfromtodata"
      datf
  }

## Replace censored states by state with highest probability that they
## could represent. Used in msm.check.model to check consistency of
## data with transition parameters

msm.impute.censored <- function(fromstate, tostate, Pmat, cmodel)
  {
    ## e.g. cmodel$censor 99,999;  cmodel$states 1,2,1,2,3;  cmodel$index 1, 3, 6
    ## Both from and to are censored
    wb <- which ( fromstate %in% cmodel$censor & tostate %in% cmodel$censor)
    for (i in wb) {
        si <- which(cmodel$censor==fromstate[i])
        fc <- cmodel$states[(cmodel$index[si]) : (cmodel$index[si+1]-1)]
        ti <- which(cmodel$censor==tostate[i])
        tc <- cmodel$states[(cmodel$index[ti]) : (cmodel$index[ti+1]-1)]
        mp <- which.max(Pmat[fc, tc])
        fromstate[i] <- fc[row(Pmat[fc, tc])[mp]]
        tostate[i] <- tc[col(Pmat[fc, tc])[mp]]
    }
    ## Only from is censored
    wb <- which(fromstate %in% cmodel$censor)
    for (i in wb) {
        si <- which(cmodel$censor==fromstate[i])
        fc <- cmodel$states[(cmodel$index[si]) : (cmodel$index[si+1]-1)]
        fromstate[i] <- fc[which.max(Pmat[fc, tostate[i]])]
    }
    ## Only to is censored 
    wb <- which(tostate %in% cmodel$censor)
    for (i in wb) {
        si <- which(cmodel$censor==tostate[i])
        tc <- cmodel$states[(cmodel$index[si]) : (cmodel$index[si+1]-1)]
        tostate[i] <- tc[which.max(Pmat[fromstate[i], tc])]
    }
    list(fromstate=fromstate, tostate=tostate)
  }

### CHECK IF TRANSITION PROBABILITIES FOR DATA ARE ALL NON-ZERO
### (e.g. check for backwards transitions when the model is irreversible)
### obstype 1 must have unitprob > 0
### obstype 2 must have qunit != 0, and unitprob > 0. 
### obstype 3 must have unitprob > 0 

msm.check.model <- function(fromstate, tostate, obs, subject, obstype=NULL, qmatrix, cmodel)
{
    n <- length(fromstate)
    Pmat <- MatrixExp(qmatrix)
    Pmat[Pmat < 1e-16] <- 0
    imputed <- msm.impute.censored(fromstate, tostate, Pmat, cmodel)
    fs <- imputed$fromstate; ts <- imputed$tostate
    unitprob <- apply(cbind(fs, ts), 1, function(x) { Pmat[x[1], x[2]] } )
    qunit <- apply(cbind(fs, ts), 1, function(x) { qmatrix[x[1], x[2]] } )

    if (identical(all.equal(min(unitprob, na.rm=TRUE), 0),  TRUE))
      {
          badobs <- min (obs[unitprob==0], na.rm = TRUE)
          warning ("Data inconsistent with transition matrix for model without misclassification:\n",
                   "individual ", if(is.null(subject)) "" else subject[obs==badobs], " moves from state ", fromstate[obs==badobs],
                   " to state ", tostate[obs==badobs], " at observation ", badobs, "\n")
      }
    if (any(qunit[obstype==2]==0)) {
        badobs <- min (obs[qunit==0 & obstype==2], na.rm = TRUE)
        warning ("Data inconsistent with intensity matrix for observations with exact transition times and no misclassification:\n",
                 "individual ", if(is.null(subject)) "" else subject[obs==badobs], " moves from state ", fromstate[obs==badobs],
                 " to state ", tostate[obs==badobs], " at observation ", badobs)
    }
    absorbing <- absorbing.msm(qmatrix=qmatrix)
    absabs <- (fromstate %in% absorbing) & (tostate %in% absorbing)
    if (any(absabs)) {
          badobs <- min( obs[absabs] )
          warning("Absorbing - absorbing transition at observation ", badobs)
      }
    invisible()
}


## Extract covariate information from a formula.
## Find which columns and which rows to keep from the original data 
## Change in R-2.2.0, model.matrix now drops NAs.  mm returns 
 
# msm.form.covdata <- function(covariates, data, center=TRUE)
# {
#     if (!inherits(covariates, "formula")) stop(deparse(substitute(covariates)), " should be a formula")
#     mm <- as.data.frame(model.matrix(covariates, data=data, na.action=NULL)[,-1,drop=FALSE])
#     mf <- model.frame(covariates, data=data)
#     covlabels <- colnames(mm)
#     ncovs <- length(covlabels)
#     droprows <- as.numeric(attr(mf, "na.action"))
#     covrows.kept <- setdiff(seq(length=nrow(mm)), droprows)
#     mm <- mm[covrows.kept,]
#     ## centre the covariates about their means
#     covmeans <- if (center) apply(mm, 2, mean) else rep(0, ncovs)
#     covfactor <- sapply(mf, is.factor)
#     if (ncovs > 0) mm <- sweep(mm, 2, covmeans)
#     colnames(mm) <- covlabels # ( ) in names are converted into . in sweep, breaks factor covs
#     covdata <- list(covlabels=covlabels, ncovs=ncovs, covmeans=covmeans,
#                     covfactor=covfactor,
#                     covmat=mm,
#                     covrows.kept=covrows.kept)
#     class(covdata) <- "msmcovdata"
#     dput(covdata, "/tmp/cov-new.R")
#     covdata
# }

msm.form.covdata <- function(covariates, data, center=TRUE)
{
    if (!inherits(covariates, "formula")) stop(deparse(substitute(covariates)), " should be a formula")
    mm <- as.data.frame(model.matrix(covariates, data=data))
    n <- nrow(model.frame(covariates, data=data, na.action=NULL))
    covlabels <- names(mm)[-1]
    ncovs <- length(covlabels)
    mm <- subset(mm, select=-1)
    mf <- model.frame(covariates, data=data)
    droprows <- as.numeric(attr(mf, "na.action"))
    covrows.kept <- (1:n)[! ((1:n) %in% droprows)]
    ## centre the covariates about their means
    covmeans <- if (center) apply(mm, 2, mean) else rep(0, ncovs)
    covfactor <- sapply(mf, is.factor)
    if (ncovs > 0) mm <- sweep(mm, 2, covmeans)
    colnames(mm) <- covlabels # ( ) in names are converted into . in sweep, breaks factor covs
    covdata <- list(covlabels=covlabels, ncovs=ncovs, covmeans=covmeans,
                    covfactor=covfactor,
                    covmat=mm,
                    covrows.kept=covrows.kept)
    class(covdata) <- "msmcovdata"
    covdata
}

### Aggregate the data by distinct values of time lag, covariate values, from state, to state, observation type
### Result is passed to the C likelihood function (for non-hidden multi-state models)

msm.aggregate.data <- function(dat)
  {
      dat2 <- as.data.frame(dat[c("fromstate","tostate","timelag","obstype")], optional=TRUE)
      dat2$covmat <- dat$covmat
      nobsf <- length(dat2$fromstate)
      apaste <- do.call("paste", c(dat2[,c("fromstate","tostate","timelag","obstype")], dat2$covmat))
      msmdata <- dat2[!duplicated(apaste),]
      msmdata <- msmdata[order(unique(apaste)),]
      msmdata$nocc <- as.numeric(table(apaste))
      apaste2 <- msmdata[,"timelag"]
      if (dat$ncovs > 0) apaste2 <- paste(apaste2,  do.call("paste", msmdata$covmat))
      ## which unique timelag/cov combination each row of aggregated data corresponds to
      ## lik.c needs this to know when to recalculate the P matrix. 
      msmdata$whicha <- match(apaste2, unique(apaste2)) 
      msmdata <- as.list(msmdata[order(apaste2),])
      msmdata <- c(msmdata, dat[c("covdata", "hcovdata", "npts", "covlabels")])
      msmdata$nobs <- length(msmdata[[1]])
      class(msmdata) <- "msmaggdata"
      msmdata
  }

### Make indicator for which distinct from, to, timelag, covariate combination each observation corresponds to
### HMM only. This indicator is not used at the moment, but may be in the future. 

msm.aggregate.hmmdata <- function(dat)
  {
      dat2 <- msm.obs.to.fromto(dat)
      firstsubj <- dat2$firstsubj
#      dat2 <- as.data.frame(dat2[c("fromstate","tostate","timelag", dat$covlabels)], optional=TRUE)
      dat2 <- as.data.frame(c(dat2[c("fromstate","tostate","timelag")], dat2$covmat), optional=TRUE)
      apaste <- as.character(do.call("paste", dat2))
      dat$whicha <- rep(0, dat$nobs)
      dat$whicha[!firstsubj] <- match(apaste, unique(apaste))
      ## index of patient's first observation 
      dat$firstobs <- c(which(firstsubj), dat$nobs+1)
      dat
  }

### Process covariates constraints, in preparation for being passed to the likelihood optimiser
### This function is called for both sets of covariates (transition rates and the misclassification probs)

msm.form.covmodel <- function(covdata, 
                              constraint,
                              nmatrix,     # number of transition intensities / misclassification probs
                              covinits
                              )
  {
      ncovs <- covdata$ncovs
      covlabels <- covdata$covlabels
      covfactor <- covdata$covfactor
      if (is.null(constraint)) {
          constraint <- rep(list(1:nmatrix), ncovs)
          names(constraint) <- covlabels
          constr <- 1:(nmatrix*ncovs)
      }
      else {
          if (!is.list(constraint)) stop(deparse(substitute(constraint)), " should be a list")
          if (!all(sapply(constraint, is.numeric)))
            stop(deparse(substitute(constraint)), " should be a list of numeric vectors")
          if (!all(names(constraint) %in% covlabels))
            stop("covariate ", paste(setdiff(names(constraint), covlabels), collapse=", "), " in ", deparse(substitute(constraint)), " unknown")
          ## check and parse the list of constraints on covariates
          for (i in names(constraint))
            if (!(is.element(i, covlabels))){
                factor.warn <- if (covfactor[i])
                  "\n\tFor factor covariates, specify constraints using covnameCOVVALUE = c(...)"
                else ""                  
                stop("Covariate \"", i, "\" in constraint statement not in model.", factor.warn)
            }
          constr <- inits <- numeric()
          maxc <- 0
          for (i in seq(along=covlabels)){
              ## build complete vectorised list of constraints for covariates in covariates statement
              ## so. e.g. constraints = (x1=c(3,3,4,4,5), x2 = (0.1,0.2,0.3,0.4,0.4))
              ##     turns into constr = c(1,1,2,2,3,4,5,6,7,7) with seven distinct covariate effects
              if (is.element(covlabels[i], names(constraint))) {
                  if (length(constraint[[covlabels[i]]]) != nmatrix)
                    stop("\"",names(constraint)[i],"\" constraint of length ",
                         length(constraint[[covlabels[i]]]),", should be ",nmatrix)
              }
              else 
                constraint[[covlabels[i]]] <- seq(nmatrix)
              constr <- c(constr, maxc + match(constraint[[covlabels[i]]], unique(constraint[[covlabels[i]]])))
              maxc <- max(constr)
          }
      }
      inits <- numeric()
      if (!is.null(covinits)) {
          if (!is.list(covinits)) warning(deparse(substitute(covinits)), " should be a list")
          else if (!all(sapply(covinits, is.numeric)))
            warning(deparse(substitute(covinits)), " should be a list of numeric vectors")
          else if (!all(names(covinits) %in% covlabels))
            warning("covariate ", paste(setdiff(names(covinits), covlabels), collapse=", "), " in ", deparse(substitute(covinits)), " unknown")
      }
      for (i in seq(along=covlabels)) {
          if (!is.null(covinits) && is.element(covlabels[i], names(covinits))) {
              thisinit <- covinits[[covlabels[i]]]
              if (!is.numeric(thisinit)) {
                  warning("initial values for covariates should be numeric, ignoring")
                  thisinit <- rep(0, nmatrix)
              }
              if (length(thisinit) != nmatrix) {
                  warning("\"", covlabels[i], "\" initial values of length ", length(thisinit), ", should be ", nmatrix, ", ignoring")
                  thisinit <- rep(0, nmatrix)
              }
              inits <- c(inits, thisinit)
          }
          else {
              inits <- c(inits, rep(0, nmatrix))
          }
      }
      npars <- ncovs*nmatrix
      ndpars <- max(unique(constr))
      ## which covariate each distinct covariate parameter corresponds to. Used in C (FormDQCov)
      whichdcov <- rep(1:ncovs, each=nmatrix)[!duplicated(constr)]
      list(npars=npars,
           ndpars=ndpars,    # number of distinct covariate effect parameters
           ncovs=ncovs,
           constr=constr,
           whichdcov=whichdcov,        
           covlabels=covlabels,
           inits = inits,
           covmeans=covdata$covmeans
           )
  }


msm.form.dmodel <- function(death, qmodel, hmodel)
{
    nstates <- qmodel$nstates
    statelist <- if (nstates==2) "1, 2" else if (nstates==3) "1, 2, 3" else paste("1, 2, ... ,",nstates)
    if (is.logical(death) && death==TRUE)
      states <- nstates
    else if (is.logical(death) && death==FALSE)
      states <- numeric(0) ## Will be changed to -1 when passing to C
    else if (!is.numeric(death)) stop("Death states indicator must be numeric")
    else if (length(setdiff(death, 1:nstates)) > 0)
      stop("Death states indicator contains states not in ",statelist)
    else states <- death
    ndeath <- length(states)
    if (hmodel$hidden) { 
        ## Form death state info from hmmIdent parameters.
        ## Special observations in outcome data which denote death states
        ## are given as the parameter to hmmIdent()
        if (!all(hmodel$models[states] == match("identity", .msm.HMODELS)))
          stop("Death states should have the identity hidden distribution hmmIdent()")
        obs <- ifelse(hmodel$npars[states]>0, hmodel$pars[hmodel$parstate %in% states], states)
    }
    else obs <- states
    if (any (states %in% transient.msm(qmatrix=qmodel$qmatrix)))
      stop("Not all the \"death\" states are absorbing states")
    list(ndeath=ndeath, states=states, obs=obs)
}

msm.form.cmodel <- function(censor=NULL, censor.states=NULL, qmatrix)
  {
      if (is.null(censor)) {
          ncens <- 0
          if (!is.null(censor.states)) warning("censor.states supplied but censor not supplied")
      }
      else {
          if (!is.numeric(censor)) stop("censor must be numeric")
          if (any(censor %in% 1:nrow(qmatrix))) warning("some censoring indicators are the same as actual states")
          ncens <- length(censor)
          if (is.null(censor.states)) {
              if (ncens > 1) {
                  warning("more than one type of censoring given, but censor.states not supplied. Assuming only one type of censoring")
                  ncens <- 1; censor <- censor[1]
              }
              absorbing <- absorbing.msm(qmatrix=qmatrix)
              if (!length(absorbing)) {
                  warning("No absorbing state and no censor.states supplied. Ignoring censoring.")
                  ncens <- 0
              }
              else {
                  transient <- setdiff(seq(length=nrow(qmatrix)), absorbing)
                  censor.states <- transient
                  states.index <- c(1, length(censor.states)+1)
              }
          }
          else { 
              if (ncens == 1) {
                  if (!is.vector(censor.states) ||
                      (is.list(censor.states) && (length(censor.states) > 1)) )
                    stop("if one type of censoring, censor.states should be a vector, or a list with one vector element")
                  if (!is.numeric(unlist(censor.states))) stop("censor.states should be all numeric")
                  states.index <- c(1, length(unlist(censor.states))+1)
              }
              else {
                  if (!is.list(censor.states)) stop("censor.states should be a list")
                  if (length(censor.states) != ncens) stop("expected ", ncens, " elements in censor.states list, found ", length(censor.states))
                  states.index <- cumsum(c(0, lapply(censor.states, length))) + 1 
              }
              censor.states <- unlist(censor.states)
          }
      }
      if (ncens==0) censor <- censor.states <- states.index <- NULL
      ## Censoring information to be passed to C 
      list(ncens = ncens, # number of censoring states
           censor = censor, # vector of their labels in the data 
           states = censor.states, # possible true states that the censoring represents 
           index = states.index # index into censor.states for the start of each true-state set, including an extra length(censor.states)+1
           )
  }

### Observation scheme
### 1: snapshots,
### 2: exact transition times (states unchanging between observation times),
### 3: death (exact entry time but state at previous instant unknown)

msm.form.obstype <- function(obstype, nobs, state, dmodel, exacttimes)
  {
      if (!is.null(obstype)) {
          if (!is.numeric(obstype)) stop("obstype should be numeric")
          if (length(obstype) == 1) obstype <- rep(obstype, nobs)
          else if (length(obstype) != nobs) stop("obstype of length ", length(obstype), ", expected 1 or ", nobs)
          if (any(! obstype %in% 1:3)) stop("elements of obstype should be 1, 2, or 3")
      }
      else if (!is.null(exacttimes) && exacttimes)
        obstype <- rep(2, nobs)
      else { 
          obstype <- rep(1, nobs)
          if (dmodel$ndeath > 0)
            obstype[state %in% dmodel$obs] <- 3
      }
      obstype
  }

msm.form.params <- function(qmodel, qcmodel, emodel, hmodel, fixedpars)
  {
      inits <- c(qmodel$inits, qcmodel$inits)
      ni <- qmodel$npars; nc <- qcmodel$npars; 
      plabs <- c(rep("qbase",ni), rep("qcov", nc))
      nh <- sum(hmodel$npars);  nhc <- sum(hmodel$ncoveffs)
      npars <- ni + nc + nh + nhc
      inits <- c(inits, hmodel$pars, unlist(hmodel$coveffect))
      plabs <- c(plabs, hmodel$plabs, rep("hcov", nhc))
      for (lab in rownames(.msm.TRANSFORMS))        
        inits[plabs==lab] <- get(.msm.TRANSFORMS[lab,"fn"])(inits[plabs==lab])
      names(inits) <- plabs
      ## Form constraint vector for complete set of parameters 
      constr <- c(qmodel$constr, ni + qcmodel$constr, ni + nc + hmodel$constr, ni + nc + nh + hmodel$covconstr)
      constr <- match(constr, unique(constr))
      ## parameters which are always fixed and not included in user-supplied fixedpars
      auxpars <- which(plabs %in% .msm.AUXPARS)
      duppars <- which(duplicated(constr))
      naux <- length(auxpars)
      ndup <- length(duppars)
      realpars <- setdiff(seq(npars), union(auxpars, duppars))
      nrealpars <- npars - naux - ndup
      if (is.logical(fixedpars))
        fixedpars <- if (fixedpars == TRUE) seq(nrealpars) else numeric()
      if (any(! (fixedpars %in% seq(along=realpars))))
        stop ( "Elements of fixedpars should be in 1, ..., ", npars - naux - ndup)
      fixedpars <- sort(c(realpars[fixedpars], auxpars))
      notfixed <- setdiff(seq(npars), fixedpars)
      allinits <- inits
      nfix <- length(fixedpars)
      optpars <- intersect(notfixed, which(!duplicated(constr)))
      nopt <- length(optpars)
      inits <- inits[optpars]
      fixed <- (nfix + ndup == npars) # TRUE if all parameters are fixed, then no optim needed, just eval likelihood      
      names(allinits) <- plabs; names(fixedpars) <- plabs[fixedpars]; names(plabs) <- NULL
      paramdata <- list(inits=inits, plabs=plabs, allinits=allinits,
                        fixed=fixed, notfixed=notfixed, optpars=optpars,
                        fixedpars=fixedpars, constr=constr, nhcovs=hmodel$ncovs,
                        npars=npars, nfix=nfix, nopt=nopt, ndup=ndup)
      paramdata
  }

### Wrapper for the C code which evaluates the -2*log-likelihood for a Markov multi-state model with misclassification
### This is optimised by optim

likderiv.msm <- function(params, deriv=0, msmdata, qmodel, qcmodel, cmodel, hmodel, paramdata)
  {
      do.what <- deriv
      p <- paramdata
      allinits <- p$allinits;  plabs <- p$plabs
      allinits[p$optpars] <- params
      ## Untransform parameters optimized on log/logit scale
      for (lab in rownames(.msm.TRANSFORMS))        
        allinits[plabs==lab] <- get(.msm.TRANSFORMS[lab,"inv"])(allinits[plabs==lab])
      ## Replicate constrained parameters
      plabs <- plabs[!duplicated(p$constr)][p$constr] 
      allinits <- allinits[!duplicated(p$constr)][p$constr]
      ## In R, work with states / parameter indices / model indices 1, ... n. In C, work with 0, ... n-1
      msmdata$fromstate <- msmdata$fromstate - 1  
      msmdata$tostate <- msmdata$tostate - 1
      msmdata$firstobs <- msmdata$firstobs - 1
      hmodel$models <- hmodel$models - 1
      hmodel$links <- hmodel$links - 1
      qmodel$intens <- allinits[plabs=="qbase"]
      qcmodel$qcov <- allinits[plabs=="qcov"]

      lik <- .C("msmCEntry",
                as.integer(do.what),
                as.integer(as.vector(t(qmodel$imatrix))),
                as.double(qmodel$intens),
                as.double(as.vector(t(matrix(qcmodel$qcov, nrow=qmodel$npars)))),
                as.double(allinits[!(plabs %in% c("qbase", "qcov", "hcov"))]),
                as.double(allinits[plabs=="hcov"]),

                ## data for non-HMM
                as.integer(msmdata$fromstate),
                as.integer(msmdata$tostate),
                as.double(msmdata$timelag),
                as.double(unlist(msmdata$covmat)),
                as.integer(msmdata$covdata$whichcov), # this is really part of the model 
                as.integer(msmdata$nocc),
                as.integer(msmdata$whicha),
                as.integer(msmdata$obstype),

                ## data for HMM or censored
                as.integer(match(msmdata$subject, unique(msmdata$subject))), 
                as.double(msmdata$time),
                as.double(msmdata$state), # If this is a misc or censored state, this is indexed from 1. 
                as.integer(msmdata$firstobs),
                
                # HMM specification
                as.integer(hmodel$hidden),
                as.integer(hmodel$models),
                as.integer(hmodel$npars),
                as.integer(hmodel$totpars),
                as.integer(hmodel$firstpar),
                as.integer(hmodel$ncovs),
                as.integer(hmodel$whichcovh),
                as.integer(hmodel$links),                
                as.double(hmodel$initprobs),

                # various dimensions
                as.integer(qmodel$nstates),
                as.integer(qmodel$npars),
                as.integer(qmodel$ndpars),
                as.integer(qcmodel$ndpars),
                as.integer(msmdata$nobs),
                as.integer(msmdata$npts),  # HMM only
                as.integer(rep(qcmodel$ncovs, qmodel$npars)),
                
                as.integer(cmodel$ncens),
                as.integer(cmodel$censor),
                as.integer(cmodel$states),
                as.integer(cmodel$index - 1),

                ## constraints needed in C to calculate derivatives
                as.integer(qmodel$constr),
                as.integer(qcmodel$constr),
                as.integer(qcmodel$whichdcov),
                
                returned = double(if (deriv)  qmodel$ndpars + qcmodel$ndpars  else 1),
                ## so that Inf values are allowed for parameters denoting truncation points of truncated distributions
                NAOK = TRUE, 
                PACKAGE = "msm"
                )
      ## transform derivatives wrt Q to derivatives wrt log Q
      if (deriv) lik$returned[1:qmodel$ndpars] <- lik$returned[1:qmodel$ndpars]*exp(params[1:qmodel$ndpars])
      lik$returned
  }

lik.msm <- function(params, ...)
  {
      likderiv.msm(params, deriv=0, ...)
  }

deriv.msm <- function(params, ...)
  {
      likderiv.msm(params, deriv=1, ...)
  }

msm.form.output <- function(whichp, model, cmodel, p)
  {
      Matrices <- MatricesSE <- MatricesL <- MatricesU <- list()
      for (i in 0:cmodel$ncovs) {
          matrixname <- if (i==0) "logbaseline" else cmodel$covlabels[i] # name of the current output matrix.
          mat <- t(model$imatrix) # I fill matrices by row, while R fills them by column. Is this sensible...?
          if (whichp=="intens")
            parinds <- if (i==0) which(p$plabs=="qbase") else which(p$plabs=="qcov")[(i-1)*model$ndpars + 1:model$npars]
          if (whichp=="misc")
            parinds <- if (i==0) which(p$plabs=="p") else which(p$plabs=="hcov")[i + cmodel$ncovs*(1:model$npars - 1)]
          mat[t(model$imatrix)==1] <- p$params[parinds]
          mat <- t(mat)
          dimnames(mat) <- dimnames(model$imatrix)
          if (p$foundse && !p$fixed){
              intenscov <- p$covmat[parinds, parinds]
              intensse <- sqrt(diag(as.matrix(intenscov)))
              semat <- lmat <- umat <- t(model$imatrix)
              semat[t(model$imatrix)==1] <- intensse 
              lmat[t(model$imatrix)==1] <- p$ci[parinds,1]
              umat[t(model$imatrix)==1] <- p$ci[parinds,2] 
              semat <- t(semat); lmat <- t(lmat); umat <- t(umat)
              diag(semat) <- diag(lmat) <- diag(umat) <- 0
              dimnames(semat)  <- dimnames(mat)
          }
          else if (!p$fixed){
              semat <- lmat <- umat <- NULL
          }
          Matrices[[matrixname]] <- mat
          if (!p$fixed) {
              MatricesSE[[matrixname]] <- semat
              MatricesL[[matrixname]] <- lmat
              MatricesU[[matrixname]] <- umat
          }
      }
      list(Matrices=Matrices,     # list of baseline log intensities/logit misc probability matrix
                                        # and linear effects of covariates
           MatricesSE=MatricesSE,  # corresponding matrices of standard errors
           MatricesL=MatricesL,  # corresponding matrices of standard errors
           MatricesU=MatricesU  # corresponding matrices of standard errors
           )
  }

msm.form.houtput <- function(hmodel, p)
  {      
      hmodel$pars <- p$estimates.t[!(p$plabs %in% c("qbase","qcov","hcov"))]
      hmodel$coveffect <- p$estimates.t[p$plabs %in% c("hcov")]
      hmodel$fitted <- !p$fixed
      hmodel$foundse <- p$foundse
      if (hmodel$foundse) {
          hmodel$ci <- p$ci[!(p$plabs %in% c("qbase","qcov","hcov")), , drop=FALSE]
          hmodel$covci <- p$ci[p$plabs %in% c("hcov"), ]
      }
      names(hmodel$pars) <- hmodel$plabs
      ## Adjust baseline categorical probabilites.  
      anypb <- tapply(hmodel$plabs, hmodel$parstate, function(x) any(x == "pbase"))
      psum <- tapply(hmodel$pars, hmodel$parstate, function(x) sum(x[names(x) == "p"]))
      hmodel$pars[hmodel$plabs=="pbase"] <- 1 - psum[anypb]
      ## TODO - calculate SEs or CIs for these using delta method?
      if (hmodel$foundse)
        hmodel$ci[hmodel$plabs=="pbase",] <- NA 
      ## Would it be better to adjust baseline means of location parameters by taking off covariate means?
      hmodel
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
    subject <- match(subject, unique(subject))
    prevsubj <- c(NA, subject[1:(n-1)])
    previous <- c(NA, state[1:(n-1)])
    previous[prevsubj!=subject] <- NA
    ntrans <- table(previous, state)
    names(dimnames(ntrans)) <- c("from", "to")
    ntrans
}

## Calculate crude initial values for transition intensities by assuming observations represent the exact transition times

crudeinits.msm <- function(formula, subject, qmatrix, data=NULL, censor=NULL, censor.states=NULL)
  {
      cens <- msm.form.cmodel(censor, censor.states, qmatrix)
      mf <- model.frame(formula, data=data)
      state <- mf[,1]
      time <- mf[,2]
      msm.check.qmatrix(qmatrix)
      msm.check.state(nrow(qmatrix), state, cens$censor)
      n <- length(state)
      if (missing(subject)) subject <- rep(1, n)
      if (!is.null(data))
        subject <- eval(substitute(subject), as.list(data), parent.frame())
      subject <- match(subject, unique(subject))
      nocens <- (! (state %in% cens$censor) )
      state <- state[nocens]; subject <- subject[nocens]; time <- time[nocens]
      n <- length(state)
      nextsubj <- c(subject[2:n], NA)
      lastsubj <- (subject != nextsubj)
      timecontrib <- ifelse(lastsubj, NA, c(time[2:n], 0) - time)
      tottime <- tapply(timecontrib[!lastsubj], state[!lastsubj], sum) # total time spent in each state
      ntrans <- statetable.msm(state, subject, data=NULL) # table of transitions      
      nst <- nrow(qmatrix)
      estmat <- matrix(0, nst, nst)
      rownames(estmat) <- colnames(estmat) <- paste(1:nst)
      tab <- sweep(ntrans, 1, tottime, "/") 
      for (i in 1:nst) # Include zero rows for states for which there were no transitions
        for (j in 1:nst)
          if ((paste(i) %in% rownames(tab)) && (paste(j) %in% colnames(tab)))
            estmat[paste(i), paste(j)] <- tab[paste(i),paste(j)]
      estmat[qmatrix == 0] <- 0 # 
      estmat <- msm.fixdiag.qmatrix(estmat)
      estmat
  }

### Unload shared library when package is detached with unloadNamespace("msm")
.onUnload <- function(libpath) { library.dynam.unload("msm", libpath) }

