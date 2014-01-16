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
                obstrue = NULL, # for hidden Markov models, which observations represent the true state
                covariates = NULL, # formula or list of formulae specifying covariates on transition rates.
                covinits = NULL,      # initial values for covariate effects
                constraint = NULL, # which intensities have covariates on them (as in Marshall et al.)
                misccovariates = NULL, # formula or list of formulae specifying covariates on misclassification probs
                misccovinits = NULL,      # initial values for misclassification covariate effects
                miscconstraint = NULL, # which misc probs have covariates on them
                hcovariates = NULL, # list of formulae specifying covariates model for each hidden state
                hcovinits = NULL,      # initial values for covariate effects on hidden emission distribution
                hconstraint = NULL, # constraints on hidden Markov model parameters
                hranges = NULL, # range constraints for HMM parameters
                qconstraint = NULL, # constraints on equality of baseline intensities
                econstraint = NULL, # constraints on equality of baseline misc probs
                initprobs = NULL,  # initial state occupancy probabilities
                est.initprobs = FALSE, # should these be estimated, starting from the given values?
                initcovariates = NULL, # If these are specified, then assume est.initprobs=TRUE, and initprobs gives initial values with these covs set to zero.
                initcovinits = NULL,
                death = FALSE,  # 'death' states, ie, entry time known exactly, but unknown transient state at previous instant
                exacttimes = FALSE, # TRUE is shortcut for all obstype 2.
                censor = NULL,
                censor.states = NULL,
                pci = NULL,
                cl = 0.95, # width of confidence intervals in saved result components
                fixedpars = NULL, # specify which parameters to fix. TRUE for all parameters
                center = TRUE, # center covariates at their means during optimisation
                opt.method = c("optim","nlm","fisher"),
                hessian = NULL,
                use.deriv = TRUE,
                use.expm = TRUE,
                analyticp=TRUE,
                ... # options to optim or nlm
                )
{
    if (missing(formula)) stop("state ~ time formula not given")
    subject <- if (missing(subject)) NULL else eval(substitute(subject), data, parent.frame())
    obstype <- if (missing(obstype)) NULL else eval(substitute(obstype), data, parent.frame())
    obstrue <- if (missing(obstrue)) NULL else eval(substitute(obstrue), data, parent.frame())
    if (missing(data)) data <- environment(formula)

### MODEL FOR TRANSITION INTENSITIES
    if (gen.inits) {
        if (is.null(hmodel) && is.null(ematrix))
            qmatrix <- crudeinits.msm(formula, subject, qmatrix, data, censor, censor.states)
        else warning("gen.inits not supported for hidden Markov models, ignoring")
    }
    qmodel <- msm.form.qmodel(qmatrix, qconstraint, exacttimes, analyticp, use.expm)

### MISCLASSIFICATION MODEL
    if (!is.null(ematrix)) {
        emodel <- msm.form.emodel(ematrix, econstraint, initprobs, est.initprobs, qmodel)
    }
    else emodel <- list(misc=FALSE, npars=0, ndpars=0)
    if (emodel$npars==0) emodel <- list(misc=FALSE, npars=0, ndpars=0) # user supplied degenerate ematrix with no misclassification

### GENERAL HIDDEN MARKOV MODEL
    if (!is.null(hmodel)) {
        msm.check.hmodel(hmodel, qmodel$nstates)
        if (!is.null(hcovariates)) msm.check.hcovariates(hcovariates, qmodel)
        hmodel <- msm.form.hmodel(hmodel, hconstraint, initprobs, est.initprobs, qmodel)
    }
    else {
        if (!is.null(hcovariates)) stop("hcovariates have been specified, but no hmodel")
        hmodel <- list(hidden=FALSE, models=rep(0, qmodel$nstates), nipars=0, nicoveffs=0, totpars=0, ncoveffs=0) # might change later if misc
    }

### CONVERT OLD STYLE MISCLASSIFICATION MODEL TO NEW GENERAL HIDDEN MARKOV MODEL
    if (emodel$misc) {
        hmodel <- msm.emodel2hmodel(emodel, qmodel)
    }
    else emodel <- list(misc=FALSE, npars=0, ndpars=0, nipars=0, nicoveffs=0)

### EXACT DEATH TIMES. Logical values allowed for backwards compatibility (TRUE means final state has exact death time, FALSE means no states with exact death times)
    dmodel <- msm.form.dmodel(death, qmodel, hmodel)  # returns death, ndeath,
    if (dmodel$ndeath > 0 && qmodel$exacttimes) warning("Ignoring death argument, as all states have exact entry times")

### CENSORING MODEL
    cmodel <- msm.form.cmodel(censor, censor.states, qmodel$qmatrix)

### READ DATA FROM SUPPLIED FORMULAE, DROP MISSING
    msmdata.obs <- msm.form.data(formula, subject, obstype, obstrue, covariates, data,
                                 hcovariates, misccovariates, initcovariates,
                                 qmodel, emodel, hmodel, cmodel, dmodel, exacttimes, center)
### EXPAND DATA AND MODEL FOR TIME DEPENDENT INTENSITIES
    if (!is.null(pci)) {
        tdmodel <- msm.pci(pci, msmdata.obs, qmodel, cmodel, center)
        if (is.null(tdmodel)) # supplied cut points not in range of data
            pci <- NULL
        else {
            cmodel <- tdmodel$cmodel
            msmdata.obs.orig <- msmdata.obs
            names(msmdata.obs.orig)[names(msmdata.obs.orig) %in% c("covmat","covmat.orig")] <- c("cov","cov.orig") # used in bootstrap
            msmdata.obs <- tdmodel$dat
            pci <- tdmodel$tcut
        }
    }

### AGGREGATE DATA
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
        for (i in c("subject", "time", "state", "n")) msmdata[[i]] <- msmdata.obs[[i]]
        msmdata$obstype.obs <- msmdata.obs$obstype
        msmdata$firstobs <- msmdata.obs$firstobs
    }
    if (is.null(pci)) {
      msmdata.obs.orig <- NULL
      msmdata$pci.imp <- rep(0, msmdata$n)
    }
    msmdata$cov <- msmdata.obs$covmat
    msmdata$cov.orig <- msmdata.obs$covmat.orig
    msmdata$covlabels.orig <- msmdata.obs$covlabels.orig

### MODEL FOR COVARIATES ON INTENSITIES
    qcmodel <-
        if (msmdata$covdata$ncovs > 0)
            msm.form.covmodel(msmdata$covdata, constraint, qmodel$npars, covinits, msmdata)
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
            ecmodel <- msm.form.covmodel(msmdata$misccovdata, miscconstraint, emodel$npars, misccovinits, msmdata)
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
    if (!is.null(initcovariates)) {
        if (hmodel$hidden)
            hmodel <- msm.form.icmodel(hmodel, msmdata$icovdata, initcovinits)
        else warning("initprobs and initcovariates ignored for non-hidden Markov models")
    }
    else if (hmodel$hidden) {
        hmodel <- c(hmodel, list(nicovs=rep(0, hmodel$nstates-1), nicoveffs=0, cri=ecmodel$cri))
        class(hmodel) <- "hmodel"
    }
    if (hmodel$hidden && !emodel$misc) {
        hmodel$constr <- msm.form.hconstraint(hconstraint, hmodel)
        hmodel$covconstr <- msm.form.hcovconstraint(hconstraint, hmodel)
    }
    if (hmodel$hidden) hmodel$ranges <- msm.form.hranges(hranges, hmodel)
### INITIAL STATE OCCUPANCY PROBABILITIES IN HMMS
    if (hmodel$hidden) hmodel <- msm.form.initprobs(hmodel, msmdata)

### FORM LIST OF INITIAL PARAMETERS, MATCHING PROVIDED INITS WITH SPECIFIED MODEL, FIXING SOME PARS IF REQUIRED
    p <- msm.form.params(qmodel, qcmodel, emodel, hmodel, fixedpars)

### CALCULATE LIKELIHOOD AT INITIAL VALUES...
    if (is.null(hessian)) hessian <- !p$fixed
    if (p$fixed) {
        p$lik <- lik.msm(p$inits, msmdata, qmodel, qcmodel, cmodel, hmodel, p)
        p$params <- p$allinits[!duplicated(abs(p$constr))][abs(p$constr)]*sign(p$constr)
        p$params.uniq <- p$allinits[!duplicated(abs(p$constr))]
        p$opt <- list(par = p$params.uniq)
        p.unfix <- msm.unfixallparams(p)
        if (!hmodel$hidden && cmodel$ncens==0) { # Derivs only for non-hidden, no censoring
            if (all(msmdata$obstype==1)) { # Fisher info only available for panel data
                info <- information.msm(p.unfix$params.uniq, msmdata, qmodel, qcmodel, cmodel, hmodel, p.unfix)
                p$deriv <- info$deriv
                p$information <- info$info
            }
            else p$deriv <- deriv.msm(p.unfix$params.uniq, msmdata, qmodel, qcmodel, cmodel, hmodel, p.unfix)
        }
        if (hessian)
            p$opt$hessian <- optimHess(par=p.unfix$params.uniq, fn=lik.msm,
                                       gr=if (!hmodel$hidden && cmodel$ncens==0 && use.deriv) deriv.msm else NULL,
                                       msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                       cmodel=cmodel, hmodel=hmodel, paramdata=p.unfix)
        p$foundse <- FALSE
        p$covmat <- NULL
    }

### ... OR DO MAXIMUM LIKELIHOOD ESTIMATION
    else {
        p$params <- p$allinits
        gr <- if (!hmodel$hidden && cmodel$ncens==0 && use.deriv) deriv.msm else NULL
        opt.method <- match.arg(opt.method)
        optim.args <- list(...)
        if (opt.method == "optim") {
            if (is.null(optim.args$method))
                optim.args$method <- if (length(p$inits)==1) "BFGS" else "Nelder-Mead"
            optim.args <- c(optim.args, list(par=p$inits, fn=lik.msm, hessian=hessian, gr=gr,
                          msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                          cmodel=cmodel, hmodel=hmodel, paramdata=p))
            opt <- do.call("optim", optim.args)
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
        else if (opt.method == "fisher") {
            if (hmodel$hidden)
                stop("Fisher scoring not supported for hidden Markov models or censored states")
            if (cmodel$ncens > 0)
                stop("Fisher scoring not supported with censored states")
            if (any(msmdata$obstype==2))
                stop("Fisher scoring not supported with exact transition times")
            if (any(msmdata$obstype==3))
                stop("Fisher scoring not supported with exact death times")
            if (is.null(optim.args$control$reltol)) reltol <- sqrt(.Machine$double.eps)
            damp <- if (is.null(optim.args$control$damp)) 0 else optim.args$control$damp
            theta <- p$inits
            lik.old <- -lik.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                               cmodel=cmodel, hmodel=hmodel, paramdata=p)
            converged <- FALSE
            while(!converged) {
                if (!is.null(optim.args$control$trace) && optim.args$control$trace > 0)
                    cat("-2loglik=",-lik.old,", pars=",theta,"\n")
                VI <- information.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                      cmodel=cmodel, hmodel=hmodel, paramdata=p)
                V <- -VI$deriv
                Info <- VI$info + diag(damp, nrow(VI$info))
                theta <- theta + solve(Info, V)
                lik.new <- -lik.msm(theta, msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                   cmodel=cmodel, hmodel=hmodel, paramdata=p)
                if (abs(lik.old - lik.new) < reltol*(abs(lik.old) + reltol))
                    converged <- TRUE
                else lik.old <- lik.new
            }
            p$lik <- -lik.new
            p$params[p$optpars] <- theta
            opt <- list(minimum=-lik.new, estimate=theta, value=-lik.new, par=theta)
            if (hessian)
                opt$hessian <- optimHess(par=theta, fn=lik.msm,
                                         gr=if (!hmodel$hidden && cmodel$ncens==0 && use.deriv) deriv.msm else NULL,
                                         msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                         cmodel=cmodel, hmodel=hmodel, paramdata=p)
        }
        if (!hmodel$hidden && cmodel$ncens==0) {
            if (all(msmdata$obstype==1)) {
                info <- information.msm(p$params[p$optpars], msmdata, qmodel, qcmodel, cmodel, hmodel, p)
                p$deriv <- info$deriv;  p$information <- info$info
            }
            else p$deriv <- deriv.msm(p$params[p$optpars], msmdata, qmodel, qcmodel, cmodel, hmodel, p)
        }
        p$opt <- opt
        p$params <- msm.rep.constraints(p$params, p, hmodel)

        if (hessian &&
            all(!is.na(opt$hessian)) && all(!is.nan(opt$hessian)) && all(is.finite(opt$hessian)) &&
            all(eigen(opt$hessian)$values > 0))
        {
            p$foundse <- TRUE
            p$covmat <- matrix(0, nrow=p$npars, ncol=p$npars)
            p$covmat[p$optpars,p$optpars] <- solve(0.5 * opt$hessian)
            p$covmat.uniq <- p$covmat[!duplicated(abs(p$constr)),!duplicated(abs(p$constr)), drop=FALSE]
            p$covmat <- p$covmat[!duplicated(abs(p$constr)),!duplicated(abs(p$constr)), drop=FALSE][abs(p$constr),abs(p$constr), drop=FALSE]
            p$ci <- cbind(p$params - qnorm(1 - 0.5*(1-cl))*sqrt(diag(p$covmat)),
                          p$params + qnorm(1 - 0.5*(1-cl))*sqrt(diag(p$covmat)))
            p$ci[p$fixedpars,] <- NA
#            for (lab in rownames(.msm.TRANSFORMS))
#                p$ci[p$plabs==lab] <- get(.msm.TRANSFORMS[lab,"inv"])(p$ci[p$plabs==lab, ])
            for (i in 1:2)
                p$ci[,i] <- gexpit(p$ci[,i], p$ranges[,"lower",drop=FALSE], p$ranges[,"upper",drop=FALSE])
        }
        else {
            p$foundse <- FALSE
            p$covmat <- p$ci <- NULL
            if (hessian)
                warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite.")
        }
    }
    p$estimates.t <- p$params  # Calculate estimates and CIs on natural scale
    p$estimates.t <- msm.inv.transform(p$params, hmodel, p$ranges)
    ## calculate CIs for misclassification probabilities (needs multivariate transform and delta method)
    if (any(p$plabs=="p") && p$foundse){
        p.se <- p.se.msm(qmodel,emodel,hmodel,qcmodel,ecmodel,p,center, covariates = if(center) "mean" else 0)
        p$ci[p$plabs %in% c("p","pbase"),] <- as.numeric(unlist(p.se[,c("LCL","UCL")]))
    }
    ## calculate CIs for initial state probabilities in HMMs (using normal simulation method)
    if (p$foundse && any(p$plabs=="initp"))  p <- initp.ci.msm(p, cl)

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
        hmodel <- msm.form.houtput(hmodel, p, msmdata)
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
                       center = center,
                       covmat = p$covmat,
                       ci = p$ci,
                       opt = p$opt,
                       foundse = p$foundse,
                       data = msmdata,
                       data.orig = msmdata.obs.orig, # before any pci imputation, NULL if no pci
                       qmodel = qmodel,
                       emodel = emodel,
                       qcmodel = qcmodel,
                       ecmodel = ecmodel,
                       hmodel = hmodel,
                       cmodel = cmodel,
                       pci = pci,
                       paramdata=p,
                       cl=cl
                       )
    attr(msmobject, "fixed") <- p$fixed
    class(msmobject) <- "msm"
    q <- qmatrix.msm(msmobject, covariates=(if(center) "mean" else 0)) # intensity matrix on natural scale
    msmobject$Qmatrices$baseline <- q$estimates
    msmobject$QmatricesSE$baseline <- q$SE
    msmobject$QmatricesL$baseline <- q$L
    msmobject$QmatricesU$baseline <- q$U
    if (emodel$misc) {
        msmobject$Ematrices <- Ematrices
        msmobject$EmatricesSE <- EmatricesSE
        msmobject$EmatricesL <- EmatricesL
        msmobject$EmatricesU <- EmatricesU
        e <- ematrix.msm(msmobject, covariates=(if(center) "mean" else 0)) # misc matrix on natural scale
        msmobject$Ematrices$baseline <- e$estimates
        msmobject$EmatricesSE$baseline <- e$SE
        msmobject$EmatricesL$baseline <- e$L
        msmobject$EmatricesU$baseline <- e$U
    }
    ## Calculate mean sojourn times at baseline
    msmobject$sojourn <- sojourn.msm(msmobject, covariates=(if(center) "mean" else 0))
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

msm.form.qmodel <- function(qmatrix, qconstraint=NULL, exacttimes=FALSE, analyticp=TRUE, use.expm=FALSE)
{
    msm.check.qmatrix(qmatrix)
    nstates <- dim(qmatrix)[1]
    qmatrix <- msm.fixdiag.qmatrix(qmatrix)
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
    ipars <- t(imatrix)[t(lower.tri(imatrix) | upper.tri(imatrix))]
    graphid <- paste(which(ipars==1), collapse="-")
    if (analyticp && graphid %in% names(.msm.graphs[[paste(nstates)]])) {
        ## analytic P matrix is implemented for this particular intensity matrix
        iso <- .msm.graphs[[paste(nstates)]][[graphid]]$iso
        perm <- .msm.graphs[[paste(nstates)]][[graphid]]$perm
        qperm <- order(perm) # diff def in 1.2.3, indexes q matrices not vectors
    }
    else {
        iso <- 0
        perm <- qperm <- NA
    }
    qmodel <- list(nstates=nstates, iso=iso, perm=perm, qperm=qperm,
                   npars=npars, imatrix=imatrix, qmatrix=qmatrix, inits=inits,
                   constr=constr, ndpars=ndpars, exacttimes=exacttimes, expm=as.numeric(use.expm))
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

msm.form.emodel <- function(ematrix, econstraint=NULL, initprobs=NULL, est.initprobs, qmodel)
{
    msm.check.ematrix(ematrix, qmodel$nstates)
    diag(ematrix) <- 0
    imatrix <- ifelse(ematrix > 0, 1, 0)
    diag(ematrix) <- 1 - rowSums(ematrix)
    if (is.null(rownames(ematrix)))
        rownames(ematrix) <- colnames(ematrix) <- paste("State", seq(qmodel$nstates))
    else if (is.null(colnames(ematrix))) colnames(ematrix) <- rownames(ematrix)
    dimnames(imatrix) <- dimnames(ematrix)
    npars <- sum(imatrix)
    nstates <- nrow(ematrix)
    inits <- t(ematrix)[t(imatrix)==1]
    if (is.null(initprobs)) {
        initprobs <- if (est.initprobs) rep(1/qmodel$nstates, qmodel$nstates) else c(1, rep(0, qmodel$nstates-1))
    }
    else {
        if (!is.numeric(initprobs)) stop("initprobs should be numeric")
        if (is.matrix(initprobs)) {
            if (ncol(initprobs) != qmodel$nstates) stop("initprobs matrix has ", ncol(initprobs), " columns, should be number of states = ", qmodel$nstates)
            if (est.initprobs) { warning("Not estimating initial state occupancy probabilities since supplied as a matrix") }
            initprobs <- initprobs / rowSums(initprobs)
            est.initprobs <- FALSE
        }
        else {
            if (length(initprobs) != qmodel$nstates) stop("initprobs vector of length ", length(initprobs), ", should be vector of length ", qmodel$nstates, " or a matrix")
            initprobs <- initprobs / sum(initprobs)
            if (est.initprobs && any(initprobs==1)) {
                est.initprobs <- FALSE
                warning("Not estimating initial state occupancy probabilities, since some are fixed to 1")
            }
        }
    }
    nipars <- if (est.initprobs) qmodel$nstates else 0
    if (!is.null(econstraint)) {
        if (!is.numeric(econstraint)) stop("econstraint should be numeric")
        if (length(econstraint) != npars)
            stop("baseline misclassification constraint of length " ,length(econstraint), ", should be ", npars)
        constr <- match(econstraint, unique(econstraint))
    }
    else
        constr <- seq(length=npars)
    ndpars <- if(npars>0) max(constr) else 0
    emodel <- list(misc=TRUE, npars=npars, nstates=nstates, imatrix=imatrix, ematrix=ematrix, inits=inits,
                   constr=constr, ndpars=ndpars, nipars=nipars, initprobs=initprobs, est.initprobs=est.initprobs)
    class(emodel) <- "msmemodel"
    emodel
}

### Extract data from supplied arguments, check consistency, drop missing data.
### Returns dataframe of cleaned data in observation time format
### Covariates returned in covmat

msm.form.data <- function(formula, subject=NULL, obstype=NULL, obstrue=NULL, covariates=NULL, data=NULL,
                          hcovariates=NULL, misccovariates=NULL, initcovariates=NULL,
                          qmodel, emodel, hmodel, cmodel, dmodel, exacttimes, center)
{
    ## Parse the model formula of subject and time, getting missing values
    if (!inherits(formula, "formula")) stop("\"formula\" argument should be a formula")
    mf <- model.frame(formula, data=data)
    state <- mf[,1]
    if (!(hmodel$hidden || emodel$misc))
        msm.check.state(qmodel$nstates, state=state, cmodel$censor)  ## replace after splitting form.hmodel
    if (is.factor(state)) state <- as.numeric(levels(state))[state]
    time <- mf[,2]
    ## TODO for v2.0: use model.frame more intelligently for handling missing data in different formulae and variables, and length checking.
    droprows <- as.numeric(attr(mf, "na.action"))
    n <- length(c(state, droprows))
    statetimerows.kept <- (1:n)[! ((1:n) %in% droprows)]
    if (is.null(subject)) subject <- rep(1, nrow(mf))
    obstype <- msm.form.obstype(obstype, n, state, statetimerows.kept, dmodel, exacttimes)
    obstrue <- msm.form.obstrue(obstrue, nrow(mf), hmodel)
    subjrows.kept <- (1:n) [!is.na(subject)]
    otrows.kept <- statetimerows.kept[!is.na(obstype)]
    omrows.kept <- statetimerows.kept[!is.na(obstrue)]

    ## Don't drop NA covariates on the transition process at the subject's final observation, since they are not used
    lastobs <- c(which(subject[1:(n-1)] != subject[2:n]), n)
    ## Parse covariates formula and extract data
    covdata <- misccovdata <- icovdata <- list(ncovs=0, covmat=numeric(0))
    if (!is.null(covariates)) {
        if (is.list(covariates))
            covdata <- msm.form.covdata.byrate(covariates, qmodel, data, lastobs, center)
        else if (inherits(covariates, "formula"))
            covdata <- msm.form.covdata(covariates, data, lastobs, center)
        else stop(deparse(substitute(covariates)), " should be a formula or list of formulae")
    }
    if (!is.null(misccovariates) && emodel$misc) {
        misccovdata <- msm.form.covdata(misccovariates, data, NULL, center)
        hcovariates <- lapply(ifelse(rowSums(emodel$imatrix)>0, deparse(misccovariates), deparse(~1)), as.formula)
    }
    hcovdata <- vector(qmodel$nstates, mode="list")
    if (!is.null(hcovariates)) {
        for (i in seq(qmodel$nstates)) {
            if (!is.null(hcovariates) && !is.null(hcovariates[[i]]))
                hcovdata[[i]] <- msm.form.covdata(hcovariates[[i]], data, NULL, center)
            else hcovdata[[i]] <- list(ncovs=0)
        }
    }
    ## Only drop NA covariates on initprobs at initial observation
    firstobs <- c(1, which(subject[2:n] != subject[1:(n-1)]) + 1)
    if (!is.null(initcovariates)) {
        icovdata <- msm.form.covdata(initcovariates, data, !firstobs, center)
    }
    ## List of which covariates are in which model
    all.covlabels <- unique(c(covdata$covlabels, unlist(lapply(hcovdata, function(x)x$covlabels)), icovdata$covlabels)) # factors as numeric contrasts
    orig.covlabels <- unique(c(covdata$covlabels.orig, unlist(lapply(hcovdata, function(x)x$covlabels.orig)), icovdata$covlabels.orig)) # factors as single variables
    covdata$whichcov <- match(covdata$covlabels, all.covlabels)
    covdata$whichcov.orig <- match(covdata$covlabels.orig, orig.covlabels)
    if(!is.null(hcovariates)) {
        for (i in seq(along=hcovdata)){
            hcovdata[[i]]$whichcov <- match(hcovdata[[i]]$covlabels, all.covlabels)
            hcovdata[[i]]$whichcov.orig <- match(hcovdata[[i]]$covlabels.orig, orig.covlabels)
        }
    }
    if (!is.null(initcovariates))
        icovdata$whichcov <- match(icovdata$covlabels, all.covlabels)

    ## Drop missing data
    final.rows <- intersect(statetimerows.kept, subjrows.kept)
    final.rows <- intersect(final.rows, otrows.kept)
    final.rows <- intersect(final.rows, omrows.kept)
    if (covdata$ncovs > 0)
        final.rows <- intersect(final.rows, covdata$covrows.kept)
    if (!is.null(hcovariates))
        for (i in seq(along=hcovariates))
            if (hcovdata[[i]]$ncovs > 0)
                final.rows <- intersect(final.rows, hcovdata[[i]]$covrows.kept)
    if (icovdata$ncovs > 0)
        final.rows <- intersect(final.rows, icovdata$covrows.kept)

    ## Drop subjects with only one observation remaining after dropping missing data
    subj.num <- factor(match(subject,unique(subject))) # avoid problems with factor subjects with empty levels
    nobspt <- table(subj.num[final.rows])[subj.num]
    final.rows <- intersect(final.rows, which(nobspt>1))

    subject <- subject[final.rows]
    time <- subset(time, statetimerows.kept %in% final.rows)
    state <- subset(state, statetimerows.kept %in% final.rows)
    msm.check.times(time, subject, state, final.rows)
    obstype <- subset(obstype, statetimerows.kept %in% final.rows)
    obstrue <- subset(obstrue, statetimerows.kept %in% final.rows)
    covmat <- covmat.orig <- covmeans <- numeric()
    if (covdata$ncovs > 0) {
        covmat <- subset(covdata$covmat, covdata$covrows.kept %in% final.rows)
        covmat.orig <- subset(covdata$covmat.orig, covdata$covrows.kept %in% final.rows)
        covdata$covmat <- covdata$covmat.orig <- NULL
    }
    for (i in seq(along=hcovariates)) {
        if (hcovdata[[i]]$ncovs > 0) {
            hcovdata[[i]]$covmat <- subset(hcovdata[[i]]$covmat, hcovdata[[i]]$covrows.kept %in% final.rows)
            hcovdata[[i]]$covmat.orig <- subset(hcovdata[[i]]$covmat.orig, hcovdata[[i]]$covrows.kept %in% final.rows)
            covmat <- cbind(covmat, as.matrix(hcovdata[[i]]$covmat))
            covmat.orig <- cbind(covmat.orig, as.matrix(hcovdata[[i]]$covmat.orig))
            hcovdata[[i]]$covmat <- hcovdata[[i]]$covmat.orig <- NULL
        }
    }
    if (icovdata$ncovs > 0) {
        icovdata$covmat <- subset(icovdata$covmat, icovdata$covrows.kept %in% final.rows)
        icovdata$covmat.orig <- subset(icovdata$covmat.orig, icovdata$covrows.kept %in% final.rows)
        covmat <- cbind(covmat, as.matrix(icovdata$covmat))
        covmat.orig <- cbind(covmat.orig, as.matrix(icovdata$covmat.orig))
        icovdata$covmat <- icovdata$covmat.orig <- NULL
    }
    nobs <- length(final.rows)
    if (length(all.covlabels) > 0) {
        covmat <- as.data.frame(covmat, optional=TRUE)[all.covlabels]
        ## don't include last observations when centering covariates, as they
        ## don't contribute to the likelihood.
        lastobs <- c(subject[1:(nobs-1)] != subject[2:nobs], TRUE)
#        lastobs <- rep(FALSE,nobs)
        covmeans <- colMeans(covmat[!lastobs,,drop=FALSE])
        if (center) covmat <- sweep(covmat, 2, covmeans)
        covmat.orig <- as.data.frame(covmat.orig, optional=TRUE)[orig.covlabels]
    }
    nmiss <- n - nobs
    plural <- if (nmiss==1) "" else "s"
    if (nmiss > 0) warning(nmiss, " record", plural, " dropped due to missing values or subjects with only one record")
    dat <- list(state=state, time=time, subject=subject, obstype=obstype, obstrue=obstrue,
                nobs=nobs, n=nobs, npts=length(unique(subject)),
                firstobs = c(1, which(subject[2:nobs] != subject[1:(nobs-1)]) + 1, nobs+1),
                ncovs=length(all.covlabels), covlabels=all.covlabels, covlabels.orig=orig.covlabels,
                covdata=covdata, misccovdata=misccovdata, hcovdata=hcovdata, icovdata=icovdata,
                covmeans=covmeans,
                covmat=covmat, # covariates including factors as 0/1 contrasts
                covmat.orig=covmat.orig # covariates in which a factor is a single variable
                )
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

msm.check.times <- function(time, subject, state=NULL, final.rows)
{
### Check if any individuals have only one observation (after excluding missing data)
### Note this shouldn't happen after 1.2 addition to msm.form.data
    subj.num <- match(subject,unique(subject)) # avoid problems with factor subjects with empty levels
    nobspt <- table(subj.num)
    if (any (nobspt == 1)) {
        badsubjs <- unique(subject)[ nobspt == 1 ]
        andothers <- if (length(badsubjs)>3) " and others" else ""
        if (length(badsubjs)>3) badsubjs <- badsubjs[1:3]
        badlist <- paste(badsubjs, collapse=", ")
        plural <- if (length(badsubjs)==1) "" else "s"
        has <-  if (length(badsubjs)==1) "has" else "have"
        warning ("Subject", plural, " ", badlist, andothers, " only ", has, " one complete observation")
    }
### Check if observations within a subject are adjacent
    ind <- tapply(1:length(subj.num), subj.num, length)
    imin <- tapply(1:length(subj.num), subj.num, min)
    imax <- tapply(1:length(subj.num), subj.num, max)
    adjacent <- (ind == imax-imin+1)
    if (any (!adjacent)) {
        badsubjs <- unique(subject)[ !adjacent ]
        andothers <- if (length(badsubjs)>3) " and others" else ""
        if (length(badsubjs)>3) badsubjs <- badsubjs[1:3]
        badlist <- paste(badsubjs, collapse=", ")
        plural <- if (length(badsubjs)==1) "" else "s"
        stop ("Observations within subject", plural, " ", badlist, andothers, " are not adjacent in the data")
    }
### Check if observations are ordered in time within subject
    orderedpt <- ! tapply(time, subj.num, is.unsorted)
    if (any (!orderedpt)) {
        badsubjs <- unique(subject)[ !orderedpt ]
        andothers <- if (length(badsubjs)>3) " and others" else ""
        if (length(badsubjs)>3) badsubjs <- badsubjs[1:3]
        badlist <- paste(badsubjs, collapse=", ")
        plural <- if (length(badsubjs)==1) "" else "s"
        stop ("Observations within subject", plural, " ", badlist, andothers, " are not ordered by time")
    }
### Check if any consecutive observations are made at the same time, but with different states
    if (!is.null(state)){
        prevsubj <- c(-Inf, subj.num[1:length(subj.num)-1])
        prevtime <- c(-Inf, time[1:length(time)-1])
        prevstate <- c(-Inf, state[1:length(state)-1])
        sametime <- final.rows[subj.num==prevsubj & prevtime==time & prevstate!=state]
        badlist <- paste(paste(sametime-1, sametime, sep=" and "), collapse=", ")
        if (any(sametime))
            warning("Different states observed at the same time on the same subject at observations ", badlist)
    }
    invisible()
}

### Convert observation time data to from-to format

msm.obs.to.fromto <- function(dat)
{
    n <- length(dat$state)
    subj.num <- match(dat$subject, unique(dat$subject))
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
                 obstype.obs=dat$obstype, # need to keep this, e.g. for bootstrap resampling.
                 covdata=dat$covdata, hcovdata=dat$hcovdata, covmeans=dat$covmeans)
    if (datf$ncovs > 0) {
        ## match time-dependent covariates with the start of the transition
        datf$covmat <- subset(as.data.frame(dat$covmat, optional=TRUE), !lastsubj)
    } ## n.b. don't need to  use this function for misc models
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

### CHECK THAT TRANSITION PROBABILITIES FOR DATA ARE ALL NON-ZERO
### (e.g. check for backwards transitions when the model is irreversible)
### obstype 1 must have unitprob > 0
### obstype 2 must have qunit != 0, and unitprob > 0.
### obstype 3 must have unitprob > 0

msm.check.model <- function(fromstate, tostate, obs, subject, obstype=NULL, qmatrix, cmodel)
{
    n <- length(fromstate)
    qmatrix <- qmatrix / mean(qmatrix[qmatrix>0]) # rescale to avoid false warnings with small rates
    Pmat <- MatrixExp(qmatrix)
    Pmat[Pmat < 1e-16] <- 0
    imputed <- msm.impute.censored(fromstate, tostate, Pmat, cmodel)
    fs <- imputed$fromstate; ts <- imputed$tostate
    unitprob <- apply(cbind(fs, ts), 1, function(x) { Pmat[x[1], x[2]] } )
    qunit <- apply(cbind(fs, ts), 1, function(x) { qmatrix[x[1], x[2]] } )

    if (identical(all.equal(min(unitprob, na.rm=TRUE), 0),  TRUE))
    {
        badobs <- which.min(unitprob)
        warning ("Data may be inconsistent with transition matrix for model without misclassification:\n",
                 "individual ", if(is.null(subject)) "" else subject[badobs], " moves from state ", fromstate[badobs],
                 " to state ", tostate[badobs], " at non-missing observation ", obs[badobs], "\n")
    }
    if (any(qunit[obstype==2]==0)) {
        badobs <- min (obs[qunit==0 & obstype==2], na.rm = TRUE)
        warning ("Data may be inconsistent with intensity matrix for observations with exact transition times and no misclassification:\n",
                 "individual ", if(is.null(subject)) "" else subject[obs==badobs], " moves from state ", fromstate[obs==badobs],
                 " to state ", tostate[obs==badobs], " at non-missing observation ", badobs)
    }
    absorbing <- absorbing.msm(qmatrix=qmatrix)
    absabs <- (fromstate %in% absorbing) & (tostate %in% absorbing)
    if (any(absabs)) {
        badobs <- min( obs[absabs] )
        warning("Absorbing - absorbing transition at non-missing observation ", badobs)
    }
    invisible()
}


## Extract covariate information from a formula.
## Find which columns and which rows to keep from the original data

msm.form.covdata <- function(covariates, data,
                             ignore.obs=NULL, # observation IDs for which missing data should not be dropped
                             center=TRUE)
{
    if (!inherits(covariates, "formula")) stop(deparse(substitute(covariates)), " should be a formula")
    mf1 <- model.frame(covariates, data=data, na.action=NULL)
    ## We shouldn't drop NA covariates at the subject's final
    ## observation, since they are not used in the analysis, therefore
    ## we impute observed zeros when final observations are NA.
    ## TODO for v2.0: use model frames more intelligently with a new na.omit method
    mf.imp <- mf1
    for (i in names(mf.imp))
        if (!is.null(ignore.obs))
            mf.imp[ignore.obs,i][is.na(mf.imp[ignore.obs,i])] <- if (is.factor(mf.imp[,i])) levels(mf.imp[,i])[1] else 0
    mm <- na.omit(as.data.frame(model.matrix(covariates, data=mf.imp)))
    n <- nrow(mf.imp)
    covlabels <- names(mm)[-1]
    ncovs <- length(covlabels)
    mm <- subset(mm, select=-1)
    mf2 <- na.omit(mf.imp)
    droprows <- as.numeric(attr(mf2, "na.action"))
    covrows.kept <- (1:n)[! ((1:n) %in% droprows)]
    covfactor <- sapply(mf2, is.factor)
    covfactorlevels <- lapply(mf2, levels)
    colnames(mm) <- covlabels # ( ) in names are converted into . in sweep, breaks factor covs
    covdata <- list(covlabels=covlabels, ncovs=ncovs,
                    covfactor=covfactor,
                    covfactorlevels=covfactorlevels,
                    covmat=mm, # data with factors as set of 0/1 contrasts, and centered.  (for model fitting)
                    covmat.orig=mf2, # with factors kept as one variable, and not centered (for bootstrap refitting)
                    covlabels.orig=colnames(mf2),
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
    msmdata$whicha <- match(apaste2, sort(unique(apaste2)))
    msmdata <- as.list(msmdata[order(apaste2,msmdata$fromstate,msmdata$tostate),])
    ## for Fisher information: number of obs over timelag/covs starting in
    ## fromstate, replicated for all tostates.
    apaste2 <- paste(msmdata$fromstate, apaste2)
    msmdata$noccsum <- tapply(msmdata$nocc, apaste2, sum)[apaste2]
    msmdata <- c(msmdata, dat[c("covdata", "hcovdata", "npts", "covlabels","covmeans")])
    msmdata$nobs <- length(msmdata[[1]])
    ## number of disaggregated transitions
    msmdata$ntrans <- nrow(dat2)
    class(msmdata) <- "msmaggdata"
    msmdata
}

### Make indicator for which distinct from, to, timelag, covariate combination each observation corresponds to
### HMM only. This indicator is not used at the moment, but may be in the future.

msm.aggregate.hmmdata <- function(dat)
{
    dat2 <- msm.obs.to.fromto(dat)
    firstsubj <- dat2$firstsubj
    dat2 <- as.data.frame(c(dat2[c("fromstate","tostate","timelag")], dat2$covmat), optional=TRUE)
    apaste <- as.character(do.call("paste", dat2))
    dat$whicha <- rep(0, dat$nobs)
    dat$whicha[!firstsubj] <- match(apaste, unique(apaste))
    ## index of patient's first observation, plus extra element for
    ## last person's last observation plus one. This is used.
    dat$firstobs <- c(which(firstsubj), dat$nobs+1)
    dat
}

msm.check.constraint <- function(constraint, covdata){
    if (is.null(constraint)) return(invisible())
    covlabels <- covdata$covlabels
    covfactor <- covdata$covfactor
    if (!is.list(constraint)) stop(deparse(substitute(constraint)), " should be a list")
    if (!all(sapply(constraint, is.numeric)))
        stop(deparse(substitute(constraint)), " should be a list of numeric vectors")
    ## check and parse the list of constraints on covariates
    for (i in names(constraint))
        if (!(is.element(i, covlabels))){
            factor.warn <- if ((i %in% names(covfactor)) && covfactor[i])
                "\n\tFor factor covariates, specify constraints using covnameCOVVALUE = c(...)"
            else ""
            stop("Covariate \"", i, "\" in constraint statement not in model.", factor.warn)
        }
}

msm.check.covinits <- function(covinits, covdata){
    covlabels <- covdata$covlabels
    if (!is.list(covinits)) warning(deparse(substitute(covinits)), " should be a list")
    else if (!all(sapply(covinits, is.numeric)))
        warning(deparse(substitute(covinits)), " should be a list of numeric vectors")
    else if (!all(names(covinits) %in% covlabels))
        warning("covariate ", paste(setdiff(names(covinits), covlabels), collapse=", "), " in ", deparse(substitute(covinits)), " unknown")
}

### Process covariates constraints, in preparation for being passed to the likelihood optimiser
### This function is called for both sets of covariates (transition rates and the misclassification probs)

msm.form.covmodel <- function(covdata,
                              constraint,
                              nmatrix,     # number of transition intensities / misclassification probs
                              covinits,
                              fulldata
                              )
{
    if (!is.null(covdata$cri))
        return(msm.form.covmodel.byrate(covdata, constraint, nmatrix, covinits, fulldata))
    ncovs <- covdata$ncovs
    covlabels <- covdata$covlabels
    covlabels.orig <- covdata$covlabels.orig
    if (is.null(constraint)) {
        constraint <- rep(list(1:nmatrix), ncovs)
        names(constraint) <- covlabels
        constr <- 1:(nmatrix*ncovs)
    }
    else {
        msm.check.constraint(constraint, covdata)
        constr <- inits <- numeric()
        maxc <- 0
        for (i in seq(along=covlabels)){
            ## build complete vectorised list of constraints for covariates in covariates statement
            ## so. e.g. constraints = (x1=c(3,3,4,4,5), x2 = (0.1,0.2,0.3,0.4,0.4))
            ##     turns into constr = c(1,1,2,2,3,4,5,6,7,7) with seven distinct covariate effects
            ## Allow constraints such as: some elements are minus others. Use negative elements of constr to do this.
            ## e.g. constr = c(1,1,-1,-1,2,3,4,5)
            ## obtained by match(abs(x), unique(abs(x))) * sign(x)
            if (is.element(covlabels[i], names(constraint))) {
                if (length(constraint[[covlabels[i]]]) != nmatrix)
                    stop("\"",covlabels[i],"\" constraint of length ",
                         length(constraint[[covlabels[i]]]),", should be ",nmatrix)
            }
            else
                constraint[[covlabels[i]]] <- seq(nmatrix)
            constr <- c(constr, (maxc + match(abs(constraint[[covlabels[i]]]),
                                             unique(abs(constraint[[covlabels[i]]]))))*sign(constraint[[covlabels[i]]]) )
            maxc <- max(abs(constr))
        }
    }
    inits <- numeric()
    if (!is.null(covinits))
        msm.check.covinits(covinits, covdata)
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
    ndpars <- max(unique(abs(constr)))
    ## which covariate each distinct covariate parameter corresponds to. Used in C (FormDQCov)
    whichdcov <- rep(1:ncovs, each=nmatrix)[!duplicated(abs(constr))]
    list(npars=npars,
         ndpars=ndpars,    # number of distinct covariate effect parameters
         ncovs=ncovs,
         constr=constr,
         whichdcov=whichdcov,
         covlabels=covlabels, # factors as separate contrasts
         covlabels.orig=covlabels.orig, # factors as one variable
         inits = inits,
         covmeans = fulldata$covmeans[covlabels]
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
        ## Form death state info from hmmIdent parameters.
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

msm.form.obstype <- function(obstype, nobs, state, notna, dmodel, exacttimes)
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
            obstype[notna][state %in% dmodel$obs] <- 3
    }
    obstype[notna]
}

msm.form.obstrue <- function(obstrue, nobs, hmodel) {
    if (!is.null(obstrue)) {
        if (!is.numeric(obstrue) && !is.logical(obstrue)) stop("obstrue should be logical or numeric")
        else if (length(obstrue) != nobs) stop("obstrue of length ", length(obstrue), ", expected ", nobs)
        if (!hmodel$hidden) warning("Specified obstrue for a non-hidden model, ignoring.")
    }
    else if (hmodel$hidden) obstrue <- rep(0, nobs)
    else obstrue <- rep(1, nobs)
    obstrue
}

### Transform set of sets of probs {prs} to {log(prs/pr1)}
msm.mnlogit.transform <- function(pars, plabs, states){
    res <- pars
    if (any(plabs=="p")) {
        whichst <- match(states[plabs == "p"], unique(states[plabs == "p"]))
        res[plabs == "p"] <- log(pars[plabs=="p"] / pars[plabs=="pbase"][whichst])
    }
    res
}

### Transform set of sets of murs = {log(prs/pr1)} to probs {prs}
### ie psum = sum(exp(mus)),  pr1 = 1 / (1 + psum),  prs = exp(mus) / (1 + psum)
msm.mninvlogit.transform <- function(pars, plabs, states){
    res <- pars
    if (any(plabs=="p")) {
        whichst <- match(states[plabs=="p"], unique(states[plabs=="p"]))
        if (is.matrix(pars)) {# will be used when applying covariates
            for (i in unique(whichst)) {
###                psum <- colSums(exp(pars[plabs=="p",][whichst==i,,drop=FALSE]))
                psum <- colSums(exp(pars[which(plabs=="p")[whichst==i],,drop=FALSE]))
###                res[plabs=="pbase",,drop=FALSE][i,,drop=FALSE] <- 1 / (1 + psum)
                res[which(plabs=="pbase")[i],] <- 1 / (1 + psum)
###                res[plabs=="p",,drop=FALSE][whichst==i,,drop=FALSE] <-
                res[which(plabs=="p")[whichst==i],] <-
                    exp(pars[plabs=="p",,drop=FALSE][whichst==i,]) /
                        rep(1 + psum, each=sum(whichst==i))
            }
        }
        else {
            psum <- tapply(exp(pars[plabs=="p"]), states[plabs=="p"], sum)
            res[plabs=="pbase"] <- 1 / (1 + psum)
            res[plabs=="p"] <- exp(pars[plabs=="p"]) / (1 + psum[whichst])
        }
    }
    res
}

## transform parameters from natural scale to real-line optimisation scale

msm.transform <- function(pars, hmodel, ranges){
    labs <- names(pars)
    pars <- glogit(pars, ranges[,"lower"], ranges[,"upper"])
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")))
    hpars <- pars[hpinds]
    hpars <- msm.mnlogit.transform(hpars, hmodel$plabs, hmodel$parstate)
    pars[hpinds] <- hpars
    pars[labs=="initp"] <- log(pars[labs=="initp"] / pars[labs=="initpbase"])
    pars
}

## transform parameters from real-line optimisation scale to natural scale

msm.inv.transform <- function(pars, hmodel, ranges){
    labs <- names(pars)
    pars <- gexpit(pars, ranges[,"lower"], ranges[,"upper"])
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initp","initp0","initpcov")))
    hpars <- pars[hpinds]
    hpars <- msm.mninvlogit.transform(hpars, hmodel$plabs, hmodel$parstate)
    pars[hpinds] <- hpars
    ep <- exp(pars[labs=="initp"])
    pars[labs=="initp"] <- ep / (1 + sum(ep))
    pars[labs=="initpbase"] <- 1 / (1 + sum(ep))
    pars
}

## Collect all model parameters together ready for optimisation
## Handle parameters fixed at initial values or constrained to equal other parameters

msm.form.params <- function(qmodel, qcmodel, emodel, hmodel, fixedpars)
{
    ## Transition intensities
    ni <- qmodel$npars
    ## Covariates on transition intensities
    nc <- qcmodel$npars
    ## HMM response parameters
    nh <- sum(hmodel$npars)
    ## Covariates on HMM response distribution
    nhc <- sum(hmodel$ncoveffs)
    ## Initial state occupancy probabilities in HMM
    nip <- hmodel$nipars
    ## Covariates on initial state occupancy probabilities.
    nipc <- hmodel$nicoveffs
    npars <- ni + nc + nh + nhc + nip + nipc
    inits <- c(qmodel$inits, qcmodel$inits, hmodel$pars, unlist(hmodel$coveffect))
    plabs <- c(rep("qbase",ni), rep("qcov", nc), hmodel$plabs, rep("hcov", nhc))
    if (nip > 0) {
        inits <- c(inits, hmodel$initprobs)
        initplabs <- c("initpbase", rep("initp",nip-1))
        initplabs[hmodel$initprobs==0] <- "initp0" # those initialised to zero will be fixed at zero
        plabs <- c(plabs, initplabs)
        if (nipc > 0) {
            inits <- c(inits, unlist(hmodel$icoveffect))
            plabs <- c(plabs, rep("initpcov",nipc))
        }
    }
    ## store indicator for which parameters are HMM location parameters (not HMM cov effects or initial state probs)
    hmmpars <- which(!(plabs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")))
    hmmparscov <- which(!(plabs %in% c("qbase","qcov","initpbase","initp","initp0","initpcov")))
    names(inits) <- plabs
    ranges <- .msm.PARRANGES[plabs,,drop=FALSE]
    if (!is.null(hmodel$ranges)) ranges[hmmparscov,] <- hmodel$ranges
    inits <- msm.transform(inits, hmodel, ranges)
    ## Form constraint vector for complete set of parameters
    ## No constraints allowed on initprobs and their covs for the moment
    constr <- c(qmodel$constr, if(is.null(qcmodel$constr)) NULL else (ni + abs(qcmodel$constr))*sign(qcmodel$constr),
                ni + nc + hmodel$constr,
                ni + nc + nh + hmodel$covconstr, ni + nc + nh + nhc + seq(length=nip),
                ni + nc + nh + nip + seq(length=nipc))
    constr <- match(abs(constr), unique(abs(constr)))*sign(constr)
    ## parameters which are always fixed and not included in user-supplied fixedpars
    auxpars <- which(plabs %in% .msm.AUXPARS)
    duppars <- which(duplicated(abs(constr)))
    realpars <- setdiff(seq(npars), union(auxpars, duppars))
    nrealpars <- npars - length(auxpars) - length(duppars)
    ## if transition-specific covariates, then fixedpars indices generally smaller
    nshortpars <- nrealpars - sum(qcmodel$cri[!duplicated(qcmodel$constr)]==0)
    if (is.logical(fixedpars))
        fixedpars <- if (fixedpars == TRUE) seq(nshortpars) else numeric()
    if (any(! (fixedpars %in% seq(length=nshortpars))))
        stop ( "Elements of fixedpars should be in 1, ..., ", nshortpars)
    if (!is.null(qcmodel$cri)) {
        ## Convert user-supplied fixedpars indexing transition-specific covariates
        ## to fixedpars indexing transition-common covariates
        inds <- rep(1, nrealpars)
        inds[qmodel$ndpars + qcmodel$constr[!duplicated(qcmodel$constr)]] <-
            qcmodel$cri[!duplicated(qcmodel$constr)]
        inds[inds==1] <- seq(length=nshortpars)
        fixedpars <- match(fixedpars, inds)
        ## fix covariate effects not included in model to zero
        fixedpars <- sort(c(fixedpars, which(inds==0)))
    }
    fixedpars <- sort(c(realpars[fixedpars], auxpars))
    notfixed <- setdiff(seq(npars), fixedpars)
    allinits <- inits
    optpars <- intersect(notfixed, which(!duplicated(abs(constr))))
    inits <- inits[optpars]
    fixed <- (length(fixedpars) + length(duppars) == npars) # TRUE if all parameters are fixed, then no optimisation needed, just evaluate likelihood
    names(allinits) <- plabs; names(fixedpars) <- plabs[fixedpars]; names(plabs) <- NULL
    paramdata <- list(inits=inits, plabs=plabs, allinits=allinits, hmmpars=hmmpars,
                      fixed=fixed, notfixed=notfixed, optpars=optpars,
                      fixedpars=fixedpars, constr=constr,  npars=npars,
                      nfix=length(fixedpars),
                      nopt=length(optpars), ndup=length(duppars), ranges=ranges)
    paramdata
}

## Unfix and unconstrain all parameters in a paramdata object p (as
## returned by msm.form.params).  Used for calculating deriv /
## information over all parameters for models with fixed parameters.

msm.unfixallparams <- function(p) {
    p$fixed <- FALSE
    p$notfixed <- seq(along=p$allinits)
    p$optpars <- seq(along=p$allinits)
    p$fixedpars <- NULL
    p$constr <- seq(along=p$allinits)
    p$nfix <- 0
    p$nopt <- length(p$optpars)
    p$ndup <- 0
    p
}

msm.rep.constraints <- function(pars, # transformed pars
                                paramdata,
                                hmodel){
    plabs <- names(pars)
    p <- paramdata
    ## Before replication, transform probs on log(pr/pbase scale) back to pr scale,
    ## so that econstraint applies to pr not log(pr/pbase). After, transform back
    pars[p$hmmpars] <- msm.mninvlogit.transform(pars[p$hmmpars], hmodel$plabs, hmodel$parstate)
    plabs <- plabs[!duplicated(abs(p$constr))][abs(p$constr)]
    pars <- pars[!duplicated(abs(p$constr))][abs(p$constr)]*sign(p$constr)
    pars[p$hmmpars] <- msm.mnlogit.transform(pars[p$hmmpars], hmodel$plabs, hmodel$parstate)
    names(pars) <- plabs
    pars
}

## Apply covariates to HMM location parameters
## Parameters enter transformed, and exit on natural scale

msm.add.hmmcovs <- function(hmodel, pars, msmdata){
    labs <- names(pars)
    hpinds <- which(!(labs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")))
    hpars <- matrix(rep(pars[hpinds], msmdata$nobs), ncol=msmdata$nobs,
                    dimnames=list(labs[hpinds], NULL))
    ito <- 0
    for (i in which(hmodel$ncovs > 0)) {
        ifrom <- ito + 1; ito <- ito + hmodel$ncovs[i]
        hbase <- pars[hpinds][i]
        coveffs <- pars[names(pars)=="hcov"][ifrom:ito]
        wc <- hmodel$whichcovh[ifrom:ito]
        covs <- as.matrix(msmdata$covmat[,wc,drop=FALSE])
        hpars[i,] <- hpars[i,] + covs %*% coveffs
    }
    for (i in seq_along(hpinds)){
        hpars[i,] <- gexpit(hpars[i,], hmodel$ranges[i,"lower"], hmodel$ranges[i,"upper"]) # TESTME
#        if (labs[hpinds][i] %in% rownames(.msm.TRANSFORMS))  {
#            invlink <- get(.msm.TRANSFORMS[labs[hpinds][i],"inv"])
#            hpars[i,] <- invlink(hpars[i,])
#        }
    }
    hpars <- msm.mninvlogit.transform(hpars, hmodel$plabs, hmodel$parstate)
    hpars
}

## Apply covariates to transition intensities. Parameters enter this
## function already log transformed and replicated, and exit on
## natural scale.

msm.add.qcovs <- function(qmodel, pars, msmdata, agg=TRUE){
    labs <- names(pars)
    n <- if (agg) msmdata$nobs else msmdata$n
    qvec <- matrix(rep(pars[labs=="qbase"],each=n), nrow=n)
    covs <- if (agg) msmdata$covmat else msmdata$cov
    if (msmdata$covdata$ncovs > 0)
        covs <- as.matrix(covs[,msmdata$covdata$whichcov,drop=FALSE])
    beta <- matrix(pars[labs=="qcov"], ncol=qmodel$npars, byrow=TRUE)
    xb <- if (msmdata$covdata$ncovs == 0) matrix(0,nrow=n,ncol=qmodel$npars) else covs %*% beta
    qvec <- exp(qvec + xb)
    imat <- t(qmodel$imatrix); row <- col(imat)[imat==1]; col <- row(imat)[imat==1]
    qmat <- array(0, dim=c(qmodel$nstates, qmodel$nstates, n))
    for (i in 1:qmodel$npars) {
        qmat[row[i],col[i],] <- qvec[,i]
    }
    for (i in 1:qmodel$nstates)
        qmat[i,i,] <- -colSums(qmat[i,,])
    qmat
}

## Derivatives of intensity matrix Q wrt q (not log q) and beta. By
## observation with covariates applied.

msm.form.dq <- function(qmodel, qcmodel, pars, msmdata, agg=TRUE){
    labs <- names(pars)
    n <- if (agg) msmdata$nobs else msmdata$n
    qvec <- matrix(rep(pars[labs=="qbase"],each=n), nrow=n, ncol=qmodel$npars)
    covs <- if (agg) msmdata$covmat else msmdata$cov
    if (msmdata$covdata$ncovs > 0)
        covs <- as.matrix(covs[,msmdata$covdata$whichcov])
    beta <- matrix(pars[labs=="qcov"], ncol=qmodel$npars, byrow=TRUE)
    xb <- if (msmdata$covdata$ncovs == 0) matrix(0,nrow=n,ncol=qmodel$npars) else covs %*% beta
    qvec <- exp(qvec + xb)
    qrvec <- exp(xb)
    dqvec <- array(pars[labs=="qbase"], dim=c(n, qmodel$npars, qmodel$npars + qcmodel$npars))
    for (i in seq_len(qmodel$npars)) { ## can this be vectorised further?
        dqvec[,i,] = qrvec[,i] * rep(abs(qmodel$constr) == abs(qmodel$constr[i]), each=n)
        if (qcmodel$npars > 0)
            for (k in 1:qcmodel$ncovs) {
                con <- qcmodel$constr[(k-1)*qmodel$npars + 1:qmodel$npars]
                con <- match(abs(con), unique(abs(con)))*sign(con)
                dqvec[,i,k*qmodel$npars + 1:qmodel$npars] <-
                    covs[,k] * qvec[,i] * rep((abs(con)==abs(con[i]))*sign(con)*sign(con[i]), each=n)
            }
    }
    dqmat <- array(0, dim=c(qmodel$nstates, qmodel$nstates, qmodel$npars + qcmodel$npars, n))
    imat <- t(qmodel$imatrix); row <- col(imat)[imat==1]; col <- row(imat)[imat==1]
    for (i in 1:qmodel$npars)
        dqmat[row[i],col[i],,] <- t(dqvec[,i,])
    for (i in 1:qmodel$nstates)
        dqmat[i,i,,] <- -colSums(dqmat[i,,,], dims=1)
    dqmat
}


msm.initprobs2mat <- function(hmodel, pars, msmdata){
    ## Convert vector initial state occupancy probs to matrix by patient
    if (!hmodel$hidden) return(0)
    if (hmodel$est.initprobs) {
        initp <- pars[names(pars) %in% c("initpbase","initp","initp0")]
        initp <- matrix(rep(initp, each=msmdata$npts), nrow=msmdata$npts,
                        dimnames=list(NULL,
                        names(pars)[names(pars) %in% c("initpbase","initp","initp0")]))
        ## Multiply baselines (entering on mnlogit scale) by current covariate
        ## effects, giving matrix of patient-specific initprobs
        est <- which(colnames(initp)=="initp")
        ip <- initp[,est,drop=FALSE]
        if (hmodel$nicoveffs > 0) {
            ## cov effs ordered by states (excluding state 1) within covariates
            coveffs <- pars[names(pars)=="initpcov"]
            coveffs <- matrix(coveffs, nrow=max(hmodel$nicovs), byrow=TRUE)
            coveffs[,hmodel$nicovs==0] <- 0 # states with fixed initp=zero
            firstobs <- !duplicated(msmdata$subject)
            ip <- ip + as.matrix(msmdata$covmat[firstobs,hmodel$whichcovi,drop=FALSE]) %*% coveffs
        }
        initp[,est] <- exp(ip) / (1 + rowSums(exp(ip)))
        initp[,"initpbase"] <- 1 / (1 + rowSums(exp(ip)))
    }
    else if (!is.matrix(hmodel$initprobs))
        initp <- matrix(rep(hmodel$initprobs,each=msmdata$npts),nrow=msmdata$npts)
    else initp <- hmodel$initprobs
    initp
}

## Entry point to C code for calculating the likelihood and related quantities

Ccall.msm <- function(params, do.what="lik", msmdata, qmodel, qcmodel, cmodel, hmodel, paramdata)
{
    p <- paramdata
    pars <- p$allinits
    pars[p$optpars] <- params
    pars <- msm.rep.constraints(pars, paramdata, hmodel)

    ## Add covariates to hpars and q here. Inverse-transformed to natural scale on exit
    agg <- if(do.what %in% c("lik","deriv","info")) TRUE else FALSE # data as aggregate transition counts, as opposed to individual observations
    Q <- msm.add.qcovs(qmodel, pars, msmdata, agg)
    DQ <- if (do.what %in% c("deriv","info","deriv.subj","dpmat")) msm.form.dq(qmodel, qcmodel, pars, msmdata, agg) else NULL
    hpars <- if (hmodel$hidden) msm.add.hmmcovs(hmodel, pars, msmdata) else NULL
    initprobs <- msm.initprobs2mat(hmodel, pars, msmdata)

    ## In R, work with states / parameter indices / model indices 1, ... n. In C, work with 0, ... n-1
    msmdata$fromstate <- msmdata$fromstate - 1
    msmdata$tostate <- msmdata$tostate - 1
    msmdata$firstobs <- msmdata$firstobs - 1
    hmodel$models <- hmodel$models - 1
    hmodel$links <- hmodel$links - 1

    lik <- .C("msmCEntry",
              as.integer(match(do.what, .msm.CTASKS) - 1),
              as.double(Q),
              as.double(DQ),
              as.double(hpars),

              ## data for non-HMM
              as.integer(msmdata$fromstate),
              as.integer(msmdata$tostate),
              as.double(msmdata$timelag),
              as.integer(msmdata$nocc),
              as.integer(msmdata$noccsum),
              as.integer(msmdata$whicha),
              as.integer(if(do.what %in% c("lik.subj","deriv.subj","dpmat")) msmdata$obstype.obs else msmdata$obstype), # per-observation obstype

              ## data for HMM or censored
              as.integer(match(msmdata$subject, unique(msmdata$subject))),
              as.double(msmdata$time),
              as.double(msmdata$state), # If this is a misc or censored state, this is indexed from 1.
              as.integer(msmdata$firstobs),
              as.integer(msmdata$obstrue),

              ## HMM specification
              as.integer(hmodel$hidden),
              as.integer(hmodel$models),
              as.integer(hmodel$totpars),
              as.integer(hmodel$firstpar),
              as.double(initprobs),

              ## various constants
              as.integer(qmodel$nstates),
              as.integer(qmodel$iso),
              as.integer(qmodel$perm),
              as.integer(qmodel$qperm),
              as.integer(qmodel$expm),
              as.integer(qmodel$npars),
              as.integer(msmdata$nobs), # number of aggregated transitions
              as.integer(msmdata$n), # number of observations
              as.integer(msmdata$npts),  # HMM only
              as.integer(msmdata$ntrans), # number of (disaggregated) transitions
              as.integer(qcmodel$npars),

              as.integer(cmodel$ncens),
              as.integer(cmodel$censor),
              as.integer(cmodel$states),
              as.integer(cmodel$index - 1),

              returned = double(
              if (do.what %in% c("deriv","info")) qmodel$npars + qcmodel$npars
              else if (do.what=="lik.subj") msmdata$npts
              else if (do.what=="deriv.subj") msmdata$npts*(qmodel$npars + qcmodel$npars)
              else if (do.what=="dpmat") msmdata$ntrans*qmodel$nstates*(qmodel$npars + qcmodel$npars)
              else if (do.what=="viterbi") msmdata$nobs
              else 1),

              returned2 = double(
              if (do.what=="info") (qmodel$npars + qcmodel$npars) ^ 2
              else 1),

              ## so that Inf values are allowed for parameters denoting truncation points of truncated distributions
              NAOK = TRUE
#              ,
#              PACKAGE = "msm"
              )
    ## transform derivatives wrt Q to derivatives wrt log Q
    ## don't return derivs for constrained or fixed parameters
    qp <- 1:qmodel$npars
    dp <- intersect(which(!duplicated(abs(p$constr))), p$notfixed)
    tpars <- if(length(params)==0) p$allinits[qp] else (params[abs(p$constr)]*sign(p$constr))[qp]
    if (do.what %in% c("deriv","info")) {
        lik$returned[qp] <- lik$returned[qp]*exp(tpars)
        lik$returned <- lik$returned[dp]
    }
    ## Fisher information matrix
    if (do.what=="info") {
        lik$returned2 <- matrix(lik$returned2, nrow=qmodel$npars + qcmodel$npars)
        lik$returned2[qp,qp] <- lik$returned2[qp,qp]*outer(exp(tpars),exp(tpars))
        if (qcmodel$ndpars > 0) {
            qcp <- (qmodel$npars + 1):(qmodel$npars + qcmodel$npars)
            lik$returned2[qp,qcp] <- lik$returned2[qp,qcp]*exp(tpars)
            lik$returned2[qcp,qp] <- lik$returned2[qcp,qp]*rep(exp(tpars),each=qcmodel$npars)
            lik$returned2 <- lik$returned2[dp,dp]
        }
    }
    ## subject-specific derivatives, to use for score residuals
    if (do.what=="deriv.subj") {
        lik$returned <- matrix(lik$returned, nrow=msmdata$npts)
        lik$returned[,qp] <- lik$returned[,qp]*rep(exp(tpars), each=msmdata$npts)
        lik$returned <- lik$returned[,dp]
    }
    ## transition-specific derivatives of P matrix, to use for Pearson test
    if (do.what=="dpmat") {
        lik$returned <- array(lik$returned, dim=c(msmdata$ntrans, qmodel$nstates, qmodel$npars+qcmodel$npars))
        tpars <- if(length(params)==0) p$allinits[!duplicated(p$constr)][qp] else params[qp]
        lik$returned[,,qp] <- lik$returned[,,qp]*rep(exp(tpars), each=msmdata$ntrans*qmodel$nstates)
        lik$returned <- lik$returned[,,dp]
    }
    if (do.what=="info") list(deriv=lik$returned, info=lik$returned2) else lik$returned
}

lik.msm <- function(params, ...)
{
    Ccall.msm(params, do.what="lik", ...)
}

deriv.msm <- function(params, ...)
{
    Ccall.msm(params, do.what="deriv", ...)
}

information.msm <- function(params, ...)
{
    Ccall.msm(params, do.what="info", ...)
}

## Convert vector of MLEs into matrices

msm.form.output <- function(whichp, model, cmodel, p)
{
    Matrices <- MatricesSE <- MatricesL <- MatricesU <- list()
    for (i in 0:cmodel$ncovs) {
        matrixname <- if (i==0) "logbaseline" else cmodel$covlabels[i] # name of the current output matrix.
        mat <- t(model$imatrix) # state matrices filled by row, while R fills them by column.
        if (whichp=="intens")
            parinds <- if (i==0) which(p$plabs=="qbase") else which(p$plabs=="qcov")[(i-1)*model$npars + 1:model$npars]
        if (whichp=="misc")
            parinds <- if (i==0) which(p$plabs=="p") else which(p$plabs=="hcov")[i + cmodel$ncovs*(1:model$npars - 1)]
        if (any(parinds)) mat[t(model$imatrix)==1] <- p$params[parinds]
        else mat[mat==1] <- Inf ## if no parinds are "p", then there are off-diag 1s in ematrix
        mat <- t(mat)
        dimnames(mat) <- dimnames(model$imatrix)
        if (p$foundse && !p$fixed){
            intenscov <- p$covmat[parinds, parinds]
            intensse <- sqrt(diag(as.matrix(intenscov)))
            semat <- lmat <- umat <- t(model$imatrix)
            if (any(parinds)){
                semat[t(model$imatrix)==1] <- intensse
                lmat[t(model$imatrix)==1] <- p$ci[parinds,1]
                umat[t(model$imatrix)==1] <- p$ci[parinds,2]
            }
            else semat[semat==1] <- lmat[lmat==1] <- umat[umat==1] <- Inf
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

## Format hidden Markov model estimates and CIs

msm.form.houtput <- function(hmodel, p, msmdata)
{
    hmodel$pars <- p$estimates.t[!(p$plabs %in% c("qbase","qcov","hcov","initp","initpbase","initp0","initpcov"))]
    hmodel$coveffect <- p$estimates.t[p$plabs == "hcov"]
    hmodel$fitted <- !p$fixed
    hmodel$foundse <- p$foundse
    if (hmodel$nip > 0) {
        if (hmodel$foundse) {
            hmodel$initprobs <- rbind(cbind(p$estimates.t[p$plabs %in% c("initpbase","initp","initp0")],
                                            p$ci[p$plabs %in% c("initpbase","initp","initp0"),,drop=FALSE]))
            rownames(hmodel$initprobs) <- paste("State",1:hmodel$nstates)
            colnames(hmodel$initprobs) <- c("Estimate", "LCL", "UCL")
            if (any(hmodel$nicovs > 0)) {
                covnames <- names(hmodel$icoveffect)
                hmodel$icoveffect <- cbind(p$estimates.t[p$plabs == "initpcov"],  p$ci[p$plabs == "initpcov",,drop=FALSE])
                iplabs <- p$plabs[p$plabs %in% c("initp","initp0")]
                whichst <- which(iplabs == "initp") + 1  # init probs for which states have covs on them (not the zero probs)
                rownames(hmodel$icoveffect) <- paste(covnames, paste("State",whichst), sep=", ")
                colnames(hmodel$icoveffect) <- c("Estimate", "LCL", "UCL")
            }
        }
        else {
            hmodel$initprobs <- c(1 - sum(p$estimates.t[p$plabs == "initp"]),
                                  p$estimates.t[p$plabs %in% c("initp","initp0")])
            names(hmodel$initprobs) <- paste("State", 1:hmodel$nstates)
            if (any(hmodel$nicovs > 0)) {
                covnames <- names(hmodel$icoveffect)
                hmodel$icoveffect <- p$estimates.t[p$plabs == "initpcov"]
                names(hmodel$icoveffect) <- paste(covnames, paste("State",2:hmodel$nstates), sep=", ")
            }
        }
    }
    hmodel$initpmat <- msm.initprobs2mat(hmodel, p$estimates.t, msmdata)
    if (hmodel$foundse) {
        hmodel$ci <- p$ci[!(p$plabs %in% c("qbase","qcov","hcov","initpbase","initp","initp0","initpcov")), , drop=FALSE]
        hmodel$covci <- p$ci[p$plabs %in% c("hcov"), ]
    }
    names(hmodel$pars) <- hmodel$plabs
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
    mf <- model.frame(formula, data=data, na.action=NULL)
    state <- mf[,1]
    time <- mf[,2]
    n <- length(state)
    if (missing(subject)) subject <- rep(1, n)
    if (!is.null(data))
        subject <- eval(substitute(subject), as.list(data), parent.frame())
    subject <- match(subject, unique(subject))
    notna <- !is.na(subject) & !is.na(time) & !is.na(state)
    subject <- subject[notna]; time <- time[notna]; state <- state[notna]
    msm.check.qmatrix(qmatrix)
    msm.check.state(nrow(qmatrix), state, cens$censor)
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
    rownames(estmat) <- rownames(qmatrix)
    colnames(estmat) <- colnames(qmatrix)
    estmat
}

### Construct a model with time-dependent transition intensities.
### Form a new dataset with censored states and extra covariate, and
### form a new censor model, given change times in tcut

msm.pci <- function(tcut, dat, qmodel, cmodel, center)
{
    if (!is.numeric(tcut)) stop("Expected \"tcut\" to be a numeric vector of change points")
    old <- as.data.frame(dat[c("subject","time","state","obstype","obstrue")])
    if (dat$ncovs > 0)
        old[c("covmat","covmat.orig")] <- dat[c("covmat","covmat.orig")]
    ## new dataset
    ntcut <- length(tcut)
    nextra <- ntcut*dat$npts
    extra <- data.frame(subject = rep(unique(dat$subject), each=ntcut),
                        time = rep(tcut, dat$npts),
                        state = rep(NA, nextra),
                        obstype = rep(1, nextra),
                        obstrue = rep(TRUE, nextra),
                        pci.imp = 1
                        )
    old$pci.imp <- 0
    if (dat$ncovs > 0){
        extra$covmat <- as.data.frame(matrix(NA, nrow=nextra, ncol=ncol(old$covmat)))
        extra$covmat.orig <- as.data.frame(matrix(NA, nrow=nextra, ncol=ncol(old$covmat.orig)))
        rownames(old$covmat.orig) <- rownames(old$covmat) <- 1:nrow(old)
        rownames(extra$covmat.orig) <- rownames(extra$covmat) <- nrow(old) + 1:nrow(extra) # get rid of dup / null rownames errors
    }
    ## merge new and old observations
    new <- rbind(old, extra)
    new <- new[order(new$subject, new$time),]
    label <- if (cmodel$ncens > 0) max(cmodel$censor)*2 else qmodel$nstates + 1
    new$state[is.na(new$state)] <- label
    ## Only keep cutpoints within range of each patient's followup
    mintime <- tapply(old$time, old$subject, min)[as.character(unique(old$subject))]
    maxtime <- tapply(old$time, old$subject, max)[as.character(unique(old$subject))]
    nobspt <- as.numeric(table(new$subject)[as.character(unique(new$subject))])
    new <- new[new$time >= rep(mintime, nobspt) & new$time <= rep(maxtime, nobspt), ]

    ## drop imputed observations at times when there was already an observation
    ## assumes there wasn't already duplicated obs times
    prevsubj <- c(NA,new$subject[1:(nrow(new)-1)]); nextsubj <- c(new$subject[2:nrow(new)], NA)
    prevtime <- c(NA,new$time[1:(nrow(new)-1)]); nexttime <- c(new$time[2:nrow(new)], NA)
    prevstate <- c(NA,new$state[1:(nrow(new)-1)]); nextstate <- c(new$state[2:nrow(new)], NA)
    new <- new[!((new$subject==prevsubj & new$time==prevtime & new$state==label & prevstate!=label) |
                 (new$subject==nextsubj & new$time==nexttime & new$state==label & nextstate!=label))
               ,]

    ## Carry last value forward for other covariates
    if (dat$ncovs > 0) {
        eind <- which(is.na(new$covmat[,1]))
        while(length(eind) > 0){
            new$covmat[eind,] <- new$covmat[eind - 1,]
            new$covmat.orig[eind,] <- new$covmat.orig[eind - 1,]
            eind <- which(is.na(new$covmat[,1]))
        }
    }

    ## constants in dataset
    new <- as.list(new)
    new$nobs <- new$n <- n <- length(new$state)
    new$npts <- dat$npts
    new$firstobs <- c(1, which(new$subject[2:n] != new$subject[1:(n-1)]) + 1, n+1)
    new$ncovs <- dat$ncovs + ntcut

    ## Check range of cut points
    if (any(tcut <= min(dat$time)))
      warning("Time cut point", if (sum(tcut <= min(dat$time)) > 1) "s " else " ",
                paste(tcut[tcut<=min(dat$time)],collapse=","),
                " less than or equal to minimum observed time of ",min(dat$time))
    if (any(tcut >= max(dat$time)))
        warning("Time cut point", if (sum(tcut >= max(dat$time)) > 1) "s " else " ",
                paste(tcut[tcut>=max(dat$time)],collapse=","),
                " greater than or equal to maximum observed time of ",max(dat$time))
    tcut <- tcut[tcut > min(dat$time) & tcut < max(dat$time)]
    ntcut <- length(tcut)
    if (ntcut==0)
        res <- NULL # no cut points in range of data, continue with no time-dependent model
    else {
        ## Insert new covariate in data representing time period
        tcovlabel <- "timeperiod"
        while (tcovlabel %in% dat$covlabels)
            tcovlabel <- paste(tcovlabel, ".1", sep="")
        tcov <- factor(cut(new$time, c(-Inf,tcut,Inf), right=FALSE))
        levs <- levels(tcov)
        levels(tcov) <- gsub(" ","", levs) # get rid of spaces in e.g. [10, Inf) levels
        assign(tcovlabel, tcov)
        mm <- model.matrix(as.formula(paste("~", tcovlabel)))[,-1,drop=FALSE]
        new$covmat <- cbind(new$covmat, mm)
        lastobs <- c(new$subject[1:(new$nobs-1)] != new$subject[2:new$nobs], TRUE)
        new$covmeans <- colMeans(new$covmat[!lastobs,,drop=FALSE])
        new$covmat.orig <- if(is.null(new$covmat.orig)) data.frame(timeperiod=tcov) else cbind(new$covmat.orig, timeperiod=tcov)
        if (center) new$covmat <- sweep(new$covmat, 2, colMeans(new$covmat))

        ## new censoring model
        cmodel$ncens <- cmodel$ncens + 1
        cmodel$censor <- c(cmodel$censor, label)
        cmodel$states <- c(cmodel$states, 1:qmodel$nstates)
        cmodel$index <- if (is.null(cmodel$index)) 1 else cmodel$index
        cmodel$index <- c(cmodel$index, length(cmodel$states) + 1)

        ## new auxiliary information about covariates
        for (i in c("covdata","misccovdata","hcovdata","icovdata")) new[[i]] <- dat[[i]]
        new$covlabels <- c(dat$covlabels, colnames(mm))
        new$covlabels.orig <- c(dat$covlabels.orig, tcovlabel)
        new$covdata$covlabels <- c(dat$covdata$covlabels, colnames(mm))
        new$covdata$ncovs <- dat$covdata$ncovs + ntcut
        new$covdata$covmeans <- c(dat$covdata$covmeans, colMeans(mm))
        new$covdata$covfactor <- c(dat$covdata$covfactor, timeperiod=TRUE)
        new$covdata$covfactorlevels <- c(dat$covdata$covfactorlevels, list(timeperiod=levels(tcov)))
        new$covdata$covlabels.orig <- c(dat$covdata$covlabels.orig, tcovlabel)
        new$covdata$whichcov <- match(new$covdata$covlabels, new$covlabels)
        new$covdata$whichcov.orig <- match(new$covdata$covlabels.orig, new$covlabels.orig)

        res <- list(dat=new, cmodel=cmodel, tcut=tcut)
    }
    res
}


msm.check.covlist <- function(covlist) {
    check.numnum <- function(str)
        length(grep("^[0-9]+-[0-9]+$", str)) == length(str)
    num <- sapply(names(covlist), check.numnum)
    if (!all(num)) {
        badnums <- which(!num)
        plural1 <- if (length(badnums)>1) "s" else "";
        plural2 <- if (length(badnums)>1) "e" else "";
        badnames <- paste(paste("\"",names(covlist)[badnums],"\"",sep=""), collapse=",")
        badnums <- paste(badnums, collapse=",")
        stop("Name", plural1, " ", badnames, " of \"covariates\" formula", plural2, " ", badnums, " not in format \"number-number\"")
    }
    for (i in seq(along=covlist))
        if (!inherits(covlist[[i]], "formula"))
            stop("\"covariates\" should be a formula or list of formulae")
}

## Given covariates argument to msm as a list of transition-specific
## formulae, convert to formula common to all intensities, with
## associated data.  This can be processed with existing
## msm.form.covdata.  Appropriate effects will be fixed to zero in
## msm.form.params.

msm.form.covdata.byrate <- function(covlist, qemodel, data, ignore.obs=NULL, center=TRUE) {
    msm.check.covlist(covlist)
    trans <- sapply(strsplit(names(covlist), "-"), as.numeric)
    npars <- qemodel$npars
    tm <- if(inherits(qemodel,"msmqmodel")) "transition" else "misclassification"
    qe <- if(inherits(qemodel,"msmqmodel")) "qmatrix" else "ematrix"
    imat <- qemodel$imatrix
    for (i in seq(length=ncol(trans))){
        if (imat[trans[1,i],trans[2,i]] != 1)
            stop("covariates on ", names(covlist)[i], " ", tm, " requested, but this is not permitted by the ", qe, ".")
    }
    imat <- t(imat) # order named transitions / misclassifications by row
    tnames <- paste(col(imat)[imat==1],row(imat)[imat==1],sep="-")

    ## Merge transition-specific formulae into one big formula
    ter <- lapply(covlist, terms)
    vars <- unique(unlist(lapply(ter, function(x)attr(x,"term.labels"))))
    form <- as.formula(paste("~ ", paste(vars,collapse="+")))
    covdata <- msm.form.covdata(form, data, ignore.obs, center)

    ## Form indicator matrix for cov effects that will be fixed to zero
    covs <- covdata$covlabels
    covdata$cri <- matrix(0, nrow=npars, ncol=length(covs), dimnames = list(tnames, covs))
    sorti <- function(x) {
        ## converts, e.g. c("b:a:c","d:f","f:e") to c("a:b:c", "d:f", "e:f")
        sapply(lapply(strsplit(x, ":"), sort), paste, collapse=":")
    }
    for (i in 1:npars) {
        if (tnames[i] %in% names(covlist)) {
            covsi <- msm.form.covdata(covlist[[tnames[i]]], data)$covlabels
            covdata$cri[i, match(sorti(covsi), sorti(covs))] <- 1
        }
    }
    covdata
}

## Process constraints and initial values for covariates supplied as a
## list of transition-specific formulae.  Convert to form needed for a
## single covariates formula common to all transitions, which can be
## processed with msm.form.covmodel.

msm.form.covmodel.byrate <- function(covdata,
                                     constraint, # as supplied by user
                                     nmatrix,
                                     covinits,  # as supplied by user
                                     fulldata
                                     ){
    covs <- covdata$covlabels
    cri <- covdata$cri
    ## Convert short form constraints to long form
    msm.check.constraint(constraint, covdata)
    constr <- inits <- numeric()
    for (i in seq(along=covs)){
        if (covs[i] %in% names(constraint)){
            if (length(constraint[[covs[i]]]) != sum(cri[,i]))
                stop("\"",covs[i],"\" constraint of length ",
                     length(constraint[[covs[i]]]),", should be ",sum(cri[,i]))
            con <- match(constraint[[covs[i]]], unique(constraint[[covs[i]]])) + 1
            constraint[[covs[i]]] <- rep(1, nmatrix)
            constraint[[covs[i]]][cri[,i]==1] <- con
        }
        else constraint[[covs[i]]] <- seq(length=nmatrix)
    }
    ## convert short to long initial values in the same way
    if (!is.null(covinits)) msm.check.covinits(covinits, covdata)
    for (i in seq(along=covs)) {
        if (!is.null(covinits) && (covs[i] %in% names(covinits))) {
            if (!is.numeric(covinits[[covs[i]]])) {
                warning("initial values for covariates should be numeric, ignoring")
                covinits[[covs[i]]] <- rep(0, nmatrix)
            }
            thisinit <- rep(0, nmatrix)
            if (length(covinits[[covs[i]]]) != sum(cri[,i])) {
                warning("\"", covs[i], "\" initial values of length ", length(covinits[[covs[i]]]), ", should be ", sum(cri[,i]), ", ignoring")
                covinits[[covs[i]]] <- rep(0, nmatrix)
            }
            else thisinit[cri[,i]==1] <- covinits[[covs[i]]]
            covinits[[covs[i]]] <- thisinit
        }
    }
    covdata$cri <- NULL
    qcmodel <- msm.form.covmodel(covdata, constraint, nmatrix, covinits, fulldata)
    qcmodel$cri <- cri
    qcmodel
}

### Unload shared library when package is detached with unloadNamespace("msm")
.onUnload <- function(libpath) { library.dynam.unload("msm", libpath) }

