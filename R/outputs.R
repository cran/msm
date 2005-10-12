### METHODS FOR MSM OBJECTS
### NEW VERSION WITH CIs 

print.msm <- function(x, ...)
  {
      cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
      
      if (!attr(x,"fixed")) { 
          cat ("Maximum likelihood estimates: \n")
          covmessage <- if (x$qcmodel$ncovs == 0) "" else "with covariates set to their means"
          for (i in c("baseline", x$qcmodel$covlabels)) {
              title <-
                if (i == "baseline") paste("Transition intensity matrix",covmessage,"\n")
                else paste("Log-linear effects of", i, "\n")
              cat (title, "\n")
              print.ci(x$Qmatrices[[i]], x$QmatricesL[[i]], x$QmatricesU[[i]])
              cat("\n")
          }
          if (x$emodel$misc) {
              for (i in c("baseline", x$ecmodel$covlabels)) {
                  title <-
                    if (i == "baseline") paste("Misclassification matrix",covmessage,"\n")
                    else paste("Logit-linear effects of", i, "\n")
                  cat (title, "\n")
                  print.ci(x$Ematrices[[i]], x$EmatricesL[[i]], x$EmatricesU[[i]])
                  cat("\n")
              }
          }
          else if (x$hmodel$hidden) {print(x$hmodel); cat("\n")}
      }
      cat ("-2 * log-likelihood: ", x$minus2loglik, "\n")

  }

summary.msm <- function(object, # fitted model
                        times = NULL,  # times at which to compare observed and expected prevalences
                        timezero = NULL,
                        initstates = NULL,
                        covariates = "mean", 
                        misccovariates = "mean", 
                        hazard.scale = 1,
                        ...
                        )
  {
      if (!inherits(object, "msm")) stop("expected object to be a msm model")
      prevalences <- prevalence.msm(object, times, timezero, initstates, covariates, misccovariates)
      if (object$qcmodel$ncovs > 0) {
          if (missing (hazard.scale))
            hazard.scale <- rep(1, object$data$covdata$ncovs)
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
          cat("\nObserved numbers of individuals occupying states at each time\n\n")
          print(x$prevalences$Observed)
          cat("\nExpected numbers of individuals occupying states at each time\n\n")
          print(x$prevalences$Expected)
          cat("\nObserved prevalences of states (percentages of population at risk)\n\n")
          print(x$prevalences$"Observed percentages")
          cat("\nExpected prevalences of states (percentages of population at risk)\n\n")
          print(x$prevalences$"Expected percentages")
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

plot.msm <- function(x, from=NULL, to=NULL, range=NULL, covariates="mean", legend.pos=NULL, ...)
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
      if (is.null(from))
        from <- transient.msm(x)
      else {
          if (!is.numeric(from)) stop("from must be numeric")
          if (any (! (from %in% 1:x$qmodel$nstates ) ) )
            stop("from must be a vector of states in 1, ..., ", x$qmodel$nstates)
      }
      if (is.null(to))
        to <- max(absorbing.msm(x))
      else {
          if (!is.numeric(to)) stop("to must be numeric")
          if (! (to %in% absorbing.msm(x) ) ) stop("to must be an absorbing state")
      }
      if (is.null(range)) {
          rg <- range(x$data$time)
      }
      else {
          if (!is.numeric(range) || length(range)!= 2) stop("range must be a numeric vector of two elements")
          rg <- range
      }
      timediff <- (rg[2] - rg[1]) / 50
      times <- seq(rg[1], rg[2], timediff)
      pr <- numeric()
      cols <- rainbow(length(from))
      for (t in times)
        pr <- c(pr, pmatrix.msm(x, t, covariates)[from[1], to])
      plot(times, 1 - pr, type="l", xlab="Time", ylab="Fitted survival probability",
           ylim=c(0,1), lty = 1, col=cols[1],...)
      lt <- 2
      for (st in from[-1]){ 
          pr <- numeric()
          for (t in times)
            pr <- c(pr, pmatrix.msm(x, t, covariates)[st, to])
          lines(times, 1 - pr, type="l", lty = lt, col=cols[lt],...)
          lt <- lt+1
      }
      if (!is.numeric(legend.pos) || length(legend.pos) != 2)
        legend.pos <- c(max(times) - 15*timediff, 1)
      legend(legend.pos[1], legend.pos[2], legend=paste("From state",from), lty = seq(lt-1), col=cols, lwd=2)
      invisible()
  }


### Extract the transition intensity matrix at given covariate values

qmatrix.msm <- function(x, # fitted msm model
                        covariates = "mean",  # covariate values to calculate transition matrix for
                        sojourn = FALSE,      # also calculate mean sojourn times and their standard errors
                        cl=0.95
                        )
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
      qematrix.msm(x, covariates, intmisc="intens", sojourn=sojourn, cl=cl)
  }

### Extract the misclassification probability matrix at given covariate values

ematrix.msm <- function(x,
                        covariates = "mean",
                        cl=0.95)
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
      if (!x$emodel$misc) NULL
      else
        qematrix.msm(x, covariates, intmisc="misc", sojourn=FALSE, cl=cl)
  }

### Extract either intensity or misclassification matrix

qematrix.msm <- function(x, covariates="mean", intmisc="intens", sojourn=FALSE, cl=0.95)
{
    nst <- x$qmodel$nstates
    if (intmisc=="intens"){
        ni <- x$qmodel$npars
        nc <- if (!is.list(covariates) && (covariates == "mean")) 0 else x$qcmodel$ncovs
        if (nc==0){
            logest <- x$Qmatrices$logbaseline
            mat <- exp(logest)
            mat[x$qmodel$imatrix == 0] <- 0
            mat <- msm.fixdiag.qmatrix(mat)
        }
        covlabels <- x$qcmodel$covlabels
        covmeans <- x$qcmodel$covmeans
    }
    else if (intmisc=="misc"){
        ni <- x$emodel$npars
        nc <- if (!is.list(covariates) && (covariates == "mean")) 0 else x$ecmodel$ncovs
        if (nc==0){
            logest <- x$Ematrices$logitbaseline
            mat <- plogis(logest)
            mat[x$emodel$imatrix == 0] <- 0
            mat <- msm.fixdiag.ematrix(mat)
        }
        covlabels <- x$ecmodel$covlabels
        covmeans <- x$ecmodel$covmeans
    }

    se <- lse <- numeric(ni)
    if ((nc == 0) && (is.list(covariates)) && (length(covariates) > 0))
      warning(paste("Ignoring covariates - no covariates in model for",
                    if (intmisc=="intens") "transition intensities" else "misclassification probabilities"))
    if (nc > 0) {
        if (!is.list(covariates) && covariates == 0)
          {covariates <- list();  for (i in 1 : nc) covariates[[covlabels[i]]] <- 0}
        if (!is.list(covariates))
          stop("covariates argument must be 0, \"mean\", or a list of values for each named covariate")
        if (length(covariates) != nc)
          stop(paste("Supplied", length(covariates), "covariate values, expected", nc))
        for (n in names(covariates))
          if ((n != "") && (! n %in% covlabels) )
            stop(paste("Covariate \"", n, "\" not in model. For factor covariates, use, for example, covnameCOVVALUE = 1", sep=""))
        if (intmisc=="intens"){
            logest <- x$Qmatrices$logbaseline
            for (i in 1 : nc)
              logest <- logest + x$Qmatrices[[i+1]] * (covariates[[i]] - covmeans[i])
            mat <- exp(logest)
            mat[x$qmodel$imatrix == 0] <- 0
            mat <- msm.fixdiag.qmatrix(mat)
        }
        else if (intmisc=="misc"){
            logest <- x$Ematrices$logitbaseline
            for (i in 1 : nc)
              logest <- logest + x$Ematrices[[i+1]] * (covariates[[i]] - covmeans[i])
            mat <- plogis(logest)
            mat[x$emodel$imatrix == 0] <- 0
            mat <- msm.fixdiag.ematrix(mat)
        }        
    }
    if (sojourn) soj = -1 / diag(mat)
    
    if (x$foundse) {
        ## Work out standard errors. 
        ## Transformation for delta method is
        ##  exp (or expit) (x1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
        ## Use delta method to find approximate SE of the transform on log scale
        ## Work out a CI by for this assuming normal and transforming back
        coefs <- if (!is.list(covariates) && (covariates=="mean")) 1 else c(1, unlist(covariates) - covmeans)
        trsum <- if (intmisc=="intens") expsum else expitsum
        form <- as.formula(paste("~", trsum(seq(nc + 1), coefs)))
        lform <- as.formula(paste("~", lsum(seq(nc + 1), coefs)))
        ## indices into estimates vector of all intens/miscs, intens covs / misc covs
        inds <-
          if (intmisc=="intens") seq(length=x$qmodel$npars + x$qcmodel$npars)
          else x$qmodel$npars + x$qcmodel$npars + c(x$hmodel$constr[x$hmodel$plabs == "p"],
                                                    x$hmodel$totpars + seq(along=x$hmodel$covconstr))
        for (i in 1 : ni){
            ## indices into estimates vector of all intens/miscs, intens covs / misc covs for that particular fromstate-tostate. 
            parinds <-
              if (intmisc=="intens") inds[seq(i, (nc * ni + i), ni)]
              else inds[c(i, seq(ni + (i-1)*nc + 1, length=nc) )]
            ests <- x$estimates[parinds]
            cov <- x$covmat[parinds, parinds]
            se[i] <- deltamethod(form, ests, cov)
            lse[i] <- deltamethod(lform, ests, cov)
        }
        
        semat <- lsemat <- lmat <- umat <- matrix(0, nst, nst)
        ivector <- as.numeric(t(if (intmisc=="intens") x$qmodel$imatrix else x$emodel$imatrix))
        semat[ivector == 1] <- se; semat <- t(semat)
        lsemat[ivector == 1] <- lse; lsemat <- t(lsemat)
        inv.t <- if(intmisc=="intens") exp else plogis
        base.t <- if(intmisc=="intens") log else qlogis
        lmat <- inv.t(logest - qnorm(1 - 0.5*(1 - cl))*lsemat)
        umat <- inv.t(logest + qnorm(1 - 0.5*(1 - cl))*lsemat)
        imatrix <- if(intmisc=="intens") x$qmodel$imatrix else x$emodel$imatrix
        lmat[imatrix == 0] <- umat[imatrix == 0] <- 0
        ## SEs of diagonal entries
        diagse <- qematrix.diagse.msm(x, covariates, intmisc, sojourn,
                                      ni, ivector, nc, covlabels, covmeans, trsum)
        diag(semat) <- diagse$diagse
        diag(lmat) <- sign(diag(mat)[1]) * (inv.t(base.t(abs(diag(mat))) - sign(diag(mat)[1]) * qnorm(1 - 0.5*(1 - cl))*diagse$diaglse))
        diag(umat) <- sign(diag(mat)[1]) * (inv.t(base.t(abs(diag(mat))) + sign(diag(mat)[1]) * qnorm(1 - 0.5*(1 - cl))*diagse$diaglse))
        if (sojourn) {
            sojse <- diagse$sojse
            sojl <- inv.t(base.t(soj) - qnorm(1 - 0.5*(1 - cl))*diagse$sojlse)
            soju <- inv.t(base.t(soj) + qnorm(1 - 0.5*(1 - cl))*diagse$sojlse)
        }
        dimnames(semat) <- dimnames(x$qmodel$qmatrix)

    }
    else semat <- lmat <- umat <- sojse <- sojl <- soju <- NULL
    
    dimnames(mat) <-  dimnames(x$qmodel$qmatrix)
    if (sojourn)
      res <- list(estimates=mat, SE=semat, L=lmat, U=umat, sojourn=soj, sojournSE=sojse, sojournL=sojl, sojournU=soju)
    else
      res <- list(estimates=mat, SE=semat, L=lmat, U=umat)

    class(res) <- "msm.est"
    res
}

print.msm.est <- function(x, digits=NULL, ...)
  {
      print.ci(x$estimates, x$L, x$U, digits=digits)
  }

print.ci <- function(x, l, u, digits=NULL)
  { 
      est <- formatC(x, digits=digits)
      if (!is.null(l)) {
          low <- formatC(l, digits=digits)
          upp <- formatC(u, digits=digits)
          res <- paste(est, " (", low, ",", upp, ")", sep="")
          res[x==0] <- 0
      }
      else res <- est
      dim(res) <- dim(x)
      dimnames(res) <- dimnames(x)
      names(res) <- names(x)
      print(res, quote=FALSE)
  }

### Work out standard errors of diagonal entries of intensity/misc matrix, or sojourn times, using delta method

qematrix.diagse.msm <- function(x, covariates="mean", intmisc="intens", sojourn,
                                ni, ivector, nc, covlabels, covmeans, trsum)
  {
      nst <- x$qmodel$nstates
      diagse <- diaglse <- sojse <- sojlse <- numeric(nst)
      indmat <- matrix(ivector, nst, nst)
      indmat[indmat==1] <- seq(length = ni)
      indmat <- t(indmat) # matrix of indices of estimate vector 
      inds <- seq(length = ni + ni*nc)
      if (intmisc=="misc") inds <- c(x$hmodel$constr[x$hmodel$plabs=="p"], max(x$hmodel$constr) + x$hmodel$covconstr) +
        x$qmodel$npars + x$qcmodel$npars
      cur.i <- 1
      if (covariates == 0 && nc > 0)
        {covariates <- list();  for (i in 1 : nc) covariates[[covlabels[i]]] <- 0}
      coefs <- if (!is.list(covariates) && (covariates=="mean")) 1 else c(1, unlist(covariates) - covmeans)
      for (i in 1:nst){
          ## Transformation for delta method is
          ## exp(x1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... ) +
          ##  exp(x4 + x5 (cov1 - covmean1) + x6 (cov2 - covmean2) + ... ) +   (or expit(...)) 
          nir <- sum(indmat[i,-i] > 0) # number of intens/misc for current state
          if (nir > 0) {
              qf <- qematrix.diagse.formstr(nir, intmisc, inds, cur.i, ni, nc, coefs, 0, trsum)
              form <- as.formula(paste("~", paste(qf$formstr, collapse = " + ")))
              lform <- as.formula(paste("~ log (", paste(qf$formstr, collapse = " + "), ")"))
              ests <- x$estimates[qf$parinds2]
              cov <- x$covmat[qf$parinds2, qf$parinds2]
              diagse[i] <- deltamethod(form, ests, cov)
              diaglse[i] <- deltamethod(lform, ests, cov)
              if (sojourn){
                  ## Mean sojourn times are -1 / diagonal entries of q matrix. Calculate their SEs and CIs.
                  form <- as.formula(paste("~ 1 / (", paste(qf$formstr, collapse = " + "), ")"))
                  lform <- as.formula(paste("~ log ( 1 / (", paste(qf$formstr, collapse = " + "), ")", ")"))
                  sojse[i] <- deltamethod(form, ests, cov)
                  sojlse[i] <- deltamethod(lform, ests, cov)
              }
              cur.i <- cur.i + nir
          }
          else diagse[i] <- 0
      }
      list(diagse=diagse, diaglse=diaglse, sojse=sojse, sojlse=sojlse)
  }

### build a string for the delta-method formula for a diagonal intensity

qematrix.diagse.formstr <- function(nir, intmisc, inds, cur.i, ni, nc, coefs, offset, trsum)
{
    formstr <- character(nir)
    parinds <- numeric()      
    for (j in (cur.i : (cur.i + nir - 1))) {
        indj <-
          if (intmisc=="intens") seq(j, (nc * ni + j), ni)
          else c(j, seq(ni + (j-1)*nc + 1, length=nc))
        parinds <- c(parinds, inds[indj])  # first 1 5 7, then 2 5 7
    }
    ## e.g. parinds = 1 5 7 2 5 7 becomes xinds = 1 3 4 2 3 4
    parinds2 <- sort(unique(parinds))
    xinds <- rank(parinds2)[match(parinds, parinds2)] + offset
    for (j in 1:nir)
      formstr[j] <- trsum(xinds[1:(nc+1) + (j-1)*(nc+1)], coefs)
    list(formstr=formstr, parinds=parinds, parinds2=parinds2)
}

## Form a string,  exp ( x1 1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
## to be made into a formula for deltamethod

expsum <- function(inds, coefs)
{
    xseq <- paste("x", inds, sep="")
    inprod <- paste(paste(coefs, xseq, sep="*"), collapse=" + ")
    paste("exp(", inprod, ")", sep="")
}

## Form a string,  expit ( x1 1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
## to be made into a formula for deltamethod

expitsum <- function(inds, coefs)
{
    xseq <- paste("x", inds, sep="")
    inprod <- paste(paste(coefs, xseq, sep="*"), collapse=" + ")
    paste("exp(", inprod, ") / (1 + exp(", inprod, "))", sep="")
}

lsum <- function(inds, coefs)
  {
      xseq <- paste("x", inds, sep="")
      paste(paste(coefs, xseq, sep="*"), collapse=" + ")
  }


### Extract a ratio of transition intensities at given covariate values

qratio.msm <- function(x, ind1, ind2,
                       covariates = "mean", cl=0.95)
  {
      q <- qmatrix.msm(x, covariates)$estimates
      if (!is.numeric(ind1) || length(ind1) != 2 || !is.numeric(ind2) || length(ind2) != 2)
        stop("ind1 and ind2 must be numeric vectors of length 2")
      if (any (! (ind1 %in% 1 : x$qmodel$nstates))  |  any (! (ind2 %in% 1 : x$qmodel$nstates) ) )
        stop("ind1 and ind2 must be pairs of states in 1, ..., ", x$qmodel$nstates)
      if (q[ind2[1], ind2[2]] == 0)
        stop (paste("Denominator q[",ind2[1],",",ind2[2],"", "] is zero\n", sep=""))
      else if (q[ind1[1], ind1[2]] ==  0) {
          warning(paste ("Numerator q[",ind1[1],",",ind1[2],"", "] is zero\n", sep=""))
          estimate <- se <- 0
      }
      else {
          estimate <- q[ind1[1], ind1[2]]  /  q[ind2[1], ind2[2]]
          if (x$foundse) { 
              se <- qratio.se.msm(x, ind1, ind2, covariates, cl)$se 
              lse <- qratio.se.msm(x, ind1, ind2, covariates, cl)$lse
              L <- exp ( log(abs(estimate)) - sign(estimate)*qnorm(1 - 0.5*(1 - cl)) * lse ) * sign(estimate)
              U <- exp ( log(abs(estimate)) + sign(estimate)*qnorm(1 - 0.5*(1 - cl)) * lse ) * sign(estimate)
          }
          else {se <- L <- U <- NULL}
      }
      c(estimate=estimate, se=se, L=L, U=U)
  }

### Work out standard error of a ratio of intensities using delta method
### Uuugh.  What a fuss for one little number.

qratio.se.msm <- function(x, ind1, ind2, covariates="mean", cl=0.95)
  {
      nst <- x$qmodel$nstates
      ni <- x$qmodel$npars
      nc <- if (!is.list(covariates) && (covariates == "mean")) 0 else x$qcmodel$ncovs
      indmat <- t(x$qmodel$imatrix)
      indmat[indmat == 1] <- seq(length = x$qmodel$npars)
      indmat <- t(indmat) # matrix of indices of estimate vector 
      inds <- seq(length = x$qmodel$npars+x$qcmodel$npars) # identifiers for q and beta parameters
      if (covariates == 0)
        {covariates <- list();  for (i in 1 : nc) covariates[[x$qcmodel$covlabels[i]]] <- 0}
      coefs <- if (!is.list(covariates) && (covariates=="mean")) 1 else c(1, unlist(covariates) - x$qcmodel$covmeans)
      parinds <- numeric()
      indmatrow.n <- indmat[ind1[1],-ind1[1]]
      nir.n <- sum(indmatrow.n > 0)
      indmatrow.d <- indmat[ind2[1],-ind2[1]]
      nir.d <- sum(indmatrow.d > 0)
      formstr.n <- character(nir.n)
      formstr.d <- character(nir.d)

      if (ind1[1]!=ind1[2] && ind2[1]!=ind2[2]) { # both intensities are off-diagonal
          parinds <- c(inds[indmat[ind1[1],ind1[2]] - 1 + seq(1, (nc * ni + 1), ni)],
                       inds[indmat[ind2[1],ind2[2]] - 1 + seq(1, (nc * ni + 1), ni)])
          parinds2 <- sort(unique(parinds))
          xinds <- rank(parinds2)[match(parinds, parinds2)]
          formstr.n <- expsum(xinds[1:(nc+1)], coefs)
          formstr.d <- expsum(xinds[1:(nc+1) + nc+1] , coefs)
      }
      
      else if (ind1[1]!=ind1[2] && ind2[1]==ind2[2]) { # numerator off-diagonal, denom diagonal
          parinds <- inds[indmat[ind1[1],ind1[2]] - 1 + seq(1, (nc * ni + 1), ni)]
          cur.i <- min(indmatrow.d[indmatrow.d>0])
          for (j in 1:nir.d)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
          parinds2 <- sort(unique(parinds))
          xinds <- rank(parinds2)[match(parinds, parinds2)]
          formstr.n <- expsum(xinds[1:(nc+1)], coefs)
          for (j in 1:nir.d)
            formstr.d[j] <- expsum(xinds[1:(nc+1) + j*(nc+1)], coefs)
      }
      
      else if (ind1[1]==ind1[2] && ind2[1]!=ind2[2]) { # numerator diagonal, denom off-diagonal
          cur.i <- min(indmatrow.n[indmatrow.n>0])
          for (j in 1:nir.n)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
          parinds <- c(parinds, inds[indmat[ind2[1],ind2[2]] - 1 + seq(1, (nc * ni + 1), ni)])
          parinds2 <- sort(unique(parinds))
          xinds <- rank(parinds2)[match(parinds, parinds2)]
          for (j in 1:nir.n)
            formstr.n[j] <- expsum(xinds[1:(nc+1) + (j-1)*(nc+1)], coefs)
          formstr.d <- expsum(xinds[nir.n*(nc+1) + 1:(nc+1)], coefs)
      }
      
      else if (ind1[1]==ind1[2] && ind2[1]==ind2[2]) { # both intensities diagonal
          cur.i <- min(indmatrow.n[indmatrow.n>0])
          for (j in 1:nir.n)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
          cur.i <- min(indmatrow.d[indmatrow.d>0])
          for (j in 1:nir.d)
            parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])
          parinds2 <- sort(unique(parinds))
          xinds <- rank(parinds2)[match(parinds, parinds2)]
          for (j in 1:nir.n)
            formstr.n[j] <- expsum(xinds[1:(nc+1) + (j-1)*(nc+1)], coefs)
          for (j in 1:nir.d)
            formstr.d[j] <- expsum(xinds[nir.n*(nc+1) + 1:(nc+1) + (j-1)*(nc+1)], coefs)
      }

      num <- paste(formstr.n, collapse = " + ")
      denom <- paste(formstr.d, collapse = " + ")
      form <- as.formula(paste("~", "(", num, ") / (", denom, ")"))
      lform <- as.formula(paste("~ ", "log (", num, ") - log (", denom, ")"))
      ests <- x$estimates[parinds2]
      cov <- x$covmat[parinds2,parinds2]
      se <- deltamethod(form, ests, cov)
      lse <- deltamethod(lform, ests, cov)
      list(se=se, lse=lse)
  }


### Extract the transition probability matrix at given covariate values

pmatrix.msm <- function(x, # fitted msm model
                        t, # time interval
                        covariates = "mean"  # covariate values to calculate transition matrix for
                        )
  {
      if (!is.numeric(t) || (t < 0)) stop("t must be a positive number")
      q <- qmatrix.msm(x, covariates)
      p <- MatrixExp(q$estimates, t)
      colnames(p) <- rownames(p) <- rownames(q$estimates)
      p
  }

### Extract the transition probability matrix at given covariate values - where the Q matrix is piecewise-constant

pmatrix.piecewise.msm <- function(x, # fitted msm model
                                   t1, # start time
                                   t2, # stop time                                    
                                   times,  # vector of cut points 
                                   covariates # list of lists of covariates, for (, times1], (times1, times2], ...
                                        # of length one greater than times 
                                   )
  {
      ## Input checks
      if (t2 < t1) stop("Stop time t2 should be greater than or equal to start time t1")
      if (!is.numeric(times) || is.unsorted(times)) stop("times should be a vector of numbers in increasing order")
      if (length(covariates) != length(times) + 1)
        stop("Number of covariate lists must be one greater than the number of cut points")
      ## Locate which intervals t1 and t2 fall in, as indices ind1, ind2 into "times". 
      if (t1 <= times[1]) ind1 <- 1
      else if (length(times)==1) ind1 <- 2
      else {          
          for (i in 2:length(times))
            if ((t1 > times[i-1]) && (t1 <= times[i]))
              {ind1 <- i; break}
          if (t1 > times[i]) ind1 <- i+1
      }
      if (t2 <= times[1]) ind2 <- 1
      else if (length(times)==1) ind2 <- 2
      else {
          for (i in 2:length(times))
            if ((t2 > times[i-1]) && (t2 <= times[i]))
              {ind2 <- i; break}
          if (t2 > times[i]) ind2 <- i+1
      }
      
      ## Calculate accumulated pmatrix
      ## Three cases: ind1, ind2 in the same interval
      if (ind1 == ind2) {
          P <- pmatrix.msm(x, t2 - t1, covariates[[ind1]])
      }
      ## ind1, ind2 in successive intervals
      else if (ind2 == ind1 + 1) {
          P.start <- pmatrix.msm(x, times[ind1] - t1 , covariates[[ind1]])
          P.end <- pmatrix.msm(x, t2 - times[ind2-1], covariates[[ind2]])
          P <- P.start %*% P.end
      }
      ## ind1, ind2 separated by one or more whole intervals
      else {
          P.start <- pmatrix.msm(x, times[ind1] - t1 , covariates[[ind1]])
          P.end <- pmatrix.msm(x, t2 - times[ind2-1], covariates[[ind2]])
          P.middle <- diag(x$qmodel$nstates)
          for (i in (ind1+1):(ind2-1)) {
              P.middle <- P.middle %*% pmatrix.msm(x, times[i] - times[i-1], covariates[[i]])
          }
          P <- P.start %*% P.middle %*% P.end
      }
      
      P   
  }


### Extract the mean sojourn times for given covariate values

sojourn.msm <- function(x,
                        covariates = "mean", cl=0.95)
  {
      qmatrix <- qmatrix.msm(x, covariates, sojourn=TRUE, cl=cl)
      sojstates <- (1 : x$qmodel$nstates) [transient.msm(x)]
      soj <- qmatrix$sojourn[sojstates]
      names (soj) <- rownames(x$qmodel$qmatrix)[sojstates]
      if (x$foundse){
          sojse <- qmatrix$sojournSE[sojstates]
          sojl <- qmatrix$sojournL[sojstates]
          soju <- qmatrix$sojournU[sojstates]
          names(sojse) <- names(sojl) <- names(soju) <- names(soj)
          res <- data.frame(estimates=soj, SE=sojse, L=sojl, U=soju)
      }
      else res <- list(estimates=soj)
      res
  }

### Extract the coefficients

coef.msm <- function(object, ...)
  {
      if (!inherits(object, "msm")) stop("expected object to be a msm model")
      if (object$emodel$misc)
        object[c("Qmatrices", "Ematrices")]
      else object$Qmatrices
  }

### Extract the log-likelihood

logLik.msm <- function(object, ...)
  {
      if (!inherits(object, "msm")) stop("expected object to be a msm model")
      val <- 0.5 * object$minus2loglik
      attr(val, "df") <- length(object$estimates) - length(object$fixedpars)
      class(val) <- "logLik"
      val
  }

## Estimate total length of stay in a given state. 
## TODO: non-homogeneous models. 

totlos.msm <- function(x, start=1, fromt=0, tot=Inf, covariates="mean", ...)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    if (! start %in% 1 : x$qmodel$nstates) stop("start should be a state in 1, ..., ", x$qmodel$nstates)
    if (!is.numeric(fromt) || !is.numeric(tot) || length(fromt) != 1 || length(tot) != 1 || fromt < 0 || tot < 0)
      stop("fromt and tot must be single non-negative numbers")
    if (fromt > tot) stop("tot must be greater than fromt")
    if (length(absorbing.msm(x)) == 0)
      if (tot==Inf) stop("Must specify a finite end time for a model with no absorbing state")
    tr <- transient.msm(x)
    totlos <- numeric(length(tr))
    for (j in seq(along=tr)){
        f <- function(time) {
            y <- numeric(length(time))
            for (i in seq(along=y))
              y[i] <- pmatrix.msm(x, time[i], covariates)[start,tr[j]]
            y
        }        
        totlos[j] <- integrate(f, fromt, tot, ...)$value
    }
    names(totlos) <- rownames(x$qmodel$qmatrix)[tr]
    totlos
}

## Return indices of transient states (can either call for a fitted model or a qmatrix)

transient.msm <- function(x=NULL, qmatrix=NULL)
{
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- x$Qmatrices$baseline
        nst <- x$qmodel$nstates
    }
    else if (!is.null(qmatrix)) {
        nst <- nrow(qmatrix)
    }
    else stop("Neither a fitted msm model nor a qmatrix have been supplied") 
    (1:nst)[diag(qmatrix) < 0]
}

## Return indices of absorbing states (can either call for a fitted model or a qmatrix)

absorbing.msm <- function(x=NULL, qmatrix=NULL)
{
    if (!is.null(x)) {
        if (!inherits(x, "msm")) stop("expected x to be a msm model")
        qmatrix <- x$Qmatrices$baseline
        nst <- x$qmodel$nstates
    }
    else if (!is.null(qmatrix)) {
        nst <- nrow(qmatrix)
    }
    else stop("Neither a fitted msm model nor a qmatrix have been supplied") 
    (1:nst)[diag(qmatrix) == 0]
}


### Table of observed and expected prevalences (works for misclassification and non-misclassification models)

prevalence.msm <- function(x,
                           times=NULL,
                           timezero=NULL,
                           initstates=NULL,
                           covariates="mean",
                           misccovariates="mean"
                           )
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
### For general HMMs use the Viterbi estimate of the observed state.       
      state <- if (x$hmodel$hidden && !x$emodel$misc) viterbi.msm(x)$fitted else x$data$state
      nst <- x$qmodel$nstates
      time <- x$data$time
      if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)

### Estimate observed state occupancies in the data at a series of times
      cat("Calculating approximate observed state prevalences...\n")
      obs <- observed.msm(x, times)
      obstab <- obs$obstab; obsperc <- obs$obsperc; risk <- obs$risk
      
### Work out expected state occupancies by forecasting from transition probabilities
      if (is.null(timezero))  timezero <- min(time)
      cat("Forecasting expected state prevalences...\n")
      if (x$emodel$misc){
          exptab <- x$emodel$initprobs * risk[1]
          for (j in 2 : length(times) ) {
              pmat <- pmatrix.msm(x, times[j] - timezero, covariates=covariates)
              emat <- ematrix.msm(x, covariates=misccovariates)$estimates
              exptruej <- risk[j] * x$emodel$initprobs %*% pmat
              expobsj <- exptruej %*% emat
              exptab <- rbind(exptab, expobsj)
          }
      }
      else {
          if (is.null(initstates)){
              initstates <- table(state[time==timezero])
              z <- setdiff(paste(1:nst), names(initstates))
              zero <- rep(0, length(z)); names(zero) <- z; 
              initstates <- c(initstates, zero)
              initstates <- initstates[order(names(initstates))]
          }
          initprobs <- initstates / sum(initstates)
          exptab <- initstates
          for (j in 2 : length(times) ) {
              pmat <- pmatrix.msm(x, times[j] - timezero, covariates=covariates)
              expj <- risk[j] * initprobs %*% pmat
              exptab <- rbind(exptab, expj)
          }
      }
      exptab <- cbind(exptab, apply(exptab, 1, sum))
      dimnames(exptab) <- list(times, c(rownames(x$qmodel$qmatrix),"Total"))
      expperc <- 100*exptab[,1:nst] / exptab[, nst+1]

      res <- list(observed=obstab, expected=exptab, obsperc=obsperc, expperc=expperc)
      names(res) <- c("Observed", "Expected", "Observed percentages", "Expected percentages")
      res      

  }


### Estimate observed state occupancies in the data at a series of times
### Assume previous observed state is retained until next observation time
### This is rather slow, especially for lots of patients, and should probably be
### rewritten in C. 

observed.msm <- function(x, times=NULL)
{
    if (!inherits(x, "msm")) stop("expected x to be a msm model")
    state <- if (x$hmodel$hidden && !x$emodel$misc) viterbi.msm(x)$fitted else x$data$state
    time <- x$data$time
    subject <- x$data$subject
    nst <- x$qmodel$nstates
    if (is.null(times))
      times <- seq(min(time), max(time), (max(time) - min(time))/10)
    getcontrib <- function(pt, subject, time, state, times, nst, absorb){
        rows <- subject==pt
        contrib <- matrix(0, nrow=length(times), ncol=nst)
        risk <- rep(0, length(times))
        for (i in seq(along=times)){
            t <- times[i]
            endtime <- max(time[rows])
            endstate <- state[rows][length(state[rows])]
            if ( t <= endtime ) {
                risk[i] <-  1
                for (j in seq(length(time[rows])-1)){
                    if ( ( (time[rows][j+1] > t) & (time[rows][j] <= t)) |
                        ( (time[rows][j+1] == t) & (t == endtime) ) )
                      {
                          currstate <- state[rows][j]
                          contrib[i, currstate] <- 1
                      }
                }
            }
            else if (t >  endtime  &  endstate %in% absorb) {
                risk[i] <- 1
                contrib[i, endstate] <- 1
            }
        }
        list(risk, contrib)
    }
    obstab <- matrix(0, nrow=length(times), ncol=nst)
    risk <- rep(0, length(times))
    for (pt in unique(as.numeric(subject))){
        cont <- getcontrib(pt, as.numeric(subject), time, state, times, nst, absorbing.msm(x))
        risk <- risk + cont[[1]]
        obstab <- obstab + cont[[2]]
    }
    obstab <- cbind(obstab, apply(obstab, 1, sum))
    dimnames(obstab) <- list(times, c(rownames(x$qmodel$qmatrix), "Total"))
    obsperc <- 100*obstab[,1:nst] / obstab[, nst+1]
    list(obstab=obstab, obsperc=obsperc, risk=risk)
}



### Obtain hazard ratios from estimated effects of covariates on log-transition rates

hazard.msm <- function(x, hazard.scale = 1, cl = 0.95)
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
      if (length(hazard.scale) == 1) hazard.scale <- rep(hazard.scale, x$qcmodel$ncovs)
      if (length(hazard.scale) != x$qcmodel$ncovs)
        stop ("hazard.scale of length ", length(hazard.scale), ", expected ", x$qcmodel$ncovs)
      keep <- (x$qmodel$imatrix != 0)
      nst <- x$qmodel$nstates
      keepvec <- as.vector(keep)      
      fromlabs <- rep(rownames(keep), nst) [keepvec] 
      tolabs <- rep(colnames(keep), rep(nst, nst)) [keepvec] 
      if (x$qcmodel$ncovs > 0) {
          haz.list <- list()
          if (x$foundse) {
              for (i in 1:x$qcmodel$ncovs) {
                  cov <- x$qcmodel$covlabels[i]
                  haz.rat <- exp(hazard.scale[i]*x$Qmatrices[[cov]])[keepvec]
                  l95 <- exp(hazard.scale[i]*(x$Qmatrices[[cov]] - qnorm(1 - 0.5*(1 - cl))*x$QmatricesSE[[cov]]) )[keepvec]
                  u95 <- exp(hazard.scale[i]*(x$Qmatrices[[cov]] + qnorm(1 - 0.5*(1 - cl))*x$QmatricesSE[[cov]]) )[keepvec]
                  haz.tab <- cbind(haz.rat, l95, u95)
                  dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs),
                                            c("HR", "L", "U"))
                  haz.list[[cov]] <- haz.tab
              }
          }
          else {
              for (i in 1:x$qcmodel$ncovs) {
                  cov <- x$qcmodel$covlabels[i]
                  haz.tab <- as.matrix(exp(hazard.scale[i]*x$Qmatrices[[cov]])[keepvec])
                  dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs), "HR")
                  haz.list[[cov]] <- haz.tab
              }
          }
      }
      else haz.list <- "No covariates on transition intensities"
      haz.list
  }


### Obtain odds ratios from estimated effects of covariates on logit-misclassification probabilities
### TODO - equivalent for general HMMs which presents cov effects on natural scale.

odds.msm <- function(x, odds.scale = 1, cl = 0.95)
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
      if (!x$emodel$misc) stop("Requires a misclassification model specified with ematrix")
      if (length(odds.scale) == 1) odds.scale <- rep(odds.scale, x$ecmodel$ncovs)
      if (length(odds.scale) != x$ecmodel$ncovs)
        stop ("odds.scale of length ", length(odds.scale), ", expected ", x$ecmodel$ncovs)
      keep <- (x$emodel$imatrix != 0)
      nst <- x$qmodel$nstates
      keepvec <- as.vector(keep)
      truelabs <- rep(rownames(keep), nst) [keepvec] 
      obslabs <- rep(colnames(keep), rep(nst, nst)) [keepvec] 
      if (x$ecmodel$ncovs > 0) {
          odds.list <- list()
          if (x$foundse) {
              for (i in 1:x$ecmodel$ncovs) {
                  cov <- x$ecmodel$covlabels[i]
                  odds.rat <- exp(odds.scale[i]*x$Ematrices[[cov]])[keepvec]
                  l95 <- exp(odds.scale[i]*(x$Ematrices[[cov]] - qnorm(1 - 0.5*(1 - cl))*x$EmatricesSE[[cov]]) )[keepvec]
                  u95 <- exp(odds.scale[i]*(x$Ematrices[[cov]] + qnorm(1 - 0.5*(1 - cl))*x$EmatricesSE[[cov]]) )[keepvec]
                  odds.tab <- cbind(odds.rat, l95, u95)
                  dimnames(odds.tab) <- list(paste("Obs", obslabs, "|", truelabs),
                                             c("OR", "L", "U"))
                  odds.list[[cov]] <- odds.tab
              }
          }
          else {
              for (i in 1:x$ecmodel$ncovs) {
                  cov <- x$ecmodel$covlabels[i]
                  odds.tab <- as.matrix(exp(odds.scale[i]*x$Ematrices[[cov]])[keepvec])
                  dimnames(odds.tab) <- list(paste("Obs", obslabs, "|", truelabs), "OR")
                  odds.list[[cov]] <- odds.tab
              }
          }
      }
      else odds.list <- "No covariates on misclassification probabilities"
      odds.list
  }



### Viterbi algorithm for reconstructing the most likely path through underlying states
### This is all done in C

viterbi.msm <- function(x)
  {
      if (!inherits(x, "msm")) stop("expected x to be a msm model")
      if (x$hmodel$hidden) 
        {
            do.what <- 2

            x$data$fromstate <- x$data$fromstate - 1  
            x$data$tostate <- x$data$tostate - 1
            x$data$firstobs <- x$data$firstobs - 1
            x$hmodel$models <- x$hmodel$models - 1
            x$hmodel$links <- x$hmodel$links - 1

            vit <- .C("msmCEntry",                
                      as.integer(do.what),

                      as.integer(as.vector(t(x$qmodel$imatrix))),
                      as.double(as.vector(t(x$Qmatrices$baseline)[t(x$qmodel$imatrix)==1])),
                      as.double(unlist(lapply(x$Qmatrices[x$qcmodel$covlabels], function(x) t(x)[t(misc.msm$qmodel$imatrix)==1]))),
                      as.double(x$hmodel$pars),
                      as.double(x$hmodel$coveffect),   # estimates

                      ## data for non-HMM
                      as.integer(x$data$fromstate),
                      as.integer(x$data$tostate),
                      as.double(x$data$timelag),
                      as.double(unlist(x$data$covmat)),
                      as.integer(x$data$covdata$whichcov), # this is really part of the model 
                      as.integer(x$data$nocc),
                      as.integer(x$data$whicha),
                      as.integer(x$data$obstype),

                      ## data for HMM / censored
                      as.integer(match(x$data$subject, unique(x$data$subject))), 
                      as.double(x$data$time),
                      as.double(x$data$state),
                      as.integer(x$data$firstobs),
                      
                      ## HMM specification
                      as.integer(x$hmodel$hidden),
                      as.integer(x$hmodel$models),
                      as.integer(x$hmodel$npars),
                      as.integer(x$hmodel$totpars),
                      as.integer(x$hmodel$firstpar),
                      as.integer(x$hmodel$ncovs),
                      as.integer(x$hmodel$whichcovh),
                      as.integer(x$hmodel$links),                
                      as.double(x$hmodel$initprobs),

                      ## various dimensions
                      as.integer(x$qmodel$nstates),
                      as.integer(x$qmodel$npars),
                      as.integer(x$data$nobs),
                      as.integer(x$data$npts),  # HMM only
                      as.integer(rep(x$qcmodel$ncovs, x$qmodel$npars)),
                      
                      as.integer(x$cmodel$ncens),
                      as.integer(x$cmodel$censor - 1),
                      as.integer(x$cmodel$states - 1),
                      as.integer(x$cmodel$index - 1),

                      fitted = double (x$data$nobs),
                      PACKAGE = "msm",
                      NAOK = TRUE
                      )
            
            fitted <- vit$fitted + 1
        }
      else fitted <- x$data$state
      data.frame(subject = x$data$subject,  
                 time = x$data$time,
                 observed = x$data$state,
                 fitted = fitted)
  }
