### METHODS FOR MSM OBJECTS

print.msm <- function(x,...)
  {
      cat ("\n Multi-state Markov models in continuous time \n\n")

      printmatlist <- function(x, cpt, matrixtype, transformtype) {
          M <- x[[cpt]]; MSE <- x[[paste(cpt,"SE",sep="")]]
          base <- if (is.list(M)) M[[1]] else M
          nmatrix <- dim(base)[1]
          covmessage <- if(length(M) == 2) "" else "with covariates set to their means"
          cat ("  * Matrix of", matrixtype, covmessage," \n\n")
          print (M[["baseline"]], na.print="."); cat("\n")
          cat ("    corresponding standard errors \n\n")
          ses <- MSE[["baseline"]]
          print (ses, na.print="."); cat("\n")
          covlabels <- if(cpt=="Qmatrices") x$data$covlabels else x$data$misccovlabels
          if (length(M) == 2)
            cat("  * No covariates on", matrixtype, "\n\n")
          else {
              for (cov in covlabels) {
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
          printmatlist (x, "Qmatrices", "transition intensities", "log")
          if (x$misc)
            printmatlist (x, "Ematrices", "misclassification probabilities", "logit")
          covmessage <- if (x$data$ncovs == 0) "" else "with covariates set to their means"
          cat ("  * Mean sojourn times in transient states", covmessage, "\n\n")
          print( x$sojourn ); cat("\n")
      }

      cat (" -2 * log-likelihood: ", x$minus2loglik, "\n")
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
      cat("Calculating tables of observed and expected stage occupancy...\n")
      prevalences <- prevalence.msm(object, times, timezero, initstates, covariates, misccovariates)
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
          print(x$prevalences$Observed)
          cat("\nExpected numbers of individuals occupying stages at each time\n\n")
          print(x$prevalences$Expected)
          cat("\nObserved prevalences of stages (percentages of population at risk)\n\n")
          print(x$prevalences$"Observed percentages")
          cat("\nExpected prevalences of stages (percentages of population at risk)\n\n")
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

plot.msm <- function(x, from=NULL, to=NULL, covariates="mean", legend.pos=NULL, ...)
  {
      if (is.null(to))
        to <- max(absorbing.msm(x))
      if (is.null(from))
        from <- transient.msm(x)
      rg <- range(x$data$time)
      timediff <- (rg[2] - rg[1]) / 50
      times <- seq(rg[1], rg[2], timediff)
      pr <- numeric()
      for (t in times)
        pr <- c(pr, pmatrix.msm(x, t, covariates)[from[1], to])
      plot(times, 1 - pr, type="l", xlab="Time", ylab="Fitted survival probability",
           ylim=c(0,1), lty = 1, ...)
      lt <- 2
      for (st in from[-1]){ 
          pr <- numeric()
          for (t in times)
            pr <- c(pr, pmatrix.msm(x, t, covariates)[st, to])
          lines(times, 1 - pr, type="l", lty = lt)
          lt <- lt+1
      }
      if (length(legend.pos) != 2)
        legend.pos <- c(max(times) - 15*timediff, 1)
      legend(legend.pos[1], legend.pos[2], legend=paste("From stage",from), lty = seq(lt-1))
      invisible()
  }


### Extract the transition intensity matrix at given covariate values

qmatrix.msm <- function(x, # fitted msm model
                        covariates = "mean",  # covariate values to calculate transition matrix for
                        sojourn = FALSE      # also calculate mean sojourn times and their standard errors
                        )
  {
      qematrix.msm(x, covariates, which="intens", sojourn)
  }

### Extract the misclassification probability matrix at given covariate values

ematrix.msm <- function(x,
                        covariates = "mean")
  {
      if (!x$misc)
        "Model without misclassification"
      else
        qematrix.msm(x, covariates, which="misc", sojourn=FALSE)
  }

### Extract either intensity or misclassification matrix

qematrix.msm <- function(x, covariates="mean", which="intens", sojourn=FALSE)
{
    nst <- x$model$nstates
    if (which=="intens"){
        ni <- x$model$nintens
        nieffs <- x$model$nintenseffs
        nc <- if (!is.list(covariates) && (covariates == "mean")) 0 else x$data$ncovs
        if (nc==0){
            mat <- exp(x$Qmatrices$logbaseline)
            qmatrix.ind <- t(matrix(x$model$qvector, nrow=x$model$nstates))
            mat[qmatrix.ind == 0] <- 0
            diag(mat) <- 0;  diag(mat) <- - apply(mat, 1, sum)
        }
        constrvec <- x$model$constrvec
        baseconstrvec <- x$model$baseconstrvec
        covlabels <- x$data$covlabels
        covmeans <- x$data$covmeans
    }
    else if (which=="misc"){
        ni <- x$model$nmisc
        nieffs <- x$model$nmisceffs
        nc <- if (!is.list(covariates) && (covariates == "mean")) 0 else x$data$nmisccovs
        if (nc==0){
            mat <- expit(x$Ematrices$logitbaseline)
            ematrix.ind <- t(matrix(x$model$evector, nrow=x$model$nstates))
            mat[ematrix.ind == 0] <- 0
            diag(mat) <- 0;  diag(mat) <- 1 - apply(mat, 1, sum)
        }
        constrvec <- x$model$miscconstrvec
        baseconstrvec <- x$model$basemiscconstrvec
        covlabels <- x$data$misccovlabels
        covmeans <- x$data$misccovmeans
    }

    se <- numeric(ni)
    if ((nc == 0) && (is.list(covariates)) && (length(covariates) > 0))
      warning(paste("Ignoring covariates - no covariates in model for",
                    if (which=="intens") "transition intensities" else "misclassification probabilities"))
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
        if (which=="intens"){
            logq <- x$Qmatrices$logbaseline
            for (i in 1 : nc)
              logq <- logq + x$Qmatrices[[i+1]] * (covariates[[i]] - covmeans[i])
            mat <- exp(logq)
            qmatrix.ind <- t(matrix(x$model$qvector, nrow=x$model$nstates))
            mat[qmatrix.ind == 0] <- 0
            diag(mat) <- 0;  diag(mat) <- - apply(mat, 1, sum)
        }
        else if (which=="misc"){
            logite <- x$Ematrices$logitbaseline
            for (i in 1 : nc)
              logite <- logite + x$Ematrices[[i+1]] * (covariates[[i]] - covmeans[i])
            mat <- expit(logite)
            ematrix.ind <- t(matrix(x$model$evector, nrow=x$model$nstates))
            mat[ematrix.ind == 0] <- 0
            diag(mat) <- 0;  diag(mat) <- 1 - apply(mat, 1, sum)
        }        
    }
    if (sojourn) soj = -1 / diag(mat)
    
    if (x$foundse) {
        ## Work out standard errors. 
        ## Transformation for delta method is
        ##  exp (or expit) (x1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... )
        coefs <- if (!is.list(covariates) && (covariates=="mean")) 1 else c(1, unlist(covariates) - covmeans)
        trsum <- if (which=="intens") expsum else expitsum
        form <- as.formula(paste("~", trsum(seq(nc + 1), coefs)))
        inds <-
          if (which=="intens")
            c(baseconstrvec, x$model$nintenseffs + constrvec)
          else
            x$model$nintenseffs + x$model$ncoveffs + c(baseconstrvec, x$model$nmisceffs + constrvec) 
        for (i in 1 : ni){
            parinds <- inds[seq(i, (nc * ni + i), ni)]
            ests <- x$estimates[parinds]
            cov <- x$covmat[parinds, parinds]
            se[i] <- deltamethod(form, ests, cov)
        }
        
        semat <- matrix(0, nst, nst)
        ivector <- if (which=="intens") x$model$qvector else x$model$evector
        semat[ivector == 1] <- se; semat <- t(semat)
        ## SEs of diagonal entries
        diagse <- qematrix.diagse.msm(x, covariates, which, sojourn,
                                      ni, nieffs, ivector, nc, constrvec, baseconstrvec, covlabels, covmeans, trsum)
        diag(semat) <- diagse$diagse
        sojse <- diagse$sojse
        dimnames(semat) <- list(paste("Stage", 1:nst), paste("Stage", 1:nst))
    }
    else semat <- sojse <- NULL
    
    dimnames(mat) <-  list(paste("Stage", 1:nst), paste("Stage", 1:nst))
    if (sojourn)
      list(estimates=mat, SE=semat, sojourn=soj, sojournSE=sojse)
    else
      list(estimates=mat, SE=semat)
    
}

### Work out standard errors of diagonal entries of intensity/misc matrix, or sojourn times, using delta method

qematrix.diagse.msm <- function(x, covariates="mean", which="intens", sojourn,
                                ni, nieffs, ivector, nc, constrvec, baseconstrvec, covlabels, covmeans, trsum)
  {
      nst <- x$model$nstates
      diagse <- sojse <- numeric(nst)
      indmat <- matrix(ivector, nst, nst)
      indmat[indmat==1] <- baseconstrvec
      indmat <- t(indmat) # matrix of indices of estimate vector 
      inds <- c(baseconstrvec, nieffs + constrvec)
      if (which=="misc") inds <- inds + x$model$nintenseffs + x$model$ncoveffs
      cur.i <- 1
      if (covariates == 0 && nc > 0)
        {covariates <- list();  for (i in 1 : nc) covariates[[covlabels[i]]] <- 0}
      coefs <- if (!is.list(covariates) && (covariates=="mean")) 1 else c(1, unlist(covariates) - covmeans)
      for (i in 1:nst){
          ## Transformation for delta method is
          ## exp(x1 + x2 (cov1 - covmean1) + x3 (cov2 - covmean2) + ... ) +
          ##  exp(x4 + x5 (cov1 - covmean1) + x6 (cov2 - covmean2) + ... ) +   (or expit(...)) 
          nir <- sum(indmat[i,-i] > 0)
          if (nir > 0) {
              qf <- qematrix.diagse.formstr(nir, inds, cur.i, ni, nc, coefs, 0, trsum)
              form <- as.formula(paste("~", paste(qf$formstr, collapse = " + ")))
              ests <- x$estimates[qf$parinds2]
              cov <- x$covmat[qf$parinds2, qf$parinds2]
              diagse[i] <- deltamethod(form, ests, cov)
              if (sojourn){
                  ## Mean sojourn times are -1 / diagonal entries of q matrix. Calculate their SEs.
                  form <- as.formula(paste("~ 1 / (", paste(qf$formstr, collapse = " + "), ")"))
                  sojse[i] <- deltamethod(form, ests, cov)
              }
              cur.i <- cur.i + nir
          }
          else diagse[i] <- 0
      }
      list(diagse=diagse, sojse=sojse)
  }

### build a string for the delta-method formula for a diagonal intensity

qematrix.diagse.formstr <- function(nir, inds, cur.i, ni, nc, coefs, offset, trsum)
{
    formstr <- character(nir)
    parinds <- numeric()
    for (j in 1:nir)
      parinds <- c(parinds, inds[cur.i - 1 + seq(j, (nc * ni + j), ni)])  # first 1 5 7, then 2 5 7
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

### Extract a ratio of transition intensities at given covariate values

qratio.msm <- function(x, ind1, ind2,
                       covariates = "mean")
  {
      q <- qmatrix.msm(x, covariates)$estimate
      if (q[ind2[1], ind2[2]] == 0)
        stop (paste("Denominator q[",ind2[1],",",ind2[2],"", "] is zero\n", sep=""))
      else if (q[ind1[1], ind1[2]] ==  0) {
          warning(paste ("Numerator q[",ind1[1],",",ind1[2],"", "] is zero\n", sep=""))
          estimate <- se <- 0
      }
      else {
          estimate <- q[ind1[1], ind1[2]]  /  q[ind2[1], ind2[2]]
          se <- if (x$foundse) qratio.se.msm(x, ind1, ind2, covariates) else NULL
      }
      list(estimate=estimate, se=se)
  }

### Work out standard error of a ratio of intensities using delta method
### Uuugh.  What a fuss for one little number. 

qratio.se.msm <- function(x, ind1, ind2, covariates="mean")
  {
      nst <- x$model$nstates
      ni <- x$model$nintens
      nieffs <- x$model$nintenseffs
      nc <- if (!is.list(covariates) && (covariates == "mean")) 0 else x$data$ncovs
      indmat <- matrix(x$model$qvector, nst, nst)
      indmat[indmat==1] <- x$model$baseconstrvec
      indmat <- t(indmat) # matrix of indices of estimate vector 
      inds <- c(x$model$baseconstrvec, nieffs + x$model$constrvec) # identifiers for q and beta parameters
      if (covariates == 0)
        {covariates <- list();  for (i in 1 : nc) covariates[[x$data$covlabels[i]]] <- 0}
      coefs <- if (!is.list(covariates) && (covariates=="mean")) 1 else c(1, unlist(covariates) - x$data$covmeans)
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
      ests <- x$estimates[parinds2]
      cov <- x$covmat[parinds2,parinds2]
      se <- deltamethod(form, ests, cov)
      se
  }


### Extract the transition probability matrix at given covariate values

pmatrix.msm <- function(x, # fitted msm model
                        t, # time interval
                        covariates = "mean"  # covariate values to calculate transition matrix for
                        )
  {
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
      if (t2 < t1) stop("Stop time t2 earlier than start time t1")
      if (length(covariates) != length(times) + 1)
        stop("Number of covariate lists must be one greater than the number of cut points")
      ## Locate which intervals t1 and t2 fall in, as indices ind1, ind2 into "times". 
      if (t1 <= times[1]) ind1 <- 1
      else {
          for (i in 2:length(times))
            if ((t1 > times[i-1]) && (t1 <= times[i]))
              {ind1 <- i; break}
          if (t1 > times[i]) ind1 <- i+1
      }
      if (t2 <= times[1]) ind2 <- 1
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
          P.middle <- diag(x$model$nstates)
          for (i in (ind1+1):(ind2-1)) {
              P.middle <- P.middle %*% pmatrix.msm(x, times[i] - times[i-1], covariates[[i]])
          }
          P <- P.start %*% P.middle %*% P.end
      }
      
      P   
  }


### Extract the mean sojourn times for given covariate values (NEW VERSION)

sojourn.msm <- function(x,
                        covariates = "mean")
  {
      qmatrix <- qmatrix.msm(x, covariates, sojourn=TRUE)
      sojstages <- (1 : x$model$nstates) [transient.msm(x)]
      soj <- qmatrix$sojourn[sojstages]
      names (soj) <- paste("Stage",sojstages)
      if (x$foundse){
          sesoj <- qmatrix$sojournSE[sojstages]
          names (sesoj) <- paste("Stage",sojstages)
      }
      else sesoj <- NULL
      list(estimate=soj, SE=sesoj)
  }

### Extract the coefficients

coef.msm <- function(object, ...)
  {
      if (object$misc)
        object[c("Qmatrices", "Ematrices")]
      else object$Qmatrices
  }


## Estimate total length of stay in a given state. 
## TODO: non-homogeneous models. 

totlos.msm <- function(x, start=1, fromt=0, tot=Inf, covariates="mean", ...)
{
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
    names(totlos) <- paste("Stage",tr)
    totlos
}

## Return indices of transient states

transient.msm <- function(x)
{
    q <- x$Qmatrices$baseline
    nst <- x$model$nstates
    seq(nst)[diag(q) < 0]
}

## Return indices of absorbing states

absorbing.msm <- function(x)
{
    q <- x$Qmatrices$baseline
    nst <- x$model$nstates
    seq(nst)[diag(q) == 0]
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
      if (x$data$fromto) stop("Need observation time data for a prevalence table (not \"fromstate-tostate\" style)")
      time <- x$data$time
      nst <- x$model$nstates
      if (is.null(times))
        times <- seq(min(time), max(time), (max(time) - min(time))/10)

### Estimate observed state occupancies in the data at a series of times
      cat("Calculating approximate observed state prevalences...\n")
      obs <- observed.msm(x, times)
      obstab <- obs$obstab; obsperc <- obs$obsperc; risk <- obs$risk
      
### Work out expected state occupancies by forecasting from transition probabilities
      if (is.null(timezero))  timezero <- min(time)
      cat("Forecasting expected state prevalences...\n")
      if (x$misc){
          exptab <- x$model$initprobs * risk[1]
          for (j in 2 : length(times) ) {
              pmat <- pmatrix.msm(x, times[j] - timezero, covariates=covariates)
              emat <- ematrix.msm(x, covariates=misccovariates)$estimates
              exptruej <- risk[j] * x$model$initprobs %*% pmat
              expobsj <- exptruej %*% emat
              exptab <- rbind(exptab, expobsj)
          }
      }
      else {
          if (is.null(initstates)){
              initstates <- table(x$data$state[time==timezero])
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
      dimnames(exptab) <- list(times, c(paste("Stage",1:nst),"Total"))
      expperc <- 100*exptab[,1:nst] / exptab[, nst+1]

      res <- list(observed=obstab, expected=exptab, obsperc=obsperc, expperc=expperc)
      names(res) <- c("Observed", "Expected", "Observed percentages", "Expected percentages")
      res      

  }


### Estimate observed state occupancies in the data at a series of times
### Assume previous observed state is retained until next observation time
### This is rather slow, especially for lots of patients, and should probably be
### rewritten in C. 

observed.msm <- function(x, times)
{
    time <- x$data$time; state <- x$data$state; subject <- x$data$subject
    nst <- x$model$nstates
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
    dimnames(obstab) <- list(times, c(paste("Stage",1:nst), "Total"))
    obsperc <- 100*obstab[,1:nst] / obstab[, nst+1]
    list(obstab=obstab, obsperc=obsperc, risk=risk)
}



### Obtain hazard ratios from estimated effects of covariates on log-transition rates

hazard.msm <- function(x, hazard.scale = 1)
  {
      keep <- (x$Qmatrices[["logbaseline"]] != 0)
      nst <- x$model$nstates
      keepvec <- as.vector(keep)
      fromlabs <- rep(row.names(keep), nst) [keepvec]
      tolabs <- rep(dimnames(keep)[[2]], rep(nst, nst)) [keepvec]
      if (x$data$ncovs > 0) {
          haz.list <- list()
          if (x$foundse) {
              for (cov in x$data$covlabels) {
                  haz.rat <- exp(hazard.scale*x$Qmatrices[[cov]])[keepvec]
                  l95 <- exp(hazard.scale*(x$Qmatrices[[cov]] - qnorm(0.975)*x$QmatricesSE[[cov]]) )[keepvec]
                  u95 <- exp(hazard.scale*(x$Qmatrices[[cov]] + qnorm(0.975)*x$QmatricesSE[[cov]]) )[keepvec]
                  haz.tab <- cbind(haz.rat, l95, u95)
                  dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs),
                                            c("HR", "L95", "U95"))
                  haz.list[[cov]] <- haz.tab
              }
          }
          else {
              for (cov in x$data$covlabels) {
                  haz.tab <- as.matrix(exp(hazard.scale*x$Qmatrices[[cov]])[keepvec])
                  dimnames(haz.tab) <- list(paste(fromlabs, "-", tolabs), "HR")
                  haz.list[[cov]] <- haz.tab
              }
          }
      }
      else haz.list <- "No covariates on transition intensities"
      haz.list
  }


### Obtain odds ratios from estimated effects of covariates on logit-misclassification probabilities

odds.msm <- function(x, odds.scale = 1)
  {
      if (!x$misc) stop("Not a misclassification model")
      keep <- (x$Ematrices[["logitbaseline"]] != 0)
      nst <- x$model$nstates
      keepvec <- as.vector(keep)
      truelabs <- rep(row.names(keep), nst) [keepvec]
      obslabs <- rep(dimnames(keep)[[2]], rep(nst, nst)) [keepvec]
      if (x$data$nmisccovs > 0) {
          odds.list <- list()
          if (x$foundse) {
              for (cov in x$data$misccovlabels) {
                  odds.rat <- exp(odds.scale*x$Ematrices[[cov]])[keepvec]
                  l95 <- exp(odds.scale*(x$Ematrices[[cov]] - qnorm(0.975)*x$EmatricesSE[[cov]]) )[keepvec]
                  u95 <- exp(odds.scale*(x$Ematrices[[cov]] + qnorm(0.975)*x$EmatricesSE[[cov]]) )[keepvec]
                  odds.tab <- cbind(odds.rat, l95, u95)
                  dimnames(odds.tab) <- list(paste("Obs", obslabs, "|", truelabs),
                                             c("OR", "L95", "U95"))
                  odds.list[[cov]] <- odds.tab
              }
          }
          else {
              for (cov in x$data$misccovlabels) {
                  odds.tab <- as.matrix(exp(odds.scale*x$Ematrices[[cov]])[keepvec])
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
      if (!x$misc) stop("Viterbi algorithm is for models with misclassification")
      ematrix <- t(matrix(x$model$evector, nrow=x$model$nstates))
      esum <- apply(ematrix, 1, sum)
      nms <- max ( seq(along=esum)[esum > 0])
      do.what <- 2
      vit <- .C("msmCEntry",
                as.integer(do.what),
                as.double (x$estimates.t),
                as.double (x$estimates.t),
                as.integer(1),
                as.integer(length(x$estimates.t)),
                as.integer (x$data$subject),
                as.double (x$data$time),
                as.integer (x$data$state - 1),
                as.integer(NULL),
                as.integer(1),
                as.integer(x$model$qvector),
                as.integer(x$model$evector),
                as.double (x$data$covvec),
                as.integer (x$model$constrvec),
                as.double (x$data$misccovvec),
                as.integer(x$model$miscconstrvec),
                as.integer(x$model$baseconstrvec),
                as.integer(x$model$basemiscconstrvec),
                as.double (x$model$initprobs),
                as.integer (x$model$nstates),
                as.integer (nms),
                as.integer (x$model$nintens),
                as.integer (x$model$nintenseffs),
                as.integer (x$model$nmisc),
                as.integer (x$model$nmisceffs),
                as.integer (x$data$nobs),
                as.integer (x$data$npts),
                as.integer (x$data$ncovs),
                as.integer (x$model$ncoveffs),
                as.integer (x$data$nmisccovs),
                as.integer (x$model$nmisccoveffs),
                as.integer (x$model$covmatch),
                as.integer (x$model$ndeath),                
                as.integer (x$model$death),
                as.integer (x$model$exacttimes),
                as.integer(0),
                as.integer(-1),
                fitted = double (x$data$nobs),
                PACKAGE = "msm"
                )
      
      fitted <- vit$fitted + 1
      data.frame(subject = x$data$subject,
                 time = x$data$time,
                 observed = x$data$state,
                 fitted = fitted)
  }
