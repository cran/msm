### Take a bootstrap sample from the data contained in a fitted msm
### model. Sample pairs of consecutive observations, i.e. independent
### transitions.  Not applicable if model is hidden or some states are
### censored.

bootdata.trans.msm <- function(x) {
  subj.num <- match(x$data$subject, unique(x$data$subject))
  nextsubj <- c(subj.num[2:length(subj.num)], Inf)
  lastsubj <- subj.num != nextsubj
  inds <- sample(which(!lastsubj), replace=TRUE)
  data.boot <- matrix(nrow=length(inds)*2, ncol=length(x$data$covlabels) + 3)
  subj.name <- deparse(x$call$subject)
  state.name <- deparse(as.list(x$call$formula)[[2]])
  time.name <- deparse(as.list(x$call$formula)[[3]])
  colnames(data.boot) <- c(subj.name, time.name, state.name, x$data$covlabels)
  data.boot[,state.name] <- as.vector(rbind(x$data$state[inds], x$data$state[inds+1]))
  # in the bootstrap data, label each transition as being from a different subject
  data.boot[,subj.name] <- rep(seq(along=inds), each=2)
  data.boot[,time.name] <- as.vector(rbind(x$data$time[inds], x$data$time[inds+1]))
  for (j in x$data$covlabels) {
    data.boot[seq(1, 2*length(inds)-1, 2), j] <- x$data$cov[inds, j] + x$data$covdata$covmeans[j]
    data.boot[seq(2, 2*length(inds), 2), j] <- 0 # this is ignored
  }
  as.data.frame(data.boot)
}

### Take a bootstrap sample from the data contained in a fitted msm
### model. Sample subjects. Used for hidden models or models with
### censoring, in which the transitions within a subject are not
### independent.

bootdata.subject.msm <- function(x) {
  subj.num <- match(x$data$subject, unique(x$data$subject))
  subjs <- sample(unique(subj.num), replace=TRUE)
  inds <- new.subj <- NULL
  for (i in seq(along=subjs)) {
    subj.inds <- which(subj.num == subjs[i])
    inds <- c(inds, subj.inds)
    new.subj <- c(new.subj, rep(i, length(subj.inds)))
  }
  data.boot <- matrix(nrow=length(inds), ncol=length(x$data$covlabels) + 3)
  subj.name <- deparse(x$call$subject)
  state.name <- deparse(as.list(x$call$formula)[[2]])
  time.name <- deparse(as.list(x$call$formula)[[3]])
  colnames(data.boot) <- c(subj.name, time.name, state.name, x$data$covlabels)
  data.boot[,state.name] <- x$data$state[inds]
  data.boot[,subj.name] <- new.subj
  data.boot[,time.name] <- x$data$time[inds]
  for (j in x$data$covlabels) {
    data.boot[, j] <- x$data$cov[inds, j] + x$data$covdata$covmeans[j]
  }
  as.data.frame(data.boot)
}

### Given a fitted msm model, draw a bootstrap dataset, refit the
### model, and optionally compute a statistic on the refitted model.
### Repeat B times, store the results in a list.
### msm objects tend to be large, so it is advised to compute a statistic on them by specifying "stat", instead
### of using this function to return a list of refitted msm objects.
### To compute more than one statistic, specify, e.g. stat=function(x)list(stat1(x),stat2(x))

### Some of the arguments to the msm call might be user-defined objects.
### e.g. qmatrix, ematrix, hmodel, ...
### Put in help file that these must be in the working environment.

boot.msm <- function(x, stat=pmatrix.msm, B=500, file=NULL){
  boot.list <- vector(B, mode="list")
  for (i in 1:B) {
    boot.data <- if (x$hmodel$hidden || x$cmodel$ncens) bootdata.subject.msm(x) else bootdata.trans.msm(x)
    x$call$data <- substitute(boot.data)
    boot.list[[i]] <- try(eval(x$call))
    if (!is.null(stat))
      boot.list[[i]] <- stat(boot.list[[i]])
    if (!is.null(file)) save(boot.list, file=file)
  }
  boot.list
}

### Utilities for calculating CIs for particular statistics e.g. pmatrix.
### Possibly do also for expected prevalences

pmatrix.ci.msm <- function(x, t, covariates="mean", cl=0.95, B=500) {
  p.list <- boot.msm(x, function(x)pmatrix.msm(x, t, covariates), B)
  p.array <- array(unlist(p.list), dim=c(dim(p.list[[1]]), length(p.list)))
  p.ci <- apply(p.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
  aperm(p.ci, c(2,3,1))
}

totlos.ci.msm <- function(x, start=1, fromt=0, tot=Inf, covariates="mean", cl=0.95, B=500, ...) {
  t.list <- boot.msm(x, function(x)totlos.msm(x, start, fromt, tot, covariates), B)
  t.array <- do.call("rbind", t.list)
  apply(t.array, 2, function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
}

expected.ci.msm <- function(x,
                            times=NULL,
                            timezero=NULL,
                            initstates=NULL,
                            covariates="mean",
                            misccovariates="mean",
                            piecewise.times=NULL,
                            piecewise.covariates=NULL,
                            risk=NULL,
                            cl=0.95, B=500) {
    if(is.null(risk)) risk <- observed.msm(x)$risk
    e.list <- boot.msm(x, function(x){
        expected.msm(x, times, timezero, initstates, covariates, misccovariates, piecewise.times, piecewise.covariates, risk)
    }, B)
    e.tab.array <- array(unlist(lapply(e.list, function(x)x[[1]])), dim=c(dim(e.list[[1]][[1]]), length(e.list)))
    e.perc.array <- array(unlist(lapply(e.list, function(x)x[[2]])), dim=c(dim(e.list[[1]][[2]]), length(e.list)))
    e.tab.ci <- apply(e.tab.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    e.perc.ci <- apply(e.perc.array, c(1,2), function(x)(quantile(x, c(0.5 - cl/2, 0.5 + cl/2))))
    res <- list(aperm(e.tab.ci, c(2,3,1)),  aperm(e.perc.ci, c(2,3,1)))
    names(res) <- c("Expected", "Expected percentages")
    res
}
