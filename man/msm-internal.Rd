\name{msm-internal}
\title{Internal msm functions}
\alias{msm-internal}
\alias{lik.msm}
\alias{lik.msmmisc}
\alias{getobs.msm}
\alias{msm.process.covs}
\alias{msm.form.output}
\alias{msm.check.consistency}
\alias{msm.check.times}
\alias{msm.check.tunit}
\description{Auxiliary R functions used in the msm package. These are not
  meant to be called by the user}
\usage{

lik.msm(params, allinits, misc, subject, time, state, tostate, fromto,
        qvector, evector, covvec, constrvec, misccovvec, miscconstrvec, initprobs,
        nstates, nintens, nmiscs, nobs, npts, ncovs, ncoveffs,
        nmisccovs, nmisccoveffs, covmatch, death, tunit, exacttimes, fixedpars, plabs)

getobs.msm(sim, obstimes, death, tunit)

msm.process.covs(covariates, data, constraint, nobs, nmatrix)

msm.form.output(matrix, nstates, nmatrix, npars, ncovs, constrvec, fixedpars, fixed, covlabels,
                params, covmat, ses, foundstderr, params.offset)

msm.check.consistency(qmatrix, misc, fromto=F, subject=NULL, state,
                      tostate=NULL, time, death=FALSE, exacttimes=FALSE, tunit=NULL)

msm.check.times(time, subject)

msm.check.tunit(fromto, time, subject, timelag, tunit)

}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{internal}
