\name{msm-internal}
\title{Internal msm functions}
\alias{msm-internal}
\alias{lik.msm}
\alias{lik.msmmisc}
\alias{getobs.msm}
\alias{observed.msm}
\alias{msm.process.covs}
\alias{msm.form.output}
\alias{msm.check.qmatrix}
\alias{msm.check.ematrix}
\alias{msm.check.consistency}
\alias{msm.check.state}
\alias{msm.check.model}
\alias{msm.check.times}
\alias{msm.check.tunit}
\alias{msm.check.qmatrix}
\alias{qematrix.msm}
\alias{qematrix.diagse.formstr}
\alias{qematrix.diagse.msm}
\alias{expsum}
\alias{expitsum}
\alias{qratio.se.msm}
\alias{transient.msm}
\alias{absorbing.msm}
\description{Auxiliary R functions used in the msm package. These are not
  meant to be called by the user}
\usage{

lik.msm(params, allinits, misc, subject, time, state, tostate, fromto,
        qvector, evector, covvec, constrvec, misccovvec, miscconstrvec,
        baseconstrvec, basemiscconstrvec, initprobs,
        nstates, nintens, nintenseffs, nmisc, nmisceffs, nobs, npts, ncovs, ncoveffs,
        nmisccovs, nmisccoveffs, covmatch, death, tunit, exacttimes, fixedpars, plabs)

getobs.msm(sim, obstimes, death, tunit)

observed.msm(x, times)

msm.process.covs(covariates, data, constraint, nobs, nmatrix)

msm.form.output(matrix, nstates, nmatrix, nmatrixeffs, npars, ncovs,
                constrvec, baseconstrvec, fixedpars, fixed, covlabels,
                params, covmat, foundse, params.offset)

msm.check.qmatrix(qmatrix)

msm.check.ematrix(ematrix, qmatrix)

msm.check.consistency(qmatrix, misc, fromto=F, subject=NULL, state,
                      tostate=NULL, time, death=FALSE, exacttimes=FALSE, tunit=NULL)

msm.check.state(nstates, state, fromto=FALSE, tostate=NULL)

msm.check.model(state, subject, qmatrix, fromto=NULL, tostate=NULL)

msm.check.times(time, subject)

msm.check.tunit(fromto, time, subject, timelag, tunit)

qematrix.msm(x, covariates="mean", which="intens", sojourn)

qematrix.diagse.msm(x, covariates="mean", which="intens", sojourn, ni,
nieffs, ivector, nc, constrvec, baseconstrvec, covlabels, covmeans, trsum)

qematrix.diagse.formstr(nir, inds, cur.i, ni, nc, coefs, offset, trsum)

expsum(inds, coefs)

expitsum(inds, coefs)

qratio.se.msm(x, ind1, ind2, covariates="mean")

transient.msm(x)

absorbing.msm(x)

}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{internal}
