\name{msm.summary}
\title{Summarise a fitted multi-state model}
\alias{summary.msm}
\alias{print.summary.msm}
\description{
  Summary method for fitted \code{\link{msm}} models.
  For models without misclassification, produces a table of observed and
  expected state prevalences for each time. For all models, prints
  hazard ratios with confidence intervals for covariate effects.
}
\usage{
\method{summary}{msm}(object, prevtimes=NULL, timezero=NULL, initstates=NULL,
hazard.scale=NULL,\dots)
\method{print}{summary.msm}(x,\dots)
}
\arguments{
  \item{object}{A fitted multi-state model object, as returned by
    \code{\link{msm}}.}

  \item{prevtimes}{A sequence of times at which to compare observed
    and expected prevalences of each state. Defaults to
    \code{seq(min(times), max(times), (max(times) - min(times))/10).
    }
  }
  \item{timezero}{Initial time of the Markov process. Expected values
    are forecasted from here. Defaults to the minimum of the observation
    times given in the data. }
  \item{initstates}{ Optional vector of the same length as the number of
    states. Gives the numbers of individuals occupying each state at
    \code{timezero}. The default is that all individuals are in state 1.  }
  \item{hazard.scale}{Vector with same elements as number of covariates
    on transition rates. Corresponds to the increase in each covariate
    used to calculate its hazard ratio. Defaults to all 1.}  
  \item{x}{An object of class \code{summary.msm} usually, a result of a
    call to \code{\link{summary.msm}}}
  \item{...}{further arguments passed to or from other methods}
}

\note {
  The prevalence table is not available if the data were specified in
  \code{fromto} format in the call to {\code{\link{msm}}}. Expected
  prevalences only make sense for processes where all individuals
  start at a common time.
}

\value{
  A list of class \code{summary.msm}, with components:
  
  \item{prevalences}{Output from \code{\link{prevalence.msm}} or \code{\link{prevalencemisc.msm}}}

  \item{hazard}{Output from \code{\link{hazard.msm}}}.

  \item{hazard.scale}{Value of the \code{hazard.scale} argument}
    
}
\seealso{
  \code{\link{msm}}, \code{\link{print.summary.msm}}
  \code{\link{prevalence.msm}}, \  \code{\link{prevalencemisc.msm}}, \code{\link{hazard.msm}}
}

\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{}