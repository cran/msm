\name{hazard.msm}
\title{Calculate tables of hazard ratios for covariates on transition intensities}
\alias{hazard.msm}
\description{
  Hazard ratios are computed by exponentiating the
  estimated covariate effects on the log-transition intensities.
  This function is called by \code{\link{summary.msm}}.
}
\usage{
hazard.msm(msm, hazard.scale = 1)
}
\arguments{
  \item{msm}{Output from \code{\link{msm}} representing a fitted
    multi-state model.}

  \item{hazard.scale}{Vector with same elements as number of covariates
    on transition rates. Corresponds to the increase in each covariate
    used to calculate its hazard ratio. Defaults to all 1.}
}
\value{

  A list of tables containing hazard ratio estimates, one table for each covariate.
  Each table has three columns, containing the hazard ratio, upper 95\% and lower 95\%
  confidence limit respectively, for each Markov chain transition
  intensity. 

}
\seealso{
  \code{\link{msm}}, \code{\link{summary.msm}}
}
\author{C. H. Jackson \email{chris.jackson@ic.ac.uk}}
\keyword{models}