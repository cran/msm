\name{hazard.msm}
\title{Calculate tables of hazard ratios for covariates on transition intensities}
\alias{hazard.msm}
\description{
  Hazard ratios are computed by exponentiating the
  estimated covariate effects on the log-transition intensities.
  This function is called by \code{\link{summary.msm}}.
}
\usage{
hazard.msm(x, hazard.scale = 1)
}
\arguments{
  \item{x}{Output from \code{\link{msm}} representing a fitted
    multi-state model.}

  \item{hazard.scale}{Vector with same elements as number of covariates
    on transition rates. Corresponds to the increase in each covariate
    used to calculate its hazard ratio. Defaults to all 1.}
}
\value{

  A list of tables containing hazard ratio estimates, one table for each covariate.
  Each table has three columns, containing the hazard ratio, and an
  approximate upper 95\% and lower 95\% confidence limit respectively
  (assuming normality of the estimator), for each Markov chain transition intensity. 

}
\seealso{
  \code{\link{msm}}, \code{\link{summary.msm}}, \code{\link{odds.msm}}
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
