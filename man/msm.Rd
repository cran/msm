\name{msm}
\title{Multi-state Markov models}
\alias{msm}
\description{
  Fit a continuous-time Markov multi-state model by maximum
  likelihood. Observations of the process
  can be made at arbitrary times, or the exact times of
  transition between states can be known. 
  Covariates can be fitted to the transition intensities.
  When the true state is observed with error, a model can be fitted
  which simultaneously estimates the state transition intensities and
  misclassification probabilities, with optional covariates on both
  processes.
  
}
\usage{
msm ( formula, qmatrix, misc = FALSE, ematrix, inits, subject,
      covariates = NULL, constraint = NULL, misccovariates = NULL,
      miscconstraint = NULL, qconstraint=NULL, econstraint=NULL,
      covmatch = "previous", initprobs = NULL, 
      data = list(), fromto = FALSE, fromstate, tostate, timelag,
      death = FALSE, tunit = 1.0, exacttimes = FALSE,
      fixedpars = NULL, ... )
}
\arguments{
  \item{formula}{ A formula giving the vectors containing the observed
    states and the  corresponding observation times. For example,
    
    \code{states ~ times}

    See \code{fromto} for an alternative way to specify the
    data. Observed states should be in the set \code{1, \dots, n},
    where \code{n} is the number of states.
  }

  \item{qmatrix}{Matrix of indicators for the allowed transitions.
    If a transition is allowed from state \eqn{r} to state \eqn{s},
    then \code{qmatrix} should have \eqn{(r,s)} entry 1, otherwise
    it should have \eqn{(r,s)} entry 0. The diagonal of \code{qmatrix}
    is ignored. For example,\cr
    
    \code{
    rbind(
    c( 0, 1, 1 ),
    c( 1, 0, 1 ),
    c( 0, 0, 0 )
    )
  }\cr
  
    represents a 'health - disease - death' model, with transitions
    allowed from health to disease, health to death, disease to health, 
    and disease to death. 
  }

  \item{misc}{Set \code{misc = TRUE} if misclassification between
    observed and underlying states is to be modelled.}

  \item{ematrix}{
    (required when \code{misc == TRUE}) Matrix of indicators for the allowed
misclassifications. 
    The rows represent underlying states, and the columns represent
    observed states.
    If an observation of state \eqn{s} is possible when the subject
    occupies underlying state \eqn{r}, then \code{ematrix} should have
    \eqn{(r,s)} entry 1, otherwise
    it should have \eqn{(r,s)} entry 0. The diagonal of \code{ematrix}
    is ignored. For example, \cr

    \code{
      rbind(
    c( 0, 1, 0 ),
    c( 1, 0, 1 ),
    c( 0, 1, 0 )
    )
  }\cr
  
    represents a model in which misclassifications are only permitted
    between adjacent states. 
  }

  \item{inits}{(required) Vector of initial parameter estimates for the
    optimisation. These are given in the order

    - transition intensities (reading across first rows of intensity
    matrix, then second row ... )

    - covariate effects on log transition intensities

    - misclassification probabilities (reading across first row of
    misclassification matrix, then second row ...)

    - covariate effects on logit misclassification probabilities

    Covariate effects are given in the following order,

    - effects of first covariate on transition/misclassification matrix elements (reading across
    first row of matrix, then second row ...)

    - effects of second covariate ... 
  }

  \item{subject}{Vector of subject identification numbers, when the data
    are specified by \code{formula}. If missing, then all observations
    are assumed to be on the same subject. These must be sorted so that
    all observations on the same subject are adjacent. Ignored if
    \code{fromto == TRUE}.}

  \item{covariates}{Formula representing the covariates on the
  transition intensities, for example,

    \code{~ age + sex + treatment}
    
  }

  \item{constraint}{A list of one vector for each named covariate. The
  vector indicates which covariate effects on intensities are
  constrained to be equal. Take, for example, a model with five
  transition intensities and two covariates. Specifying\cr
    
    \code{constraint = list (age = c(1,1,1,2,2),  treatment = c(1,2,3,4,5))}\cr

    constrains the effect of age to be equal for the first three
    intensities, and equal for the fourth and fifth. The effect of
    treatment is assumed to be different for each intensity. Any vector of
    increasing numbers can be used as indicators. The intensity parameters are
    assumed to be ordered by reading across the rows of the
    transition matrix, starting at the first row.

    For categorical covariates, defined using \code{factor(covname)},
    specify constraints as follows:\cr

    \code{list(..., covnameVALUE1 = c(...), covnameVALUE2 = c(...), ...)}\cr
    
    where \code{VALUE1}, \code{VALUE2}, ... are the levels of the factor.
    Make sure the \code{contrasts} option is set appropriately, for
    example, the default \code{options(contrasts=c(contr.treatment,
      contr.poly))} sets the first (baseline) level of unordered factors to
    zero.

    To assume no covariate effect on a certain transition, set its
    initial value to zero and use the \code{fixedpars} argument to fix
    it during the optimisation.
    
  }

  \item{misccovariates}{A formula representing the covariates on the 
  misclassification probabilities, analogously to \code{covariates}.
  }

  \item{miscconstraint}{A list of one vector for each named covariate on
    misclassification probabilities. The vector indicates which
    covariate effects on misclassification probabilities are
    constrained to be equal, analogously to \code{constraint}.
  }

  \item{qconstraint}{A vector of indicators specifying which baseline
    transition intensities are equal. For example,
    
    \code{qconstraint = c(1,2,3,3)}

    constrains the third and fourth intensities to be equal, in a model
    with four allowed instantaneous transitions. 
  }

  \item{econstraint}{A similar vector of indicators specifying which
    baseline misclassification probabiliities are constrained to be
    equal. }
  
  \item{covmatch}{If \code{"previous"}, then time-dependent covariate
    values are taken from the observation at the start of the
    transition. If \code{"next"}, then the covariate value is taken from
    the end of the transition.}

  \item{initprobs}{Currently only used in misclassification
    models.  Vector of assumed underlying state occupancy
    probabilities at each individual's first observation. Defaults to
    \code{c(1, rep(0, nstates-1))}, that is, in state 1 with a
    probability of 1.  
  }

  \item{data}{Optional data frame in which to interpret the state, time,
    subject ID, covariate, fromstate, tostate and timelag vectors.}

  \item{fromto}{If \code{TRUE}, then the data are given as three vectors, 
    \emph{from-state, to-state, time-difference},
    representing the set of observed transitions between states, and the
    time taken by each one. Otherwise, the data are given by \code{formula},
    containing observation times, corresponding observed states, and a
    vector of corresponding subject identification numbers \code{subject}.
  }

  \item{fromstate}{ Starting states for the observed transitions
    (required if \code{fromto == TRUE} ). }

  \item{tostate}{ Finishing states for the observed transitions
    (required if \code{fromto == TRUE} ). }

  \item{timelag}{ Time difference between observing \code{fromstate}
    and \code{tostate}  (required if \code{fromto == TRUE} ). }

  \item{death}{Vector of indices of the death states.  A death state is
    an absorbing state whose time of entry is known exactly, but the
    individual is assumed to be in an unknown transient state ("alive")
    at the previous instant.  This is the usual situation for times of death in
    chronic disease monitoring data.  For example, if you specify
    \code{death = c(4, 5)} then states 4 and 5 are assumed to be death
    states.

    \code{death = TRUE} indicates that the final state is a death state,
    and \code{death = FALSE} (the default) indicates that there is no
    death state.
  }

  \item{tunit}{No longer used, since \pkg{msm} version 0.3.2. The entry time
    into death states is now assumed to be exact, rather than known
    within a day.}

  \item{exacttimes}{By default, the transitions of the Markov process
  are assumed to take place at unknown occasions in between the
  observation times.  If \code{exacttimes} is
    set to \code{TRUE}, then the observation times are assumed to
    represent the exact and complete times of transition of the process. 
  }

  \item{fixedpars}{Vector of indices of parameters whose values will be
    fixed at their initial values during the optimisation. These
    correspond to indices of the \code{inits} vector, whose order is
    specified above. This can be useful for building complex models
    stage by stage. }

  \item{...}{Optional arguments to the general-purpose R
    optimization routine \code{\link{optim}}. Useful options include
    \code{method="BFGS"} for using a quasi-Newton optimisation
    algorithm, which can often be faster than the default Nelder-Mead.
    If the optimisation fails to converge, consider normalising the
    problem using, for example, \code{control=list(fnscale = 2500)}, for
    example, replacing 2500 by a number of the order of magnitude of the
    likelihood. If 'false' convergence is reported and the standard
    errors cannot be calculated due to a non-positive-definite Hessian,
    then consider tightening the tolerance criteria for convergence. If
    the optimisation takes a long time, intermediate steps can be
    printed using the \code{trace} argument of the control list. See
    \code{\link{optim}} for details.
  }
}
\value{
  A list of class \code{msm}, with components:

  \item{misc}{Logical value indicating whether the model includes
    misclassification.}
  \item{Qmatrices}{A list of matrices. The first component, labelled
    \code{logbaseline}, is a matrix containing the estimated
    transition intensities on the log scale with any covariates fixed at
    their means in the data. Each remaining component is a matrix giving the linear
    effects of the labelled covariate on the matrix of log
    intensities.  To extract an estimated intensity matrix on the natural
    scale, at an arbitrary combination of covariate values, use the
    function \code{\link{qmatrix.msm}}.
  }
  \item{QmatricesSE}{The standard error matrices corresponding to
    \code{Qmatrices}.
  }
  \item{Ematrices}{A list of matrices. The first component, labelled
    \code{logitbaseline}, is the estimated misclassification probability
    matrix  with any covariates fixed at their means
    in the data. Each remaining component is a matrix giving the linear
    effects of the labelled covariate on the matrix of logit
    misclassification probabilities.   To extract an estimated misclassification
    probability matrix on the natural scale, at an arbitrary combination
    of covariate values, use the function \code{\link{ematrix.msm}}.}
  \item{EmatricesSE}{The standard error matrices corresponding to \code{Ematrices}.}
  \item{sojourn}{
    A list with components:

    \code{mean} = estimated mean sojourn times in the transient states,
    with covariates fixed at their means.

    \code{se} = corresponding standard errors.
  }
  \item{minus2loglik}{Minus twice the maximised log-likelihood.}
  \item{foundse}{Logical value indicating whether the Hessian was positive-definite
    at the supposed maximum of the likelihood.  If not, the covariance matrix of the
    parameters is unavailable. Furthermore, it is doubtful whether the optimisation
    has converged to the maximum. 
  }
  \item{estimates}{Vector of untransformed maximum likelihood estimates
    returned from \code{\link{optim}}, with intensities on the log scale
    and misclassification probabilities on the logit scale.}
  \item{estimates.t}{Vector of transformed maximum likelihood estimates
    with intensities and probabilities on their natural scales.}
  \item{covmat}{Covariance matrix corresponding to \code{estimates}.}
  \item{data}{A list of constants and vectors giving the data, for use
    in post-processing.}
  \item{model}{A list of constants and vectors specifying the model, 
    for use in post-processing.}
}
\details{
  For full details about the methodology behind the \pkg{msm} package,
  refer to the PDF manual \file{msm-manual.pdf} in the \file{doc}
  subdirectory of the package. This includes a tutorial in the typical
  use of \pkg{msm}.

  For models without misclassification,
  the likelihood is calculated in terms of the transition intensity
  matrix \eqn{Q}. When the data consist of observations of the Markov
  process at arbitrary times, the exact transition times are not known.
  Then the likelihood is calculated using the transition probability
  matrix \eqn{P(t) = exp(tQ)}. If state \eqn{i} is observed at time
  \eqn{t} and  state \eqn{j} is observed at time \eqn{u}, then the
  contribution to the likelihood from this pair of observations is
  the \eqn{i,j} element of \eqn{P(u - t)}. See, for example, Kay (1986),
  or Gentleman \emph{et al.} (1994).

  For models with misclassification, the likelihood for an individual
  with \eqn{k} observations is calculated by summing over the unknown
  state at each time,  producing a product of \eqn{k} matrices. The
  calculation, adapted from Satten and Longini (1996), is given by
  Jackson and Sharples (2002).

  There must be enough information in the data on each state to estimate
  each transition rate, otherwise the likelihood will be flat and the
  maximum will not be found.  It may be appropriate to reduce the
  number of states in the model, or reduce the number of covariate
  effects, to ensure convergence.

  Choosing an appropriate set of initial values for the optimisation can
  also be important.  For flat likelihoods, 'informative' initial values
  will often be required.
}
\references{
  Jackson, C.H., Sharples, L.D., Thompson, S.G. and Duffy, S.W. and
  Couto, E.  Multi-state Markov models for disease progression with
  classification error. \emph{The Statistician}, 52(2): 193--209 (2003)
    
  Jackson, C.H. and Sharples, L.D. Hidden Markov models for the
  onset and progresison of bronchiolitis obliterans syndrome in lung
  transplant recipients \emph{Statistics in Medicine}, 21(1): 113--128
  (2002).
    
  Kay, R.  A Markov model for analysing cancer markers and disease
  states in survival studies.  \emph{Biometrics} (1986) 42: 855--865.
  
  Gentleman, R.C., Lawless, J.F., Lindsey, J.C. and Yan, P.  Multi-state
  Markov models for analysing incomplete disease history data with
  illustrations for HIV disease.  \emph{Statistics in Medicine} (1994) 13(3):
  805--821.

  Satten, G.A. and Longini, I.M.  Markov chains with measurement error:
  estimating the 'true' course of a marker of the progression of human
  immunodeficiency virus disease (with discussion) \emph{Applied
    Statistics} 45(3): 275-309 (1996)
}

\seealso{
  \code{\link{simmulti.msm}}, \code{\link{print.msm}}, \code{\link{plot.msm}},
  \code{\link{summary.msm}}, \code{\link{qmatrix.msm}},
  \code{\link{pmatrix.msm}},  \code{\link{sojourn.msm}}, 
}
\examples{
data(aneur)
print(aneur$from[1:10])
print(aneur$to[1:10])
print(aneur$dt[1:10])
### four states corresponding to increasing disease severity,
### with progressive transitions only 
qmat <- rbind( c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1), c(0, 0, 0,
0))

aneurysm.msm <- msm(data=aneur, fromto=TRUE, fromstate=from, tostate=to,
                    qmatrix=qmat, timelag=dt, inits=c(0.001, 0.03, 0.3),
                    method="BFGS", control=list(trace=2))

print(aneurysm.msm)
qmatrix.msm(aneurysm.msm) # Extract only the transition intensities from the printed results
pmatrix.msm(aneurysm.msm, 10) # Extract the 10 year transition probability matrix
sojourn.msm(aneurysm.msm) # Extract the mean sojourn times
}
\author{C. H. Jackson \email{chris.jackson@imperial.ac.uk}}
\keyword{models}
