\name{heart}
\alias{heart}
\title{Heart transplant monitoring data}
\description{
  A series of approximately yearly angiographic examinations of heart transplant
  recipients.  The state at each time is a grade of coronary artery
  vasculopathy (CAV),  a deterioration of the arterial walls.  
}
\usage{data(heart)}

\format{
  A data frame containing 2846 rows.  There are 622 patients, the rows
  are grouped by patient number and ordered by years after transplant,
  with each row representing an examination and containing additional
  covariates.
  
  \tabular{rll}{
    PTNUM   \tab (numeric) \tab Patient identification number \cr
    age     \tab (numeric) \tab Recipient age at examination (years) \cr
    years   \tab (numeric) \tab Examination time (years after transplant)\cr
    dage    \tab (numeric) \tab Age of heart donor (years) \cr
    sex   \tab (character) \tab sex (0=male, 1=female) \cr
    pdiag \tab (character) \tab Primary diagnosis (reason for
    transplant) IHD=ischaemic heart disease, IDC=idiopathic dilated
    cardiomyopathy. \cr
    cumrej     \tab (numeric) \tab Cumulative number of acute rejection
    episodes \cr
    state \tab (numeric) \tab State at the examination. Grade 1 represents no CAV, grade 2 is mild/moderate CAV and grade 3 is severe
  CAV.  Grade 4 is recorded at the time of death.   \cr
  }}

\keyword{datasets}
