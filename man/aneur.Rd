\name{aneur}
\alias{aneur}
\title{Aortic aneurysm progression data}
\description{
  This dataset contains a list of transitions between
  grades of aortic aneurysms, measured by ultrasound examination of the
  diameter of the aorta. 
  
}
\usage{data(aneur)}

\format{
  A data frame containing 2662 rows. Columns correspond to
  'starting state', 'end state', 'time interval', 'age at beginning of
  interval'.  The stages represent successive degrees of aneurysm
  severity. 

  \tabular{rll}{
    Stage 1 \tab Aneurysm-free \tab < 30 cm \cr
    Stage 2 \tab Mild aneurysm \tab 30-44 cm \cr
    Stage 3 \tab Moderate aneurysm \tab 45-54 cm \cr
    Stage 4 \tab Severe aneurysm \tab > 55 cm \cr
  }
  
  These were obtained from 4337 ultrasound scans from 838 men over 65.
  683 of these men were aneurysm-free at age 65 and were re-screened every two
  years.  The remaining men were aneurysmal at entry and had successive
  screens with frequency depending on the stage of the aneurysm.
  Stage 4 aneurysms are repaired by surgery.
}

\source{
  The Chichester, U.K. randomised controlled trial of screening for
  abdominal aortic aneurysms by ultrasonography.
}

\references{
  Couto, E. and Duffy, S. W. and Ashton, H. A. and Walker, N. M.
  and Myles, J. P. and Scott, R. A. P. and Thompson, S. G. (2002)
  \emph{Probabilities of progression of aortic aneurysms: estimates
    and implications for screening policy}
  Journal of Medical Screening 9(1):40--42
}

\keyword{datasets}
