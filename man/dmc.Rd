\name{dmc}
\alias{dmc}
\docType{data}
\title{
dmc 2009 insurance variables
}
\description{
five insurance variables from the dmc 2009 dataset, which have a ordinal structure which has been lost somehow. Can we find it again?
}
\usage{data(dmc)}
\format{
  A data frame with 693 observations on the following 6 variables.
  \describe{
    \item{\code{eiw_scr}}{a factor with levels \code{5} \code{6} \code{4} \code{3} \code{2} \code{1}}
    \item{\code{eih_scr}}{a factor with levels \code{6} \code{3} \code{5} \code{1} \code{4} \code{2}}
    \item{\code{ifi_scr}}{a factor with levels \code{4} \code{3} \code{5} \code{2} \code{1} \code{6}}
    \item{\code{tec_scr}}{a factor with levels \code{5} \code{1} \code{3} \code{2} \code{4} \code{6}}
    \item{\code{klv_scr}}{a factor with levels \code{2} \code{5} \code{6} \code{1} \code{3} \code{4}}
    \item{\code{Freq}}{a numeric vector}
  }
}
\details{
The Data Mining Cup (dmc) is a competition for students.
}

\examples{
data(dmc)
}
\keyword{datasets}
