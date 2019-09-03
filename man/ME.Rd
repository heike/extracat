\name{ME}
\alias{ME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measure of Effectiveness
}
\description{
Computes the measure of effectiveness for a table, a matrix or an array.
}
\usage{
ME(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix, table or array.
}
}

\value{
The ME value.
}

\seealso{
\link[extracat]{optME}
}
\examples{
a <- arsim(2000,c(8,9,10),3,0.2)
ME(a)
a2<-optME(a)
ME(a2)
}
