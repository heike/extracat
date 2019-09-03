\name{itab}
\alias{itab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Independence Table
}
\description{
Computes the independence table for a data table with non-negative entries.
The entries of the independence table are defined by the multiplication of all corresponding marginal relative frequencies.
}
\usage{
itab(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data table of any dimension but with non-negative entries.
}
}

\value{
A data table of the same dimension as the input table.}

\author{
Alexander Pilhoefer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
A <- optile(arsim(4000,c(13,17),4,0.1),iter=10)
fluctile(A)
fluctile(itab(A))
D <- A-itab(A)
G <- (D)^2/itab(A)
fluctile(G, tile.col = c(2,4)[1+ (D>0)])

}

