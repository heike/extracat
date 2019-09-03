\name{CBCI}
\alias{CBCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Conditional Independence Bertin Classification Index
}
\description{
Computes the Conditional Independence Bertin Classification Index which uses conditional independence as a reference for normalization. High values indicate that the BCC is not far from the expectation if we know the two marginal 2D BBC values.
}
\usage{
CBCI(x, r = 1,  joint.order = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The 3D table with non-negative entries.
}
  \item{r}{
The index of the conditioning variable, e.g. \code{r = 1} uses the table with variables 2 and 3 conditionally independent given 1 for normalization.
}
 \item{joint.order}{
Whether or not to use a joint ordering for all variables. Otherwise the pairwise values are computed using separate reorderings.
}
}
\details{
The \link[extracat]{BCI} of a 3D table but instead of the total independence case the conditional independence case is used for normalization.}
\value{
Numeric value in [0,1].
}

\author{
Alexander Pilhoefer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{BCI}, \link[extracat]{JBCI}, \link[extracat]{WBCI}
}
\examples{
\dontrun{
A <- optile(arsim(10000, c(11,12,13), 4,0.1))

BCI(A)

CBCI(A,1,TRUE)
CBCI(A,1,FALSE)

}
}

