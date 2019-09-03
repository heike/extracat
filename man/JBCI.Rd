\name{JBCI}
\alias{JBCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Joint Bertin Classification Index
}
\description{
Computes the Joint Bertin Classification Criterion which uses joint independence as a reference for normalization.
}
\usage{
JBCI(x, r = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The 3D table with non-negative entries.
}
  \item{r}{
The index of the variable which is tested for joint independence of the other two.
}
}
\details{
The \link[extracat]{BCI} of a 3D table but instead of the total independence case the joint independence case is used for normalization. With an optimal reordering we have \code{JBCI(x) >= BCI(x)}.
}
\value{
Numeric value in [0,1].
}

\author{
Alexander Pilhoefer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{BCI}, \link[extracat]{CBCI}, \link[extracat]{WBCI}
}
\examples{
\dontrun{
A <- optile( arsim(144*5*20,c(12,12),6,0.1) , iter = 1000)

p1 <- 0.1 + runif(5)
p1 <- p1/sum(p1)

A2 <- apply(A,1:2,function(z) rmultinom(1,z,p1))
A2 <- optile(A2, iter = 1000,return.type="table")

BCI(A)
BCI(A2)

DA2 <- subtable(A2,1:3)
names(DA2) <- c("X","Y","Z","Freq")

rmb(~Y+Z+X,data=DA2)

JBCI(A2,3)
}
}

