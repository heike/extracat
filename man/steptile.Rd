\name{steptile}
\alias{steptile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
stepwise reordering
}
\description{
Starts with the first \code{k+1} variables and applies \link[extracat]{optile} to the corresponding subtable. Then one additional variable at a time is reordered using the subtable defined by this variable and the last \code{k} variables. Only the current variable is reordered and the others are fixed since they have been reordered in the previous steps.
}
\usage{
steptile(x, k = 1, cpcp = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The \code{data.frame} (which is better for high-dimensional data) or data table.
}
  \item{k}{
The number of preceding variables used for the reordering. E.g. if \code{k = 3} then variable \code{6} is reordered using the variables \code{3, 4, 5, 6}.
}
\item{cpcp}{
If \code{TRUE} a special version of the algorithm which minimizes crossings in CPCP plots (e.g. \link[extracat]{scpcp} ) is used. This modification works with aggregations of the last \code{k} variables and is much faster than the standard procedure if \code{k > 1}.
}
\item{\dots}{
Arguments passed to \link[extracat]{optile}.
}
}
\details{
The optile function also offers stepwise reordering via the argument \code{method = "sw"} but always starts with the first pair of variables and then considers the complete past: for the reordering of variable \code{i} all variables \code{1...(i-1)} are considered. The stepwise algorithms are applicable to high-dimensional problems with a large number of variables where the multivariate techniques fail.
Even if \code{k} is high (i.e. the subtables are also high-dimensional) the procedure is very fast since it ca use the following trick: instead of applying optile to the multidimensional table it is applied to a 2D-table with one dimension defined by the variable that is reordered and the other dimension defined by the
(ordered) combinations of all other variables. This way only combinations which appear at least once in the dataset matter and all empty entries (the majority in high-dimensional tables) can be left aside. The maximum possible size of such a table is therefore \eqn{N * max(n_i)} when \eqn{N} is the number of observations and \eqn{n_i} is the number of categories in dimension \code{i}.
}
\value{
The reordered data either as a \code{table} or \code{data.frame} depending on the input type.
}

\author{
Alexander Pilhoefer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[extracat]{optile}
}
\examples{
\dontrun{

# scaled numeric variables from the olives data
# and 20 k-means solutions
so <- scale(olives[,3:10])
rr <- replicate(20,{
	kmeans(so,8)$cluster
})

# par(mfrow=c(3,1))

# initial cluster orders
x <- as.data.frame(cbind(olives[,1:2],rr))
require(scales)
scpcp(x, sel = "data[,1]",
		 sel.palette="rgb", col.opt=list(alpha=0.5))

# reordering using steptile. 
# optile does not work for the complete table since it has 9*3*2^60 > 3E19 entries
# colors by the first unordered example:

x2 <- steptile(x, k = 4)
scpcp(x2, sel = "match(data[,1],levels(.GlobalEnv$x[,1]))", 
		sel.palette="rgb", col.opt=list(alpha=0.5))


# additionally reordering the variables ... cmat takes about 20-30 seconds
CM <- cmat(x[,3:22])
require(seriation)
sM <- get_order(seriate(1-CM))

x3 <- steptile(x2[, c(1,2,2+sM,23)], k = 4)

scpcp(x3, sel = "match(data[,1],levels(.GlobalEnv$x[,1]))",
		 sel.palette="rgb", col.opt=list(alpha=0.5))

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
