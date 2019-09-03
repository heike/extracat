\name{rmbmat}
\alias{rmbmat}
\title{
Pairwise RMB-Plots
}
\description{
This function generates a matrix with RMB-plots of all pairs of variables with a specified target variable. Both categorical and numerical variables are accepted and the latter will be binned. This makes the graphic useful for a mixture of variable types and the binning avoids overplotting and color mash as it occurs in (colored) scatterplots of large datasets.
}
\usage{
rmbmat(x, tv, cut = 20, freqvar = NULL, plot.tv = FALSE, num.mode = TRUE,
 mode = "circular", eqwidth = FALSE, freq.trans = "sqrt", innerval = 1,
  allocation = I,max.scale = 1, use.na = FALSE, expected = FALSE,
  model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl", col.opt = list(),
  label = FALSE, label.opt = list(), diag.opt = list(), lower.opt = list(), 
  upper.opt = list(), rc.opt = list(), factor.opt = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Anything that can be converted to a \code{data.frame} via \code{as.data.frame}.
}
  \item{tv}{
The index of the target variable. The target variable will not be plotted unless \code{plot.tv} is \code{TRUE}.
}
  \item{cut}{
The number of intervals into which numeric variables will be cut.
}
  \item{freqvar}{
An optional frequency variable. \code{"Freq"} is handled automatically.
}
  \item{plot.tv}{
Whether or not to include the target variable(s) in the plot.
}
  \item{num.mode}{
Whether or not to use the numeric mode (no gaps and a numeric axis) for numeric variables.
}
  \item{mode}{
One of \code{"circular"}, \code{"pie"}, \code{"piechart"}, \code{"p"} or \code{"c"} for piecharts, \code{"spine"} or \code{"s"} for spineplots, \code{"bars"}, \code{"bar"} or \code{"b"} for barcharts. 

NOT YET IMPLEMENTED: \code{"rect"} or \code{"r"} for nested rectangles.
 \code{"nested.circles"} are abbreviated by \code{"nc"} or \code{"ncircles"}.
}
  \item{eqwidth}{
See  \link[extracat]{rmb}.
}
  \item{freq.trans}{
See  \link[extracat]{rmb}.
}
  \item{innerval}{
See  \link[extracat]{rmb}.
}
  \item{allocation}{
The widths and heights for the plots are proportional to \code{allocation(nlevels(x))}.
}
  \item{max.scale}{
See  \link[extracat]{rmb}.
}
  \item{use.na}{
See  \link[extracat]{rmb}.
}
  \item{expected}{
See  \link[extracat]{rmb}.
}
  \item{model.opt}{
See  \link[extracat]{rmb}.
}
  \item{gap.prop}{
See  \link[extracat]{rmb}.
}
  \item{gap.mult}{
See  \link[extracat]{rmb}.
}
  \item{col}{
See  \link[extracat]{rmb}.
}
  \item{col.opt}{
See  \link[extracat]{rmb}.
}
  \item{label}{
See  \link[extracat]{rmb}.
}
  \item{label.opt}{
See  \link[extracat]{rmb} and details.
}
  \item{diag.opt}{
A list with \link[extracat]{rmb} parameters. These overwrite the general parameters for all plots on the diagonal.
}
  \item{lower.opt}{
The same as  \code{diag.opt} but for the lower triangular matrix. Additionally it is possible to define a second target variable, e.g.
\code{lower.opt = list(tv2 = 3, ... )}.
}
  \item{upper.opt}{
The same as  \code{diag.opt} but for the upper triangular matrix. Additionally it is possible to define a second target variable, e.g.
\code{upper.opt = list(tv2 = 3, ... )}.
}
  \item{rc.opt}{
A list with which it is possible to define parameters for single matrix cells (plots), columns or rows. This will overwrite all other parameters for the specified plots. It works like this:

\code{ rc.opt = list( r2c12 = list(spine = FALSE), r1 = list(col="rgb"), c4 = list(col="seq") )}
where the plot in row 2 and column 14 is a spineplot, the first row uses RGB colors and the fourth column a sequential color palette.
Later arguments overwrite the preceding ones. For instance in the example the plot in row 1 and column 4 will use the sequential color palette.

}
  \item{factor.opt}{
The same as \code{diag.opt}, \code{lower.opt}, \code{upper.opt} but for all pairs of two categorical variables. This overwrites the other option lists.
}
  \item{\dots}{
Further parameters.
}
}
\details{
Creates a matrix of all pairwise rmb-plots using all possible rmb parameters except \code{cat.ord}, \code{expected =list()} and \code{residuals}.
The parameters are applied to all plots and afterwards possibly overwritten by one of the parameter lists.
}
\value{
An \code{environment} with the parameter lists and matrices. This can be used to update (parts of) the plot without a complete new construction.
The \code{update.rmbmat} function is under development.
}

\author{
Alexander Pilhoefer
}



\seealso{
\link[extracat]{rmb}, \link[graphics]{pairs}
}
\examples{
data(olives)


\dontrun{

# mode = "c" piecharts are currently slow

rmbmat(olives, tv=2, mode = "s")

rmbmat(olives[,1:5], tv=2, col ="div", plot.tv = TRUE,
 lower.opt = list(tv2 = 1, col ="rgb"))

rmbmat(olives[,c(1:5,11)], tv=2,
upper.opt=list(mode="s", eqwidth = TRUE), 
rc.opt = list( c5 = list(eqwidth=FALSE,mode="s"),
r5 = list(eqwidth=TRUE, mode="s")),allocation=NULL)
}

}

