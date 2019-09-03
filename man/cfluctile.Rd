\name{cfluctile}
\alias{cfluctile}
\title{
Pseudo-Diagonal Partitioning for two-way tables
 }
\description{
Identifies a diagonal of block-clusters in a two-way table using a top-down-partitioning algorithm then plots the table and adds the clusters as rectangles.
}
\usage{
cfluctile(x, tau0 = NULL, method = "Kendall", nsplit = NULL,
                 maxsplit = NULL, trafo = I, gap.prop = 0.2, floor = 0,
                 rev.y = FALSE, add = FALSE, shape = "r", just = "c",
                 dir = "b", plot = TRUE, rect.opt = list(), border =
                 NULL, label = TRUE, lab.opt = list(), tile.col =
                 hsv(0.1, 0.1, 0.1, alpha = 0.6),
                 tile.border = NA, bg.col = "lightgrey", \dots)
}
\arguments{
  \item{x}{
A 2-way table or matrix.
}
  \item{tau0}{
The minimum acceptable value of Kendall's tau, Cohen's Kappa or WBCI. Defaults to the criterion of the input matrix \code{x}.
}
\item{method}{
Either \code{"Kendall"} for Kendall's Tau, \code{"Cohen"} for Cohen's Kappa,  \code{"WBCI"} for the Weighted Bertin Classification Criterion and \code{"s"} for the minimum residual method.
}
\item{nsplit}{
The number of splits to make. \code{tau0} is ignored.
}
\item{maxsplit}{
The maximum number of splits.
}
 \item{trafo}{
A transformation of the table entries for the plot, but not for the computation of the splits.
E.g. \code{trafo = function(z) log(1+z)}.
}
  \item{gap.prop}{
proportion of the gaps between rows/columns.
}
  \item{floor}{
floor censored zooming: all cases will be plotted but only those with a frequency of at least \code{floor} will be considered for the clustering.
}
  \item{rev.y}{
revert the y axis.
}
\item{add}{
	Whether to make a new plot or to add to an existing one.
}
\item{shape}{
	The shape of the objects. See \link[extracat]{fluctile}.
}
\item{just}{
	See \link[extracat]{fluctile}.
}
\item{dir}{
	See \link[extracat]{fluctile}.
}
\item{plot}{
	Whether or not to create a plot via \code{fluctile}.
}

 \item{ rect.opt }{A list with optional parameters for the rectangles. Possible parameters are:
 	\tabular{rlll}{
         \code{col}   \tab The rectangle color.\cr
      	 \code{lwd} \tab The line width. Default is \code{"red"}.	\cr
         \code{lty} \tab The line type. Default is 1 (solid).	\cr
         \code{fill} \tab The color to fill the rectangles. Defaults is \code{NULL}.
          A sensible choice is for instance \code{alpha(col, 0.1)}.	\cr
        }
}
\item{border}{
	The white margins around the plot which are also used for the labels.
	Must be a vector of length 1, 2 or 4 with values in [0, 1]. Default is \code{border = 0.05}.
}
  \item{label}{
Whether or not to draw labels.
}
  \item{lab.opt}{
Label options, see \link[extracat]{fluctile}.
}
  \item{tile.col}{
Color(s) for the tiles, see \link[extracat]{fluctile}.
}
  \item{tile.border}{
Border color for the tiles. Can also be a matrix.
}
  \item{bg.col}{
Color for the background of the cells, see \link[extracat]{fluctile}.
}


  \item{\dots}{
dots}
}
\details{
This function calls \link[extracat]{fluctile} to create a 2-way fluctuation diagram and then adds cluster rectangles to it. The cluster rectangles are computed in the following way: 

The algorithm cuts the data matrix once horizontally and once vertically and computes a criterion for the 2x2 table consisting of the sums of the four parts that resulted from the cuts. This is done for all possible horizontal and vertical cuts and the best combination is chosen.
Then the same procedure is applied to the bottom right submatrix and the top left submatrix. The algorithms stops if no cut yields a better criterion value than \code{tau0}.

}
\value{
\code{invisible(TRUE)}
}

\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}
\note{
This was part of the Google Summer of Code 2011.
}


\seealso{
\link[extracat]{optile}, \link[extracat]{sortandcut}, \link[extracat]{tfluctile}
}

\examples{

M <- arsim(10000,c(30,40),8, noise = 0.4)
cfluctile( M2 <- optile(M,iter=20) )

cfluctile( M3 <- sortandcut(M) )

cfluctile( M3, nsplit = 4 )

cfluctile( M3, maxsplit = 12 )

cfluctile( M3, tau0 = 0.8 )


}
