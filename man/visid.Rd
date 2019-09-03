\name{visid}
\alias{visid}
\alias{visna}
\alias{visdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualizing Indicator Matrices and Missing Values
}
\description{
This function aggregates a binary dataframe or matrix using \link[extracat]{subtable} and visualizes the combinations along with the marginal distributions using \link[extracat]{fluctile}. Options include reordering of rows and columns, filtering the most important rows and columns aas well as ceiling censored zooming for the marginals.\cr
\cr Missing values can be visualized using \link{is.na} and datasets with categorical variables can be transoformed via \link[extracat]{idat}.
\code{visna(x)} is a shortcut for \code{visid( is.na(x) + 0 )}.
\code{visdf(x, freqvar)} is a shortcut for \code{visid( idat(x, allcat = TRUE, keep = freqvar) )}.
}
\usage{
visid(x, freqvar = "Freq", tp = FALSE, fr = 1, fc = 1, sort = "n",
 sort.method = "count", col = "w",
  mar.col = c(alpha("black", 0.7), alpha("darkred", 0.8), "red", "green"),
   s = Inf, pmax = 1, opts = list(), plot = TRUE, return.data = !plot, \dots)
visna(x, freqvar = "Freq", tp = FALSE, fr = 1, fc = 1, sort = "n",
 sort.method = "count", col = "w",
  mar.col = c(alpha("black", 0.7), alpha("darkred", 0.8), "red", "green"),
   s = Inf, pmax = 1, opts = list(), plot = TRUE, return.data = !plot, \dots)
visdf(x, freqvar = "Freq", tp = FALSE, fr = 1, fc = 1, sort = "n",
 sort.method = "count", col = "w",
  mar.col = c(alpha("black", 0.7), alpha("darkred", 0.8), "red", "green"),
   s = Inf, pmax = 1, opts = list(), plot = TRUE, return.data = !plot, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A binary dataframe or matrix. \link{is.na} and \link[extracat]{idat} can be used to bring in missing values and categorical variables.
}
  \item{freqvar}{
An optional frequency variable. If this is not found the data is aggregated using\link[extracat]{subtable}.
}
  \item{tp}{
Logical. Whether or not to transpose the indicator matrix for the visualization.
}
  \item{fr}{
This controls the row filtering: Only the \code{fr} most frequent rows are kept. Values below 1 are interpreted as proportions and only the minmum number of rows covering at least \code{fr} percent of the observations are kept.

}
  \item{fc}{
See \code{fr}.
}
\item{sort}{
One of \code{"n"} (no reordering), \code{"r"} (row reordering), \code{"c"} (column reordering) or \code{"b"} (row and column reordering).
}
  \item{sort.method}{
The default is reordering by frequency (rows) and average (columns). Other options include \code{"ME"} for reordering by the measure of effectiveness (see \link[extracat]{optME}) and
\code{"optile"} which enables all reorderings offered by \link[extracat]{optile}.
}
  \item{col}{
The color palette. For the basic indicator case only a single color is meaningful. Transfrmations via \link[extracat]{idat} feature coloring by variable.
}
  \item{mar.col}{
Colors for the marginals: rows, columns, row markers, column markers (see \code{s} and \code{pmax}).
}
  \item{s}{
Ceiling censored zooming for the rows. The maximum of the scale is \code{s} times the second largest value.
}
  \item{pmax}{
The maximum for the average/percentage scale used for the column marginal plot.
}
  \item{opts}{
A list of options. Currently:
\tabular{rlll}{
         \code{gap.prop} \tab gaps proportion used in \link[extracat]{fluctile}). \cr 
         \tab The default is to use no gaps if more than \cr
         \tab 40 rows or columns are involved.\cr
         
      	 \code{mar} \tab Vector controlling the size of the marginal plots. \cr
	 \code{border} \tab Borders around the plot. \cr
	  \code{shape} \tab Shape of the symbols. Default is rectangles. \cr   
	    \code{bg.col} \tab Background color center/right/bottom \cr    
	      \code{abbrev} \tab Label abbreviation. \cr         
	  }
}
  \item{plot}{
Whether or not to draw the plot.
}
  \item{return.data}{
Whether or not to return the data after filtering and reordering took place.
}
  \item{\dots}{
Further arguments passed to either \link[extracat]{optME} or \link[extracat]{optile}. Can for instance be used to choose the reordering method.
}
}

\value{
The filtered and reordered data or \code{invisible(TRUE)} (default).
}

\author{
Alexander Pilhoefer
}

\seealso{
\link[extracat]{fluctile}, \link[extracat]{optile}, \link[extracat]{idat}
}
\examples{

\dontrun{
require(reshape2)
require(scales)

MJ <- read.table(
"http://www.rosuda.org/~pilhoefer/MJnew.txt",
header=T,sep="\t",quote="")


MJS <- MJ[,13:105]

visid(MJS)

# sort by count/percentage
visid(MJS, sort="b")

# sort via ME
visid(MJS, sort = "b", sort.method="ME")

# only rows, only columns
visid(MJS, sort = "r", sort.method="ME")
visid(MJS, sort = "c", sort.method="ME")

# sort via optile
visid(MJS, sort = "b", sort.method="optile")

visid(MJS, sort = "b", sort.method="optile", iter=10)

visid(MJS, sort = "b", sort.method="optile",fun="ca")

# 24 rows
visid(MJS, sort = "r", sort.method="optile", fr=24)


# 24 rows, >= 40% columns
visid(MJS, sort = "r", sort.method="optile", fr=24, fc = 0.4)


# zoom y marginal
visid(MJS, sort = "r", sort.method="optile", fr=24, s=1)


# zoom x marginal
visid(MJS, sort = "r", sort.method="optile", fr=24, pmax=0.1, s =0.5)

}

# NA-example: GeneEx

visna(GeneEx, sort = "b", sort.method="optile", fr=50, pmax=0.05, s = 2)

require(MASS)
visdf(housing)
visdf(housing,opts=list(var.col="w"))

visdf(housing,opts=list(var.col="w"), sort="r")






}

