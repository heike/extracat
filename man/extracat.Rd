\name{extracat}
\alias{extracat}

\docType{package}
\title{
Categorical Data Analysis and Visualization
}
\description{
This package offers a variety of functions that can be used for categorical data analysis or at least have to do with categorical data. 

Among the most interesting features are
\tabular{ll}{
\link[extracat]{rmb} \tab RMB plots visualize contingency tables.\cr
 \tab The function offers different visualizations of conditional distributions and their corresponding weights (frequencies)\cr
 \tab including multiple barcharts, spineplots and piecharts.  There are different ways of displaying the residuals from statistical models.\cr
 \link[extracat]{scpcp} \tab  A static version of CPCP plots with labeling and color options.\cr
\link[extracat]{optile} \tab The OPTILE interface was developed for the Google Summer of Code 2011.\cr
\tab It offers a variety of reordering techniques for contingency tables. \cr
\tab The reordering of the categories not only improves visualizations.\cr
\link[extracat]{rmbmat} \tab A matrix of RMB plots not unlike a Scatterplot matrix (SPLOM) is produced. \cr
\tab It uses binning for continuous data.\cr
\link[extracat]{hexpie} \tab After a hexagonal binning of x and y a third categorical target variable \cr
\tab is displayed via piecharts (or embedded hexagons) within the hexagons. 
This avoids problems with overplotting.\cr
\link[extracat]{quickfechner} \tab A very fast implementation of the fechnerian scaling technique, \cr
\tab which computes a fechnerian distance matrix from a (dis.)similarity matrix.\cr
\tab The technique is equivalent to a shortest path algorithm.\cr
\link[extracat]{fluctile} \tab Offers two- or multidimensional fluctuation diagrams and multiple barcharts.\cr
\link[extracat]{arsim} \tab Simulates categorical data arrays of any dimension. \cr
\tab The number of observations, the number of block clusters and the noise level and type can be chosen. \cr
}

}
\details{
\tabular{ll}{
Package: \tab extracat\cr
Type: \tab Package\cr
Version: \tab 1.6-4\cr
Date: \tab 2013-12-11\cr
License: \tab -\cr
LazyLoad: \tab yes\cr
}

}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany

Maintainer: Alexander Pilhoefer \email{alexander.pilhoefer@math.uni-augsburg.de}
}
\references{
Alexander Pilhoefer, Antony Unwin (2013). New Approaches in Visualization of Categorical Data: R Package
  extracat. Journal of Statistical Software, 53(7), 1-25. URL \url{http://www.jstatsoft.org/v53/i07/}
}

\keyword{ categorical }
\keyword{ parallel coordinates }
\keyword{ multiple barcharts }
\keyword{ mosaicplots }
\keyword{ interactive }

