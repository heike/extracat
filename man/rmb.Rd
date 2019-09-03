\name{rmb}
\alias{rmb}
\alias{rmb.formula}
\alias{rmb.table}
\alias{rmb.ftable}
\alias{rmb.glm}
\concept{Multiple Barchart for relative frequencies and generalized Spineplots}
\title{
  Multiple Barchart for relative frequencies and generalized Spineplots
}
\description{
  The \code{rmb} function basically produces a Multiple Barchart for the relative frequencies of some target categories within each combination of the explanatory variables. The weights of those combinations (that is the absolute frequencies) are represented in the total with of the corresponding barchart.
  The result is a graphic which allows to read the conditional target distributions exactly from the graphic without losing the information about the importance (in the sense of the number of observations) of the different combinations.
 
 Additionally the \code{rmb} function allows to draw spineplots instead of the barcharts within each explanatory combination. On that score it can be seen as a generalization of Spineplots.
}
\usage{
\method{rmb}{formula}(formula, data, col.vars = NULL, spine = FALSE,
 circular = FALSE, eqwidth = FALSE, cat.ord = NULL,   cut = NULL,
 innerval = 1,   freq.trans = NULL, num.mode = FALSE, max.scale = 1,
 use.na = FALSE,  expected = NULL, residuals = NULL, model.opt = list(),
 gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,
  label.opt = list(),  vp = NULL, \dots)
\method{rmb}{ftable}(x, col.vars = NULL, spine = FALSE, circular = FALSE,
eqwidth = FALSE, cat.ord = NULL,  freq.trans = NULL, max.scale = 1,
use.na = FALSE, expected = NULL, residuals = NULL, model.opt = list(),
gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,
label.opt = list(), vp = NULL, \dots)
}
\arguments{
 \item{x}{ Either a table or a model of class \code{"glm"} and family \code{"poisson"} or \code{"binomial"}.
 A table must be either of class \code{table} or of class \code{ftable}. The latter also implicitly defines the the order in which the variables will be added to the plot. The arguments formula and data will be omitted. Please note that the model based version is still beta and will be improved in a future release.
 }
  \item{formula}{The \code{\link[stats]{formula}} specifying the variables in their given order with the last variable being the target variable. The left hand side defines a weighting variable. If the weights are frequencies in a variable called \code{"Freq"} this is detected automatically if no other variable is defined. }
  \item{data}{The dataset as a  data.frame or \code{\link[stats]{ftable}}.}
  
  \item{col.vars}{Logical vector with split directions where \code{TRUE} stands for horizontal splitting. The last (target) variable is always arranged on the x-axis. }
  \item{spine}{If \code{TRUE} a spineplot will be drawn instead of each barchart. This is recommended for binary target variables.}
 \item{circular}{ If \code{TRUE} a piechart will be drawn instead of each barchart. \code{spine} is set to \code{FALSE}.}
  \item{eqwidth}{If \code{TRUE} the bar length of the multiple barchart in the background no longer restricts the width of the barcharts/spineplots for the relative frequencies of the target variable.}
   \item{cat.ord}{A vector specifying the categories of the target variable which will be visualized in the specified order. The default is to use all categories.}
    \item{cut}{Numeric variables will be cut into this number of intervals. May also be a vector with specifications for each variable.}
        \item{innerval}{The function \link[extracat]{innerval} is used to reduce numeric variables to an interval which is symmetric around the median  contains the specified proportion of observations (or as close to this as possible).}
   \item{freq.trans}{This parameter allows to transform the absolute frequencies used for the underlying multiple barchart. Possible values are \code{"log"}, \code{"sqrt"} or \code{list("sqrt",k)}. The latter stands for the k-th root transformation.}
   \item{num.mode}{ In the numeric mode the gaps are removed and axes typical for numeric variables are drawn. Ignored for factor variables.}
  \item{max.scale}{The maximum value of the probability (y-axis) scale for each combination. Unsurprisingly the default is 1. The axis will be drawn if \code{yaxis} is \code{TRUE}.}
      
  \item{use.na}{If \code{TRUE} missing values will be changed to a level \code{"N/A"} and else (which is the default) the function \link[stats]{na.omit} will be called to reduce the dataset to complete cases only.}
  \item{expected}{ There are three possibilities how to specify this parameter:\cr
  
  1. A list of integer vectors denoting the interaction terms in the poisson or proportional odds model, e.g. \code{list( c(1,2,3), c(1,4) )} for all interactions between variables 1,2 and 3 as well as between 1 and 4.\cr
  2. A logical indicating whether or not to use a model (logit independence model).\cr
  3. A vector with expected values, e.g. from a model. If \code{residuals} remains undefined the response residuals will be plotted.\cr
  \cr
  
  If undefinded or set to \code{FALSE} only the observed values will be plotted.}
 \item{residuals}{ If \code{expected} is a vector with expected values it is also possible to specify residuals. This is used internally by \code{rmb.glm}.}
 \item{ model.opt }{A list with optional parameters for model specifications. Possible parameters are:
 	\tabular{rlll}{
         \code{use.ecpected.values}   \tab A logical specifying whether or not to use the frequencies \cr
         \tab predicted by the model instead of the observed frequencies.\cr
      	 \code{mod.type} \tab Either \code{"poisson"} or \code{"polr"}.\cr
	 \tab See \link[stats]{glm} and \link[MASS]{polr}.	\cr
         \code{resid.type} \tab \code{"pearson"}, \code{"deviance"}, \code{working}, \code{partial} \cr
         \tab or  \code{"response"}. For \code{polr} models only the latter is available.		\cr
         \code{resid.display} \tab One of \code{"values", "color"} or \code{"both"}. \code{"values"} will result in bars \cr 
         \tab or wedges for both expected and observed frequencies. Hence the raw residuals are shown in the graphic.\cr
         \tab  \code{"color"} will set the \code{col} argument aside and use colors on a red-blue-scale \cr
         \tab to represent (pearson) residuals. \code{"both"} does both.		\cr
         \code{max.rat} \tab 	If a model is specified and \code{resid.display = "both"} the x-scales \cr
         \tab will not be reduced to less than \code{1/max.rat}:\cr
         \tab The x-scales are reduced whenever an observed frequency exceeds the maximal scale.	\cr
        }
}

  \item{gap.prop}{The maximum proportion of the total plot width which is used for the gaps.}
  \item{gap.mult}{The incremental multiplier for the gaps of different dimensions. The gaps corresponding to any one variable are \code{gap.mult} times larger than those corresponding to the next variable on the same axis.}
  \item{col}{Either a vector defining the colors of the bars or a name specifying a palette: \code{"hsv"} and \code{"rgb"} for hsv-based rainbow colors, \code{"hcl"} for hcl-based rainbow colors (default), \code{"div"} or \code{"diverge"} for hcl-based diverging colors and finally \code{"seq"} or \code{"sequential"}
  for hcl-based sequential colors.Additional arguments can be specified via the \code{col.opt} argument according to the underlying functions in the colorspace package, e.g. \link[colorspace]{rainbow_hcl}. For the hsv-based colors see \link[grDevices]{rainbow}. Specifying a color or palette has no effect if an \code{expected} model is defined.}
  \item{col.opt}{Further options for the color palettes. See e.g. \link[colorspace]{rainbow_hcl} or \link[grDevices]{rainbow}.
  Other parameters are:
  \tabular{rlll}{
  \code{col2}\tab for the color of the background/weight bars, \cr
  \code{line.col}\tab for the color of all lines (bars, rectangles), \cr
  \code{bg}\tab for the background color of the whole graphic, \cr
  \code{bgs}\tab for the background color of each tile \cr
  }
  }

  \item{label}{Either a logical specifying whether or not to draw labels or a numeric vector defining which variables shall be labelled. }
 \item{ label.opt }{A list with optional parameters for label specifications. Possible parameters are:
 	\tabular{rlll}{
         \code{yaxis}   \tab If \code{TRUE} a vertical axis will be drawn at both sides of the plot. \cr
         \tab This is recommended when changing the \code{max.scale} argument.\cr
      	 \code{boxes} \tab Should the labels be surrounded by boxes?	\cr
         \code{lab.tv} \tab Should the target variable be included in the labeling?		\cr
         \code{varnames} \tab Should the variable names be shown as labels?	\cr
         \code{abbrev} \tab 	An single integer value or a vector of integer values specifying the number \cr
         \tab of characters to which the labels will automatically be abbreviated. \cr
         	 \code{lab.cex} \tab The fontsize multplier.	\cr
        }
}
 \item{vp}{An optional viewport to plot in. \code{vp = c(i, j)} can be used as a shortcut to \code{viewport(layout.pos.row = i, layout.pos.col = j)} }
   \item{\dots}{further arguments. Usually not necessary.}
}
\details{
  A similar way to regard the graphic is the following: A Multiple Barchart of the explanatory variables is drawn with bars in horizontal direction. Then within each of the resulting bars a barchart of the conditional distribution of the target variable is drawn with bars in vertical direction.
}
\value{
 \code{invisible(TRUE)}
}
\references{
Alexander Pilhoefer, Antony Unwin (2013). New Approaches in Visualization of Categorical Data: R Package
  extracat. Journal of Statistical Software, 53(7), 1-25. URL \url{http://www.jstatsoft.org/v53/i07/}
}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}

\seealso{
  \code{\link[graphics]{mosaicplot}} 
}
\examples{
  require(MASS)
    # simple example
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, gap.mult = 2,
        col.vars = c(FALSE,TRUE,TRUE,FALSE), label.opt = list(abbrev = 3))
    
    # with sqrt-transformation and horizontal splits only
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, gap.mult = 2,
        col.vars = c(TRUE,TRUE,TRUE,TRUE), freq.trans = "sqrt",
           label.opt = list(abbrev = 3) )
    
    # a generalized spineplot with the first category highlighted
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, spine = TRUE, 
        cat.ord = 1, mult = 2, col.vars = c(1,3,4), 
        freq.trans = list("sqrt",3),  label.opt = list(abbrev = 2))
  \dontrun{  
    # a generalized spineplot with all categories highlighted 
    # in a changed order
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, spine = TRUE,
        cat.ord = c(3,1,2), gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE),
        freq.trans = "sqrt",  label.opt = list(abbrev = 3))
    
    # the barchart version only for categories 1 and 3
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, 
        cat.ord = c(1,3), gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE),
        freq.trans = "sqrt",  label.opt = list(abbrev = c(4,1,1,1)))
        
        
    # with equal widths
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, eqwidth = TRUE,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
         label.opt = list(abbrev = 2, lab.tv = TRUE))
    
    # ----- models and residuals ----- #
    # using the logistic model: Sat by Type only
    
    #   residual shadings and expected values
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), expected = list(c(1,2,3),c(1,4)),
        model.opt = list(use.expected.values = TRUE, resid.display = "color") )
       
    #   residual values without shadings
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), expected = list(c(1,2,3),c(1,4)),
        model.opt = list( resid.display = "values") )
    
    #   residual shadings and expected values
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), expected = list(c(1,2,3),c(1,4)),
        model.opt = list(use.expected.values = TRUE, resid.display = "color") )
   
    #   barcharts with residual shadings and values
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), expected = list(c(1,2,3),c(1,4)) )
        
    #   spineplots with residual shadings and values 
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, spine = TRUE,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), expected = list(c(1,2,3),c(1,4)) )
        
    #   piecharts with residual shadings and values
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, circular = TRUE,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), expected = list(c(1,2,3),c(1,4)) )
        
     # ----- using an ftable to create the plot ----- #
     tt = xtabs(Freq~Type+Cont+Infl+Sat, data = housing)
     ft = ftable(tt, col.vars= c(1,4))
      rmb(tt)
      rmb(ft)

    # ----- using a glm model ----- #
    mod1 <- glm(Freq ~ Type*Infl*Cont + Type*Sat, data = housing, family = poisson)
    rmb(mod1, circular = TRUE,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), 
        label.opt = list(abbrev = 3), model.opt = list(use.expected.values = TRUE) )
        
        
      # ----- the numeric mode and cuts ----- #   
     data(olives)
     # only three cuts to show how it works
     rmb(~palmitoleic+stearic+Region, data = olives, cut = c(3,3,0))
     
     require(ggplot2)
     data(diamonds)
     diamonds$lprice <- log(diamonds$price)
     # a minority of extreme observations mess the display up:
     rmb(~depth+table+lprice, data = diamonds, eqwidth = TRUE, spine = TRUE, 
     cut = c(36,36,5), col = "seq", num.mode = TRUE)
     
     # we can zoom in via innerval:
      rmb(~depth+table+lprice, data = diamonds,  circular = TRUE, 
      cut = c(36,36,5), col = "div", innerval = 0.95, 
      num.mode = TRUE, freq.trans ="log")
      
           # price, carat and color
        diamonds$lprice <- log(diamonds$price)
        diamonds$lcarat <- log(diamonds$carat)     
      rmb(~lcarat+lprice+color, data = diamonds,
      cut = c(24,24,0), col = "rgb", num.mode = TRUE,
       freq.trans="sqrt", eqwidth=TRUE, max.scale=0.5)
       }
}

