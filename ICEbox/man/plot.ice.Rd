\name{plot.ice}
\alias{plot.ice}
\title{
Plotting of \code{ice} objects.
}
\description{
Plotting of \code{ice} objects.
}
\usage{
\method{plot}{ice}(x, plot_margin = 0.05, frac_to_plot = 1, plot_orig_pts_preds = TRUE, 
     pts_preds_size = 1.5, colorvec, color_by = NULL, x_quantile = FALSE, plot_pdp = TRUE,
     centered = FALSE, rug = TRUE, prop_range_y = TRUE,
     centered_percentile = 0.01, ...)
}

\arguments{
  \item{x}{
Object of class \code{ice} to plot.
}
  \item{plot_margin}{
Extra margin to pass to \code{ylim} as a fraction of the range of \code{x$ice_curves}.
}
  \item{frac_to_plot}{
If \code{frac_to_plot} is less than 1, randomly plot \code{frac_to_plot} fraction of the
curves in \code{x$ice_curves}.
}
  \item{plot_orig_pts_preds}{
If \code{TRUE}, marks each curve at the location of the observation's actual fitted value. If \code{FALSE},
no mark is drawn.
}
  \item{pts_preds_size}{
Size of points to make if \code{plot_origin_pts_preds} is \code{TRUE}.
}
  \item{colorvec}{
Optional vector of colors to use for each curve.
}
  \item{color_by}{
Optional variable name (or column number) in \code{Xice} to color curves by. If the \code{color_by}
variable has 10 or fewer unique values, a discrete set of colors is used for each value and a legend is 
printed and returned. If there are more values, curves are colored from light to dark corresponding 
to low to high values of the variable specified by \code{color_by}.
}
  \item{x_quantile}{
If \code{TRUE}, the plot is drawn with the x-axis taken to be \code{quantile(gridpts)}. If \code{FALSE},
the predictor's original scale is used.
}
  \item{plot_pdp}{
If \code{TRUE}, the PDP is plotted and highlighted in yellow.
}
  \item{centered}{
If \code{TRUE}, all curves are re-centered to be 0 at the quantile given by \code{centered_percentile}.
See Goldstein et al (2013) for details and examples. If \code{FALSE}, the original \code{ice_curves} are plotted. 
}
  \item{rug}{
If \code{TRUE}, marks are drawn on the x-axis indicating the presence of actual observed values.
Forced to \code{FALSE} when \code{x_quantile} is set to \code{TRUE}.
}
  \item{prop_range_y}{
When \code{TRUE} and \code{centered=TRUE} as well, the right vertical axis displays the
centered values as a fraction of the response's range.  See Goldstein et al (2013) for details
and examples.
}
  \item{centered_percentile}{
The percentile of \code{predictor} for which all \code{ice_curves} are "pinched together" and set to be 0.
Default is .01.
}
  \item{\dots}{
Other arguments to be passed to the \code{plot} function.
}
}
\value{
A list with the following elements.
  \item{plot_points_indices}{Row numbers of \code{Xice} of those observations presented in the plot.}
  \item{legend_text}{If the \code{color_by} argument was used, 
		a legend describing the map between the \code{color_by} predictor
        and curve colors.}
}

\seealso{
ice
}
\examples{
\dontrun{
require(ICEbox)
require(randomForest)
require(MASS) #has Boston Housing data, Pima

data(Boston) #Boston Housing data
X = Boston
y = X$medv
X$medv = NULL

## build a RF:
bhd_rf_mod = randomForest(X, y)

## Create an 'ice' object for the predictor "age":
bhd.ice = ice(object = bhd_rf_mod, X = X, y = y, predictor = "age",
            frac_to_build = .1) 

## plot
plot(bhd.ice, x_quantile = TRUE, plot_pdp = TRUE, frac_to_plot = 1) 

## centered plot
plot(bhd.ice, x_quantile = TRUE, plot_pdp = TRUE, frac_to_plot = 1
            centered = TRUE) 

## color the curves by high and low values of 'rm'.
# First create an indicator variable which is 1 if the number of 
# rooms is greater than the median:
median_rm = median(X$rm)
bhd.ice$Xice$I_rm = ifelse(bhd.ice$Xice$rm > median_rm, 1, 0)  

plot(bhd.ice, frac_to_plot = 1, centered = TRUE, prop_range_y = TRUE,  
            x_quantile = T, plot_orig_pts_preds = T, color_by = "I_rm")
}
}
