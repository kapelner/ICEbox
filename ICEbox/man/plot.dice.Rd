\name{plot.dice}
\alias{plot.dice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a plot of a \code{dice} object.
}
\description{
Plotting of \code{dice} objects.
}
\usage{
\method{plot}{dice}(x, plot_margin = 0.05, frac_to_plot = 1, 
      plot_sd = TRUE, plot_orig_pts_deriv = TRUE, pts_preds_size = 1.5, 
      colorvec, color_by = NULL, x_quantile = FALSE, plot_dpdp = TRUE, 
      rug = TRUE, ...)
}
\arguments{
  \item{x}{
Object of class \code{dice} to plot.
}
  \item{plot_margin}{
Extra margin to pass to \code{ylim} as a fraction of the range of \code{x$d_ice_curves}.
}
  \item{frac_to_plot}{
If \code{frac_to_plot} is less than 1, randomly plot \code{frac_to_plot} fraction of the
curves in \code{x$d_ice_curves}.
}
  \item{plot_sd}{
If \code{TRUE}, plot the cross-observation sd of partial derivatives below the derivative plots.
}
  \item{plot_orig_pts_deriv}{
If \code{TRUE}, marks each curve at the location of the derivative estimate at the
location of \code{predictor} actually occurring in the data. If \code{FALSE}
no mark is drawn.
}
  \item{pts_preds_size}{
Size of points to make if \code{plot_orig_pts_deriv} is \code{TRUE}.
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
  \item{plot_dpdp}{
If \code{TRUE}, the estimated derivative of the PDP is plotted and highlighted in yellow.
}
  \item{rug}{
If \code{TRUE}, marks are drawn on the x-axis indicating the presence of actual observed values.
Forced to \code{FALSE} when \code{x_quantile} is set to \code{TRUE}.
}
  \item{\dots}{
Additional plotting arguments.
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
dice
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
bhd.ice = ice(object = bhd_rf_mod, X = X, y = y, predictor = "age", frac_to_build = .1) 

# estimate derivatives, then plot.
bhd.dice = dice(bhd.ice)
plot(bhd.dice)
}
}
