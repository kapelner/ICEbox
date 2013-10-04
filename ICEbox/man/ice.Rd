\name{ice}
\alias{ice}

\title{
Creates an object of class \code{ice}.
}
\description{
Creates an \code{ice} object with individual conditional expectation curves
for the passed model object, \code{X} matrix, predictor, and response.  See
Goldstein et al (2013) for further details.
}
\usage{
ice(object, X, y, predictor, predictfcn, verbose = TRUE, frac_to_build = 1, 
             indices_to_build = NULL, num_grid_pts, logodds = F, ...)
}

\arguments{
  \item{object}{
The fitted model to estimate ICE curves for.
}
  \item{X}{
The design matrix we wish to estimate ICE curves for. Rows are observations, columns are
predictors.  Typically this is taken to be \code{object}'s training data, but this is not
strictly necessary.
}
  \item{y}{
Optional vector of the response values \code{object} was trained on.  It is used
to compute y-axis ranges that are useful for plotting.  If not passed, the range
of predicted values is used and a warning is printed.
}
  \item{predictor}{
The column number or variable name in \code{X} of the predictor of interest, 
(\eqn{x_S = X[, j]}{x_S= X[, j]}).
}
  \item{predictfcn}{
Optional function that accepts two arguments, \code{object} and \code{newdata}, and
returns an \code{N} vector of \code{object}'s predicted response for data \code{newdata}. 
If this argument is not passed, the procedure attempts to find a generic \code{predict} 
function corresponding to \code{class(object)}.
}
  \item{verbose}{
If \code{TRUE}, prints messages about the procedure's progress.
}
  \item{frac_to_build}{
Number between 0 and 1, with 1 as default.  For large \code{X} matrices or fitted models
that are slow to make predictions, specifying \code{frac_to_build} less than 1 will choose
a subset of the observations to build curves for. The subset is chosen such that the remaining
observations' values of \code{predictor} are evenly spaced throughout the quantiles of the
full \code{X[,predictor]} vector.  
}
  \item{indices_to_build}{
Vector of indices, \eqn{\subset \{1, \ldots, nrow(X)\}}{each element between \code{1} and \code{nrow(X)}} specifying which observations to build ICE curves for. As this is an alternative to setting \code{frac_to_build}, both 
cannot be specified.
}
  \item{num_grid_pts}{
Optional number of values in the range of \code{predictor} at which to estimate each curve. 
If missing, the curves are estimated at each unique value of \code{predictor} 
in the \code{X} observations we estimate ICE curves for.
}
  \item{logodds}{
If \code{TRUE}, for classification creates PDPs by plotting the centered log-odds of implied by the fitted 
probabilities. We assume that the generic or passed predict function 
returns probabilities, and so the flag tells us to transform these to centered logits after 
the predictions are generated. 
}
  \item{\dots}{
Other arguments to be passed to \code{object}'s generic predict function.  
}
}
\value{
A list of class \code{ice} with the following elements.
\item{gridpts}{Sorted values of \code{predictor} at which each curve is estimated. Duplicates
are removed -- by definition, elements of \code{gridpts} are unique.}
\item{ice_curves}{Matrix of dimension \code{nrow(X)} by \code{length(gridpts)}.
Each row corresponds to an observation's ICE curve, estimated at the values of \code{predictor} in 
\code{gridpts}.}
\item{xj}{The actual values of \code{predictor} observed in the data in the order
of \code{Xice}.}
\item{actual_predictions}{Vector of length \code{nrow(X)} containing the model's 
predictions at the actual value of the predictors in the order of \code{Xice}.}
\item{xlab}{String with the predictor name corresponding to \code{predictor}. If \code{predictor}
is a column number, \code{xlab} is set to \code{colnames(X)[, predictor]}.}
\item{nominal_axis}{If \code{TRUE}, \code{length(gridpts)} is 5 or fewer; otherwise \code{FALSE}.
When \code{TRUE} the \code{plot} function treats the x-axis as if x is nominal.}
\item{range_y}{If \code{y} was passed, the range of the response.  Otherwise it is defaulted to be
\code{max(ice_curves)-min(ice_curves)} and a message is printed to the console.}
\item{Xice}{A matrix containing the subset of \code{X} for which ICE curves are estimated. 
Observations are ordered to be increasing in \code{predictor}.  This ordering is the same one 
as in \code{ice_curves}, \code{xj} and \code{actual_predictions}, meaning for all these objects
the \code{i}-th element refers to the same observation in \code{X}.}
\item{pdp}{A vector of size \code{length(gridpts)} which is a numerical approximation to the partial 
dependence function (PDP) corresponding to the estimated ICE curves. See Goldstein et al (2013) for a discussion
of how the PDP is a form of post-processing. See Friedman (2001) for a description of PDPs.}
\item{predictor}{Same as the argument, see argument description.}
\item{logodds}{Same as the argument, see argument description.}
\item{indices_to_build}{Same as the argument, see argument description.}
\item{frac_to_build}{Same as the argument, see argument description.}
\item{predictfcn}{Same as the argument, see argument description.}
}
\references{
Jerome Friedman. Greedy Function Approximation: A Gradient Boosting Machine. The Annals of Statistics,
29(5): 1189-1232, 2001.

A Goldstein, A Kapelner, J Bleich and E Pitkin. Peeking Inside the Black Box: Visualizing Statistical
Learning with Plots of Individual Conditional Expectation. ArXiv e-prints, 2013.
}

\seealso{
plot.ice, print.ice, summary.ice
}
\examples{
\dontrun{
require(ICEbox)
require(randomForest)
require(MASS) #has Boston Housing data, Pima

########  regression example
data(Boston) #Boston Housing data
X = Boston
y = X$medv
X$medv = NULL

## build a RF:
bhd_rf_mod = randomForest(X, y)

## Create an 'ice' object for the predictor "age":
bhd.ice = ice(object = bhd_rf_mod, X = X, y = y, predictor = "age", frac_to_build = .1) 

#### classification example
data(Pima.te)  #Pima Indians diabetes classification
y = Pima.te$type
X = Pima.te
X$type = NULL

## build a RF:
pima_rf_mod = randomForest(x = X, y = y)

## Create an 'ice' object for the predictor "skin":
# For classification we plot the centered log-odds. If we pass a predict
# function that returns fitted probabilities, setting logodds = TRUE instructs
# the function to set each ice curve to the centered log-odds of the fitted 
# probability.
pima.ice = ice(object = pima_rf_mod, X = X, predictor = "skin", logodds = TRUE,
                    predictfcn = function(object, newdata){ 
                         predict(object, newdata, type = "prob")[, 2]
                    }
              )

}
}
