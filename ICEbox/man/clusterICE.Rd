\name{clusterICE}
\alias{clusterICE}
\title{
Clustering of ICE and d-ICE curves by kmeans.
}
\description{
Clustering if ICE and d-ICE curves by kmeans. All curves are centered to have mean 0
and then kmeans is applied to the curves with the specified number of clusters.
}
\usage{
clusterICE(ice_obj, nClusters, plot = TRUE, plot_margin = 0.05, 
           colorvec, plot_pdp = FALSE, x_quantile = FALSE, 
           avg_lwd = 3, centered = FALSE, 
           plot_legend = FALSE, ...)
}

\arguments{
  \item{ice_obj}{
Object of class \code{ice} or \code{dice} to cluster.
}

\item{nClusters}{
Number of clusters to find.
}

\item{plot}{
If \code{TRUE}, plots the clusters.
}

\item{plot_margin}{
Extra margin to pass to \code{ylim} as a fraction of the range of cluster centers.
}

\item{colorvec}{
Optional vector of colors to use for each cluster.
}

\item{plot_pdp}{
If \code{TRUE}, the PDP (\code{ice} object) or d-PDP (\code{dice} object) 
is plotted with a dotted black line and highlighted in yellow.
}

\item{x_quantile}{
If \code{TRUE}, the plot is drawn with the x-axis taken to be \code{quantile(gridpts)}. If \code{FALSE},
the predictor's original scale is used.
}

\item{avg_lwd}{
Average line width to use when plotting the cluster means.  Line width is proportional to the cluster's
size.
}

\item{centered}{
If \code{TRUE}, all cluster means are shifted to be to be 0 at the minimum value of the predictor.
If \code{FALSE}, the original cluster means are used. 
}
\item{plot_legend}{
If \code{TRUE} a legend mapping line colors to the proportion of the data in each cluster is 
added to the plot.
}
\item{...}{Additional arguments for plotting.}
}

\value{
The ouput of the \code{kmeans} call (a list of class \code{kmeans}).
}

\seealso{
ice, dice
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
bh_rf = randomForest(X, y)

## Create an 'ice' object for the predictor "age":
bh.ice = ice(object = bh_rf, X = X, y = y, predictor = "age",
            frac_to_build = .1) 

## cluster the curves into 2 groups.
clusterICE(bh.ice, nClusters = 2, plot_legend = TRUE)

## cluster the curves into 3 groups, start all at 0.
clusterICE(bh.ice, nClusters = 3, plot_legend = TRUE, center = TRUE)
}
}
