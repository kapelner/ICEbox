#' Clustering of ICE and d-ICE curves by kmeans.
#'
#' Clustering if ICE and d-ICE curves by kmeans. All curves are centered to have mean 0
#' and then kmeans is applied to the curves with the specified number of clusters.
#'
#' @param ice_obj Object of class \code{ice} or \code{dice} to cluster.
#' @param nClusters Number of clusters to find.
#' @param plot If \code{TRUE}, plots the clusters.
#' @param plot_margin Extra margin to pass to \code{ylim} as a fraction of the range of cluster centers.
#' @param colorvec Optional vector of colors to use for each cluster.
#' @param plot_pdp If \code{TRUE}, the PDP (\code{ice} object) or d-PDP (\code{dice} object)
#'   is plotted with a dotted black line and highlighted in yellow.
#' @param x_quantile If \code{TRUE}, the plot is drawn with the x-axis taken to be \code{quantile(gridpts)}. If \code{FALSE},
#'   the predictor's original scale is used.
#' @param avg_lwd Average line width to use when plotting the cluster means.  Line width is proportional to the cluster's
#'   size.
#' @param centered If \code{TRUE}, all cluster means are shifted to be to be 0 at the minimum value of the predictor.
#'   If \code{FALSE}, the original cluster means are used. 
#' @param plot_legend If \code{TRUE} a legend mapping line colors to the proportion of the data in each cluster is 
#'   added to the plot.
#' @param num_cores Integer number of cores to use for parallel operations. Default is 1.
#' @param main Optional title for the plot.
#' @param ... Additional arguments for plotting.
#'
#' @return A list with the following elements:
#'   \item{cl}{The output of the \code{kmeans} call (a list of class \code{kmeans}).}
#'   \item{plot}{The ggplot object used for plotting (if \code{plot = TRUE}).}
#'
#' @seealso \code{\link{ice}}, \code{\link{dice}}
#' @examples
#' \dontrun{
#' require(ICEbox)
#' require(randomForest)
#' require(MASS) #has Boston Housing data, Pima
#'
#' data(Boston) #Boston Housing data
#' X = Boston
#' y = X$medv
#' X$medv = NULL
#'
#' ## build a RF:
#' bh_rf = randomForest(X, y)
#'
#' ## Create an 'ice' object for the predictor "age":
#' bh.ice = ice(object = bh_rf, X = X, y = y, predictor = "age",
#'             frac_to_build = .1) 
#'
#' ## cluster the curves into 2 groups.
#' clusterICE(bh.ice, nClusters = 2, plot_legend = TRUE)
#'
#' ## cluster the curves into 3 groups, start all at 0.
#' clusterICE(bh.ice, nClusters = 3, plot_legend = TRUE, center = TRUE)
#' }
#' @export
clusterICE = function(ice_obj, nClusters, plot = TRUE, plot_margin = 0.05, colorvec, plot_pdp = FALSE,
			x_quantile = FALSE, avg_lwd = 3, centered = FALSE, plot_legend = FALSE, main = NULL, num_cores = 1, ...){

	DEFAULT_COLORVEC = c("green", "red", "blue", "black", "green", "yellow", "pink", "orange", "forestgreen", "grey")

	if(!inherits(ice_obj, "ice") && !inherits(ice_obj, "dice")){
		stop("'ice_obj' must be of class 'ice' or 'dice'.")
	}

	objclass = class(ice_obj)

	if (missing(nClusters) || !(nClusters %% 1 == 0 ) || (nClusters <= 0)){
		stop("nClusters must be a positive integer.")
	}

	#make all curves have avg value = 0.
	if(objclass == "ice"){
		curves = rowCenter_cpp(ice_obj$ice_curves, n_cores = num_cores)  #sum(ice_curves[i,]) = 0 for all i
	}
	else{
		curves = rowCenter_cpp(ice_obj$d_ice_curves, n_cores = num_cores)
	}

	#cluster
	cl = kmeans(curves, iter.max = 20, centers = nClusters, ...)
	if (missing(colorvec)){
		colorvec = DEFAULT_COLORVEC
		if(length(colorvec) < nClusters){
			colorvec = c(colorvec, rgb(runif(nClusters - 10, 0, 0.7), runif(nClusters - 10, 0, 0.7), runif(nClusters - 10, 0, 0.7)))
		}
	}

	cluster_centers = cl$centers

	if (centered){
		cluster_centers = cluster_centers - cluster_centers[, 1]
	}

	if (plot){
		#y limits
		rg = range(cluster_centers)
		dist = rg[2] - rg[1]
		rg_min = rg[1] - plot_margin * dist
		rg_max = rg[2] + plot_margin * dist

		#x grid and xlab
		xlab = ice_obj$xlab
		grid = ice_obj$gridpts
		ecdf_fcn = NULL
		if (x_quantile){
			xlab = paste("quantile(", xlab, ")", sep = "")
			ecdf_fcn = ecdf(grid)
			grid = ecdf_fcn(grid)
		}

		cluster_order = order(cl$centers[, 1]) #use original, non-centered object only
		cluster_size = cl$size / sum(cl$size)
		total_line_width = avg_lwd * nClusters

		plot_centers = cluster_centers[cluster_order, , drop = FALSE]
		plot_sizes = cluster_size[cluster_order]
		plot_ids = seq_len(nrow(plot_centers))
		plot_labels = as.character(round(plot_sizes, 2))
		colorvec = colorvec[seq_len(nClusters)]

		centers_df = data.table(
			cluster_id = factor(rep(plot_ids, each = length(grid)), levels = plot_ids),
			x = rep(grid, times = nClusters),
			y = as.vector(t(plot_centers)),
			line_size = rep(plot_sizes * total_line_width, each = length(grid))
		)

		lwd_scale = 0.5
		p = ggplot(
			centers_df,
			aes(x = x, y = y, group = cluster_id, color = cluster_id, linewidth = line_size * lwd_scale)
		) +
			geom_line() +
			scale_linewidth_identity(guide = "none") +
			scale_color_manual(
				values = setNames(colorvec, plot_ids),
				breaks = plot_ids,
				labels = plot_labels,
				name = "Prop."
			)

		if(plot_pdp){
			# this is done after all other plotting so nothing obfuscates the PDP
			if(objclass == "ice"){
				pdp = ice_obj$pdp
			}
			else{
				pdp = ice_obj$dpdp #really want the dpdp
			}
			#now mean center it:
			pdp = pdp - mean(pdp)

			#if centered=TRUE, start everything at 0.
			if(centered){
				pdp = pdp - pdp[1]
			}

			#calculate the line thickness based on how many lines there are
			#every add'l cluster we get 0.5 more highlight up to 8
			pdp_df = data.table(x = grid, y = pdp)
			p = p +
				geom_line(
					data = pdp_df,
					aes(x = x, y = y),
					color = "yellow",
					linewidth = lwd_scale * min(4 +  nClusters * 0.5, 5),
					inherit.aes = FALSE
				) +
				geom_line(
					data = pdp_df,
					aes(x = x, y = y),
					color = "BLACK",
					linewidth = lwd_scale * 3,
					linetype = "dashed",
					inherit.aes = FALSE
				)
		}#end of pdp

		if (ice_obj$nominal_axis){
			axis_vals = sort(unique(ice_obj$xj))
			if (x_quantile){
				axis_breaks = ecdf_fcn(axis_vals)
			} else {
				axis_breaks = axis_vals
			}
			p = p + scale_x_continuous(breaks = axis_breaks, labels = axis_vals)
		}

		p = p + scale_y_continuous(limits = c(rg_min, rg_max)) +
			labs(x = xlab, y = "cluster yhat", title = main) +
			theme_classic(base_size = 11) +
			theme(
				panel.border = element_rect(fill = NA, color = "black"),
				panel.grid = element_blank()
			)

		if (plot_legend){
			p = p + theme(
				legend.position = c(0.02, 0.98),
				legend.justification = c(0, 1)
			)
		} else {
			p = p + guides(color = "none")
		}

		print(p)
		return(invisible(list(cl = cl, plot = p)))
	}

	invisible(list(cl = cl, plot = NULL))
}
