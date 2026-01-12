#' Plotting of \code{ice} objects.
#'
#' Plotting of \code{ice} objects.
#'
#' @param x Object of class \code{ice} to plot.
#' @param plot_margin Extra margin to pass to \code{ylim} as a fraction of the range of \code{x$ice_curves}.
#' @param frac_to_plot If \code{frac_to_plot} is less than 1, randomly plot \code{frac_to_plot} fraction of the
#'   curves in \code{x$ice_curves}.
#' @param plot_points_indices If not \code{NULL}, this plots only the indices of interest. If not \code{NULL}, \code{frac_to_plot} must be 1 otherwise
#'   an error is thrown. Default is \code{NULL}.
#' @param plot_orig_pts_preds If \code{TRUE}, marks each curve at the location of the observation's actual fitted value. If \code{FALSE},
#'   no mark is drawn.
#' @param pts_preds_size Size of points to make if \code{plot_origin_pts_preds} is \code{TRUE}.
#' @param colorvec Optional vector of colors to use for each curve.
#' @param color_by Optional variable name in \code{Xice}, column number in \code{Xice}, or data vector of the correct length to color curves by.
#'   If the \code{color_by}
#'   variable has 10 or fewer unique values, a discrete set of colors is used for each value and a legend is
#'   printed and returned. If there are more values, curves are colored from light to dark corresponding
#'   to low to high values of the variable specified by \code{color_by}.
#' @param x_quantile If \code{TRUE}, the plot is drawn with the x-axis taken to be \code{quantile(gridpts)}. If \code{FALSE},
#'   the predictor's original scale is used.
#' @param plot_pdp If \code{TRUE}, the PDP is plotted and highlighted in yellow.
#' @param centered If \code{TRUE}, all curves are re-centered to be 0 at the quantile given by
#'   \code{centered_percentile}.
#'   See Goldstein et al (2013) for details and examples. If \code{FALSE}, the original \code{ice_curves} are plotted.
#' @param prop_range_y When \code{TRUE} and \code{centered=TRUE} as well, the range of the right vertical axis displays the
#'   centered values as a fraction of the sd of the fitted values on actual observations if \code{prop_type}
#'   is missing or set to \code{"sd"}. If \code{prop_type} is set to \code{"range"}, the right axis displays the
#'   centered values as a fraction of the range of the fitted values over the actual observations.
#' @param rug_quantile If not \code{NULL}, tick marks are drawn on the x-axis corresponding to the vector of quantiles specified by this parameter.
#'   Forced to \code{NULL} when \code{x_quantile} is set to \code{TRUE}.
#' @param centered_percentile The percentile of \code{predictor} for which all \code{ice_curves} are "pinched together" and set to be 0.
#'   Default is 0.
#' @param point_labels If not \code{NULL}, labels to plot next to each point. Default is \code{NULL}.
#' @param point_labels_size If not \code{NULL}, size of labels to plot next to each point. Default is \code{NULL} which means it's the size of \code{pts_preds_size}.
#' @param prop_type Scaling factor for the right vertical axis in centered plots if \code{prop_range_y} is \code{TRUE}. Can be one of
#'   \code{"sd"} (default) or \code{"range"}. Ignored if \code{centered} and \code{prop_range_y} are not both \code{TRUE}.
#' @param verbose If \code{TRUE}, prints the color legend to the console.
#' @param ... Other arguments to be passed to the \code{plot} function.
#' @param num_cores Used for parallel plotting speedup. Default is 1.
#'
#' @return A list with the following elements.
#'   \item{plot_points_indices}{Row numbers of \code{Xice} of those observations presented in the plot.}
#'   \item{legend_text}{If the \code{color_by} argument was used,
#'   a legend describing the map between the \code{color_by} predictor
#'   and curve colors.}
#'
#' @seealso \code{\link{ice}}
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
#' bhd_rf_mod = randomForest(X, y)
#'
#' ## Create an 'ice' object for the predictor "age":
#' bhd.ice = ice(object = bhd_rf_mod, X = X, y = y, predictor = "age",
#'             frac_to_build = .1)
#'
#' ## plot
#' plot(bhd.ice, x_quantile = TRUE, plot_pdp = TRUE, frac_to_plot = 1)
#'
#' ## centered plot
#' plot(bhd.ice, x_quantile = TRUE, plot_pdp = TRUE, frac_to_plot = 1,
#'		centered = TRUE)
#'
#' ## color the curves by high and low values of 'rm'.
#' # First create an indicator variable which is 1 if the number of
#' # rooms is greater than the median:
#' median_rm = median(X$rm)
#' bhd.ice$Xice$I_rm = ifelse(bhd.ice$Xice$rm > median_rm, 1, 0)
#'
#' plot(bhd.ice, frac_to_plot = 1, centered = TRUE, prop_range_y = TRUE,
#'             x_quantile = T, plot_orig_pts_preds = T, color_by = "I_rm")
#' bhd.ice = ice(object = bhd_rf_mod, X = X, y = y, predictor = "age",
#'             frac_to_build = 1)
#' plot(bhd.ice, frac_to_plot = 1, centered = TRUE, prop_range_y = TRUE,
#'             x_quantile = T, plot_orig_pts_preds = T, color_by = y)
#' }
#' @export
plot.ice = function(x, plot_margin = 0.05, frac_to_plot = 1, plot_points_indices = NULL,
					plot_orig_pts_preds = TRUE, pts_preds_size = 1.5,
					colorvec, color_by = NULL, x_quantile = TRUE, plot_pdp = TRUE, centered = FALSE, 
					prop_range_y = TRUE, rug_quantile = seq(from = 0, to = 1, by = 0.1), 
					centered_percentile = 0, point_labels = NULL, point_labels_size = NULL,
					prop_type = "sd", verbose = TRUE, num_cores = 1, ...){

	DEFAULT_COLORVEC = c("firebrick3", "dodgerblue3", "gold1", "darkorchid4", "orange4", "forestgreen", "grey", "black", "deeppink", "turquoise4")
	#think of x as x. needs to be 'x' to match R's generic.

	#some argument checking
	if (!inherits(x, "ice")){
		stop("object is not of class \"ice\"")
	}
	if (frac_to_plot <= 0 || frac_to_plot > 1 ){
		stop("frac_to_plot must be in (0,1]")
	}
	if(!(prop_type %in% c("sd","range"))){
		stop("prop_type must be either 'sd' or 'range'")
	}

	#extract the grid and lines to plot
	grid = x$gridpts
	n_grid = length(grid)
	ecdf_fcn = NULL
	if (x_quantile){
		ecdf_fcn = ecdf(grid)
		grid = ecdf_fcn(grid)
	}
	ice_curves = x$ice_curves
	N = nrow(ice_curves)

	if (!is.null(point_labels)){
		if (length(point_labels) != N){
			stop("point_labels must be same length as number of ICE curves: ", N)
		}
	}

	#### figure out the colorvec.
	legend_text = NULL #default is no legend.
	
	# Case 1: Neither colorvec nor color_by provided
	if (missing(colorvec) && missing(color_by)){
		#we're going to choose dark grey and randomly alpha the lines
		colorvec = sort(rgb(rep(0.4, N), rep(0.4, N), rep(0.4, N), runif(N, 0.4, 0.8)))
	} else if (!missing(color_by)) {
		# Case 2: color_by provided (palette mode if colorvec is provided)
		
		# argument checking first:
		arg_type = class(color_by)
		if(!(arg_type %in% c("character", "numeric", "factor"))){
			stop("color_by must be a column name in X or a column index")
		}
		if(inherits(color_by, "character")){
			if(!(color_by %in% names(x$Xice))){
				stop("The predictor name given by color_by was not found in the X matrix")
			}
			x_color_by = x$Xice[[color_by]]
		} else if (length(color_by) > N){ #tell the user the thing they passed in doesn't line up
			stop("The color_by_data vector you passed in has ", length(color_by), " entries but the ICEbox object only has ", N, " curves.")
		} else if (length(color_by) == N){ #it's an actual data vector
			x_color_by = color_by
		}
		else{  #check numeric
			if( color_by < 1 || color_by > ncol(x$Xice) || (color_by%%1 !=0)){
				stop("color_by must be a column name in X or a column index")
			}
			x_color_by = x$Xice[[color_by]]
		}
		x_unique = sort(unique(x_color_by)) #if you don't sort, it can switch...
		num_x_color_by = length(x_unique)

		# if there are 10 or fewer unique values of this x value, we use the
		# same color in DEFAULT_COLORVEC for each. Otherwise, we use a rainbow.
		if (num_x_color_by <= 10){
			if (missing(colorvec)){
				palette = DEFAULT_COLORVEC[1 : num_x_color_by]
			} else {
				if (length(colorvec) < num_x_color_by) {
					stop("color vector has length ", length(colorvec), " but there are ", num_x_color_by, " unique values of color_by")
				}
				palette = colorvec[1 : num_x_color_by]
			}
			which_category = match(x_color_by, x_unique)
			colorvec = palette[which_category]

			#now make the legend.
			legend_text = as.data.frame(cbind(x_unique, palette))
			x_column_name = ifelse(length(color_by) == N, "data vector level", ifelse(is.character(color_by), color_by, paste("x_", color_by, sep = "")))
			names(legend_text) = c(x_column_name,"color")
			if (verbose){
				cat("ICE Plot Color Legend\n")
				print(legend_text, row.names = FALSE)
			}
		} else {
			if (!missing(colorvec)) {
				if (length(colorvec) < N) {
					stop("color vector has length ", length(colorvec), " but there are ", N, " lines to plot")
				}
				# Use colorvec as is (assumed to be of length N)
			} else {
				if (is.factor(x_color_by)){
					warning("color_by is a factor with greater than 10 levels: coercing to numeric.")
					x_color_by = as.numeric(x_color_by)
				}

				alpha_blend_colors = matrix(0, nrow = N, ncol = 3)
				alpha_blend_colors[, 1] = seq(from = 1, to = 0, length.out = N)
				alpha_blend_colors[, 2] = seq(from = 0, to = 1, length.out = N)
				alpha_blend_colors[, 3] = 0

				rgbs = rgb(alpha_blend_colors[, 1], alpha_blend_colors[, 2], alpha_blend_colors[, 3])
				colorvec = rgbs[order(x_color_by)]
				if (verbose) cat("ICE Plot Color Legend: red = low values of the color_by variable and green = high values\n")
			}
		}
	} else {
		# Case 3: Only colorvec provided
		if (length(colorvec) < N) {
			stop("color vector has length ", length(colorvec), " but there are ", N, " lines to plot")
		}
	}


	#pull out a fraction of the lines to plot
	if (is.null(plot_points_indices)){
		all_ice_curves = x$ice_curves #all
		plot_points_indices = sample(1 : N, round(frac_to_plot * N))
		ice_curves = ice_curves[plot_points_indices, ]
	} else {
		if (frac_to_plot < 1){
			stop("frac_to_plot has to be 1 when plot_points_indices is passed to the plot function.")
		}
		ice_curves = ice_curves[plot_points_indices, ]
		all_ice_curves = ice_curves
	}

	if (nrow(ice_curves) == 0){
		stop("no rows selected: frac_to_plot too small.")
	}
	if (centered){
		centering_vector = ice_curves[, max(ncol(ice_curves) * centered_percentile, 1)]
		ice_curves = ice_curves - centering_vector
	}
	colorvec = colorvec[plot_points_indices]

	##### now start plotting
	min_ice_curves = min(ice_curves)
	max_ice_curves = max(ice_curves)
	range_ice_curves = max_ice_curves - min_ice_curves
	min_ice_curves = min_ice_curves - plot_margin * range_ice_curves
	max_ice_curves = max_ice_curves + plot_margin * range_ice_curves


	arg_list = list(...)

	#get the xlabel if it wasn't already passed explicitly.
	if( is.null(arg_list$xlab)){
		xlab = x$xlab

		#create some smart default for the x-label
		if (x_quantile){
			xlab = paste("quantile(", xlab, ")", sep = "")
		}
		if (!missing(color_by)){
			xlab = paste(xlab, "colored by", ifelse(length(color_by) == N, "a provided data vector", color_by))
		}
	} else {
		xlab = arg_list$xlab
	}


	#same for y label
	if( is.null(arg_list$ylab)){
		if (x$logodds){
			ylab = "partial log-odds"
		} else if(x$probit){
			ylab = "partial probit"
		}else {
			ylab = paste("partial yhat", ifelse(centered, "(centered)", ""))
		}
	} else {
		ylab = arg_list$ylab
	}

	#set ylim if not passed explicitly
	if (is.null(arg_list$ylim)){
		ylim = c(min_ice_curves, max_ice_curves)
	} else {
		ylim = arg_list$ylim
	}

	cex_axis = ifelse(is.null(arg_list$cex.axis), 1, arg_list$cex.axis)
	cex_lab = ifelse(is.null(arg_list$cex.lab), 1, arg_list$cex.lab)
	base_size = 11
	plot_theme = theme_classic(base_size = base_size) +
		theme(
			panel.border = element_rect(fill = NA, color = "black"),
			panel.grid = element_blank(),
			axis.text = element_text(size = base_size * cex_axis),
			axis.title = element_text(size = base_size * cex_lab),
			axis.text.y.right = element_text(size = base_size * cex_axis),
			axis.title.y.right = element_text(size = base_size * cex_lab)
		)
	if (!is.null(arg_list$xaxt) && arg_list$xaxt == "n"){
		plot_theme = plot_theme + theme(
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank()
		)
	}

	ice_df = data.table(
		curve_id = rep(seq_len(nrow(ice_curves)), each = n_grid),
		x = rep(grid, times = nrow(ice_curves)),
		y = melt_ice_curves_cpp(ice_curves, n_cores = num_cores),
		color = rep(colorvec, each = n_grid)
	)

	lwd_scale = 0.5
	p = ggplot(ice_df, aes(x = x, y = y, group = curve_id, color = color)) +
		geom_line() +
		scale_color_identity(guide = "none")

	points_df = NULL
	if (plot_orig_pts_preds || !is.null(point_labels)){
		yhat_actual = x$actual_prediction[plot_points_indices]
		if (centered){
			yhat_actual = yhat_actual - centering_vector
		}

		if (x_quantile){
			xj = ecdf_fcn(x$xj)[plot_points_indices]
		} else {
			xj = x$xj[plot_points_indices]
		}
		points_df = data.table(
			x = xj,
			y = yhat_actual,
			color = colorvec,
			label = if (is.null(point_labels)) NA else point_labels[plot_points_indices]
		)
	}

	if (plot_orig_pts_preds){
		p = p +
			geom_point(
				data = points_df,
				aes(x = x, y = y),
				color = rgb(0.1, 0.1, 0.1),
				size = pts_preds_size,
				inherit.aes = FALSE
			) +
			geom_point(
				data = points_df,
				aes(x = x, y = y, color = color),
				size = pts_preds_size * 0.7,
				inherit.aes = FALSE
			)
	}

	if (!is.null(point_labels)){
		label_size = ifelse(is.null(point_labels_size), pts_preds_size, point_labels_size)
		nudge_x = diff(range(grid))
		if (is.na(nudge_x) || nudge_x == 0){
			nudge_x = 0
		} else {
			nudge_x = nudge_x * 0.01
		}
		p = p +
			geom_text(
				data = points_df,
				aes(x = x, y = y, label = label),
				hjust = 0,
				nudge_x = nudge_x,
				size = label_size,
				inherit.aes = FALSE
			)
	}

	if (!is.null(rug_quantile) && !x_quantile){
		rug_df = data.table(x = as.numeric(quantile(x$xj, rug_quantile)))
		p = p +
			geom_rug(
				data = rug_df,
				aes(x = x),
				sides = "b",
				color = "blue4",
				inherit.aes = FALSE
			)
	}

	#if plot_pdp is true, plot actual pdp (in the sense of Friedman '01)
	#Ensure this is done after all other plotting so nothing obfuscates the PDP
	pdp = NULL
	if (plot_pdp){
		pdp = colMeans(all_ice_curves) # pdp = average over all the columns unless the user specifically specifies
		if (centered){
			pdp = pdp - pdp[ceiling(length(pdp) * centered_percentile + 0.00001)]
		}

		#calculate the line thickness based on how many lines there are
		num_lines = length(plot_points_indices)
		pdp_df = data.table(x = grid, y = pdp)
		p = p +
			geom_line(
				data = pdp_df,
				aes(x = x, y = y),
				color = "yellow",
				linewidth = lwd_scale * min(5.5 + (num_lines / 100) * 0.75, 8),
				inherit.aes = FALSE
			) +
			geom_line(
				data = pdp_df,
				aes(x = x, y = y),
				color = "BLACK",
				linewidth = lwd_scale * 4,
				inherit.aes = FALSE
			)
	}

	if (x$nominal_axis){
		axis_vals = sort(unique(x$xj))
		if (x_quantile){
			axis_breaks = ecdf_fcn(axis_vals)
		} else {
			axis_breaks = axis_vals
		}
		p = p + scale_x_continuous(breaks = axis_breaks, labels = axis_vals)
	}

	if (centered && prop_range_y && !x$logodds && !x$probit){ #don't draw this axis for logodds since it makes no sense
		at = seq(min(ice_curves), max(ice_curves), length.out = 5)
		#we need to organize it so it's at zero
		at = at - min(abs(at))

		#check prop type.
		if(prop_type == "range"){
			labels = round(at / x$range_y, 2)  #as a fraction of range of y
		}else{
			labels = round(at / x$sd_y, 2)     #as a fraction of sd(y)
		}
		p = p + scale_y_continuous(
			limits = ylim,
			sec.axis = sec_axis(~., breaks = at, labels = labels)
		)
	} else {
		p = p + scale_y_continuous(limits = ylim)
	}

	p = p + labs(x = xlab, y = ylab) + plot_theme
	print(p)

	invisible(list(ice_curves = ice_curves, range_ice_curves = range_ice_curves, plot_points_indices = plot_points_indices, legend_text = legend_text, pdp = pdp, plot = p))
}
