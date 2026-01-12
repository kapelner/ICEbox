#' Create a plot of a \code{dice} object.
#'
#' Plotting of \code{dice} objects.
#'
#' @param x Object of class \code{dice} to plot.
#' @param plot_margin Extra margin to pass to \code{ylim} as a fraction of the range of \code{x$d_ice_curves}.
#' @param frac_to_plot If \code{frac_to_plot} is less than 1, randomly plot \code{frac_to_plot} fraction of the
#'   curves in \code{x$d_ice_curves}.
#' @param plot_sd If \code{TRUE}, plot the cross-observation sd of partial derivatives below the derivative plots.
#' @param plot_orig_pts_deriv If \code{TRUE}, marks each curve at the location of the derivative estimate at the
#'   location of \code{predictor} actually occurring in the data. If \code{FALSE}
#'   no mark is drawn.
#' @param pts_preds_size Size of points to make if \code{plot_orig_pts_deriv} is \code{TRUE}.
#' @param colorvec Optional vector of colors to use for each curve.
#' @param color_by Optional variable name (or column number) in \code{Xice} to color curves by. If the \code{color_by}
#'   variable has 10 or fewer unique values, a discrete set of colors is used for each value and a legend is
#'   printed and returned. If there are more values, curves are colored from light to dark corresponding
#'   to low to high values of the variable specified by \code{color_by}.
#' @param x_quantile If \code{TRUE}, the plot is drawn with the x-axis taken to be \code{quantile(gridpts)}. If \code{FALSE},
#'   the predictor's original scale is used.
#' @param plot_dpdp If \code{TRUE}, the estimated derivative of the PDP is plotted and highlighted in yellow.
#' @param rug_quantile If not null, tick marks are drawn on the x-axis corresponding to the vector of quantiles specified by this parameter.
#'   Forced to \code{NULL} when \code{x_quantile} is set to \code{TRUE}.
#' @param verbose If \code{TRUE}, prints the color legend to the console.
#' @param ... Additional plotting arguments.
#'
#' @return A list with the following elements.
#'   \item{plot_points_indices}{Row numbers of \code{Xice} of those observations presented in the plot.}
#'   \item{legend_text}{If the \code{color_by} argument was used,
#'   a legend describing the map between the \code{color_by} predictor
#'   and curve colors.}
#'   \item{plot}{The ggplot object used for plotting.}
#'
#' @seealso \code{\link{dice}}
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
#' bhd.ice = ice(object = bhd_rf_mod, X = X, y = y, predictor = "age", frac_to_build = .1)
#'
#' # estimate derivatives, then plot.
#' bhd.dice = dice(bhd.ice)
#' plot(bhd.dice)
#' }
#' @export
plot.dice = function(x, plot_margin = 0.05, frac_to_plot = 1, plot_sd = TRUE, plot_orig_pts_deriv = TRUE,
						pts_preds_size = 1.5, colorvec, color_by = NULL, x_quantile = TRUE, plot_dpdp = TRUE,
						rug_quantile = seq(from = 0, to = 1, by = 0.1), verbose = TRUE, ...){
	#think of x as 'dice_obj'

	DEFAULT_COLORVEC = c("firebrick3", "dodgerblue3", "gold1", "darkorchid4", "orange4", "forestgreen", "grey", "black", "deeppink", "turquoise4")

	arg_list = list(...)

	#some argument checking
	if (!inherits(x, "dice")){
		stop("object is not of class 'dice'")
	}
	if (frac_to_plot <= 0 || frac_to_plot > 1 ){
		stop("frac_to_plot must be in (0,1]")
	}
	if(!is.null(arg_list$ylim) && plot_sd == TRUE){
		stop("Cannot specify both ylim and plot_sd=TRUE.")
	}

	#extract the grid and lines to plot
	grid = x$gridpts
	n_grid = length(grid)
	ecdf_fcn = NULL
	if (x_quantile){
		ecdf_fcn = ecdf(grid)
		grid = ecdf_fcn(grid)
	}
	d_ice_curves = x$d_ice_curves
	N = nrow(d_ice_curves)

	#### figure out the colorvec.
	legend_text = NULL #default is no legend.
	
	# Case 1: Neither colorvec nor color_by provided
	if (missing(colorvec) && missing(color_by)){
		colorvec = sort(rgb(rep(0.4, N), rep(0.4, N), rep(0.4, N), runif(N, 0.4, 0.8)))
	} else if (!missing(color_by)) {
		# Case 2: color_by provided (palette mode if colorvec is provided)
		
		#argument checking first:
		arg_type = class(color_by)
		if(!(arg_type %in% c("character", "numeric"))){
			stop("color_by must be a column name in X or a column index")
		}
		if(inherits(color_by, "character")){
			if(!(color_by %in% names(x$Xice))){
				stop("The predictor name given by color_by was not found in the X matrix")
			}
		} else{  #check numeric
			if( color_by < 1 || color_by > ncol(x$Xice) || (color_by%%1 !=0)){
				stop("color_by must be a column name in X or a column index")
			}
		}
		x_color_by = x$Xice[[color_by]]
		x_unique = unique(x_color_by)
		num_x_color_by = length(x_unique)

		#if there are 10 or fewer unique values of this x value, we use the
		#same color in DEFAULT_COLORVEC for each. Otherwise, we use a rainbow.
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
			x_column_name = ifelse(is.character(color_by), color_by, paste("x_", color_by, sep = ""))
			names(legend_text) = c(x_column_name,"color")
			if (verbose){
				cat("dICE Color Legend\n")
				print(legend_text)
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
				if (verbose) cat("dICE Color Legend: red = low values of the color_by variable and green = high values\n")
			}
		}
	} else {
		# Case 3: Only colorvec provided
		if (length(colorvec) < N) {
			stop("color vector has length ", length(colorvec), " but there are ", N, " lines to plot")
		}
	}


	#pull out a fraction of the lines to plot
	plot_points_indices = which(as.logical(rbinom(N, 1, frac_to_plot)))
	d_ice_curves = d_ice_curves[plot_points_indices, ]
	if (nrow(d_ice_curves) == 0){
		stop("no rows selected: frac_to_plot too small.")
	}
	colorvec = colorvec[plot_points_indices]

	##### now start plotting
	min_dice = min(d_ice_curves)
	max_dice = max(d_ice_curves)
	range_dice = max_dice - min_dice
	min_dice = min_dice - plot_margin * range_dice
	max_dice = max_dice + plot_margin * range_dice

	#get the xlabel if it wasn't already passed explicitly.
	if( is.null(arg_list$xlab)){
		xlab = x$xlab
	} else {
		xlab = arg_list$xlab
	}
	if (x_quantile){
		xlab = paste("quantile(", xlab, ")", sep = "")
	}
	if (!missing(color_by)){
		xlab = paste(xlab, "colored by", color_by)
	}

	#same for y label
	if( is.null(arg_list$ylab)){
		if (x$logodds){
			ylab = "partial log-odds"
		} else {
			ylab = "derivative f-hat"
		}
	} else {
		ylab = arg_list$ylab
	}

	#set ylim if not passed explicitly
	if( is.null(arg_list$ylim) ){
		if(plot_sd){
			offset = 1.5 * max(x$sd_deriv)
			ylim = c(min_dice - offset, max_dice)
		}else{
			ylim = c(min_dice, max_dice)
		}
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

	dice_df = data.table(
		curve_id = rep(seq_len(nrow(d_ice_curves)), each = n_grid),
		x = rep(grid, times = nrow(d_ice_curves)),
		y = as.vector(t(d_ice_curves)),
		color = rep(colorvec, each = n_grid)
	)

	lwd_scale = 0.5
	p = ggplot(dice_df, aes(x = x, y = y, group = curve_id, color = color)) +
		geom_line() +
		scale_color_identity(guide = "none")

	if (plot_orig_pts_deriv){ #indicate the fitted values associated with observed xj values
		deriv_actual = x$actual_deriv[plot_points_indices]

		if (x_quantile){
			xj = ecdf_fcn(x$xj)[plot_points_indices]
		} else {
			xj = x$xj[plot_points_indices]
		}
		points_df = data.table(x = xj, y = deriv_actual, color = colorvec)
		p = p +
			geom_point(
				data = points_df,
				aes(x = x, y = y),
				color = "black",
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


	#if plot_dpdp is true, plot actual dpdp (in the sense of Friedman '01)
	if (plot_dpdp){
		friedman_dpdp = x$dpdp

		#calculate the line thickness based on how many lines there are
		num_lines = length(plot_points_indices)
		#every 100 lines we get 0.5 more highlight up to 8
		dpdp_df = data.table(x = grid, y = friedman_dpdp)
		p = p +
			geom_line(
				data = dpdp_df,
				aes(x = x, y = y),
				color = "yellow",
				linewidth = lwd_scale * min(5.5 + (num_lines / 100) * 0.75, 8),
				inherit.aes = FALSE
			) +
			geom_line(
				data = dpdp_df,
				aes(x = x, y = y),
				color = "BLACK",
				linewidth = lwd_scale * 4,
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

	if (plot_sd){
		abline_y = ylim[1] + offset
		sd_df = data.table(x = grid, y = x$sd_deriv + ylim[1])
		p = p +
			geom_hline(yintercept = abline_y, color = rgb(0.8,0.8,0.8)) +
			geom_line(
				data = sd_df,
				aes(x = x, y = y),
				color = "black",
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

	if (plot_sd){
		at = seq(ylim[1], ylim[1] + max(x$sd_deriv), length.out = 2)
		labels = round(seq(0, max(x$sd_deriv), length.out = 2), 1)
		p = p + scale_y_continuous(
			limits = ylim,
			sec.axis = sec_axis(~., breaks = at, labels = labels, name = "sd(deriv)")
		)
	} else {
		p = p + scale_y_continuous(limits = ylim)
	}

	p = p + labs(x = xlab, y = ylab) + plot_theme
	print(p)

	invisible(list(plot_points_indices = plot_points_indices, legend_text = legend_text, plot = p))
}
