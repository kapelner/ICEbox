#' Lineup plot for additivity
#'
#' This function creates a lineup plot to assess the additivity of a predictor's effect.
#' It uses a nonparametric bootstrap approach to generate null plots.
#'
#' @param backfit_obj An object of class \code{backfitter}.
#' @param fitMethod A function that accepts \code{X} and \code{y} and returns a fitted model.
#' @param realICE The \code{ice} object for the real data.
#' @param figs The total number of plots in the lineup (including the real one). Default is 10.
#' @param colorvecfcn Optional function to generate a color vector for the curves.
#' @param usecolorvecfcn_inreal If \code{TRUE}, use \code{colorvecfcn} for the real plot.
#' @param null_predictfcn Optional prediction function for the null models.
#' @param ... Additional arguments passed to \code{plot.ice}.
#'
#' @return An object of class \code{additivityLineup} (invisibly).
#'
#' @export
additivityLineup = function(backfit_obj, fitMethod, realICE, figs = 10, colorvecfcn, usecolorvecfcn_inreal = FALSE, null_predictfcn, ...){

	assert_class(backfit_obj, "backfitter")
	assert_function(fitMethod)
	assert_class(realICE, "ice")
	assert_count(figs, positive = TRUE)
	if (!missing(colorvecfcn)) assert_function(colorvecfcn)
	assert_flag(usecolorvecfcn_inreal)

	# predictfcn for nulls:
	if(missing(null_predictfcn)){
		if(is.null(realICE$predictfcn)){
			null_predictfcn = NULL
		} else {
			null_predictfcn = realICE$predictfcn	
		}
	}

	predictor = backfit_obj$predictor
	arg_list = list(...)
  
	# frac_to_build is not allowed
	if(!is.null(arg_list$frac_to_build)){
		cat("Cannot specify frac_to_build. Can specify frac_to_plot, which applies to the realICE, and frac_to_build for the null ICEs is then inferred to plot the same number of curves.\n")
	}
  
	centered = if (is.null(arg_list$centered)) FALSE else arg_list$centered
	centered_percentile = arg_list$centered_percentile
	if(is.null(centered_percentile) && centered == TRUE){ 
		centered_percentile = 0.01  # default in plot.ice
	}

	# and some more for frac_to_plot
	frac_to_build_null = 1

	if(!is.null(arg_list$frac_to_plot)){
		frac_to_plot = arg_list$frac_to_plot
		warning_msg = paste("'frac_to_plot' only applies to plotting 'realICE'.",
		"'frac_to_build' is set in null ICEs to ensure the same number of curves are plotted for null and real plots.", sep = "\n")  
		warning(warning_msg)
		frac_to_build_null = nrow(realICE$ice_curves) * arg_list$frac_to_plot / nrow(backfit_obj$X)

		# fix indices to plot so that the ylim's can be constrained to only those
		# curves actually plotted.
		plot_points_indices = which(as.logical(rbinom(nrow(realICE$ice_curves), 1, frac_to_plot)))
		realICE$ice_curves = realICE$ice_curves[plot_points_indices, ]
		realICE$gridpts = realICE$gridpts[plot_points_indices]
		realICE$Xice = realICE$Xice[plot_points_indices, ]
		realICE$xj = realICE$xj[plot_points_indices]
		frac_to_plot = 1
	}

	# figure out min and max of real ice object -- depends on centering  
	if(centered){
		centering_vector = realICE$ice_curves[, max(1, ceiling(ncol(realICE$ice_curves) * centered_percentile))]
		rg = range(realICE$ice_curves - centering_vector) 
	} else {
		rg = range(realICE$ice_curves)
	}
	icecurve_min = rg[1]
	icecurve_max = rg[2]
  
	additive_fit = backfit_obj$g1_of_Xs + backfit_obj$g2_of_Xc
	additive_res = backfit_obj$y - additive_fit
	
	null_additive_fits = list()
	null_ices = list()
  
	for(i in 1:(figs - 1)){
		response = additive_fit + sample(additive_res, size = length(additive_res), replace = FALSE)
		new_fit = fitMethod(X = backfit_obj$X, y = response)
		null_additive_fits[[i]] = new_fit

		if(is.null(null_predictfcn)){ # no predictfcn found, use generic
			null_ices[[i]] = ice(new_fit, X = backfit_obj$X, predictor = predictor, y = backfit_obj$y,
								frac_to_build = frac_to_build_null, verbose = FALSE)
		} else {
			null_ices[[i]] = ice(new_fit, X = backfit_obj$X, predictor = predictor, y = backfit_obj$y, 
								frac_to_build = frac_to_build_null, predictfcn = null_predictfcn, verbose = FALSE)
		}
    
		### keep track of min and max. 
		if(centered){
			centering_vector = null_ices[[i]]$ice_curves[, max(1, ceiling(ncol(null_ices[[i]]$ice_curves) * centered_percentile))]
			rg = range(null_ices[[i]]$ice_curves - centering_vector) # range for centered plot
		} else {
			rg = range(null_ices[[i]]$ice_curves)      
		}

		# update min range and max range
		icecurve_min = min(icecurve_min, rg[1])
		icecurve_max = max(icecurve_max, rg[2])
		cat("Finished null ice ", i, "\n")
	} # end loop through null ices
  
	# graphics
	num_plot_cols = ceiling(figs / 4)
	num_plot_rows = ceiling(figs / num_plot_cols)
    
    # Note: Since plot.ice uses ggplot2, sequential printing in a lineup
    # would ideally use something like gridExtra::grid.arrange.
    # For now we collect the plots and return them.
    plot_list = list()
    
	ylim = c(icecurve_min, icecurve_max)
  
	# argument list for the null plots.
	null_arg_list = arg_list
	null_arg_list$ylim = ylim
	null_arg_list$frac_to_plot = 1
  
	# randomly place the real plot somewhere...
	where_to_place = sample(1:figs, 1)
	plotted_truth = FALSE
	for(i in 1:figs){

		if(plotted_truth){
			idx_to_plot = i - 1
		} else {
			idx_to_plot = i
		}

		if((!plotted_truth) && (i == where_to_place)){
			if(usecolorvecfcn_inreal){
				colors = colorvecfcn(realICE)	
				p_obj = plot(realICE, ylim = ylim, colorvec = colors, ...)
			} else {
				p_obj = plot(realICE, ylim = ylim, ...)
			}
			plotted_truth = TRUE
		} else {
			null_arg_list$x = null_ices[[idx_to_plot]]
			if(!missing(colorvecfcn)){
				colors = colorvecfcn(null_ices[[idx_to_plot]])
				null_arg_list$colorvec = colors
			}
			p_obj = do.call(plot, null_arg_list) 
		}
        plot_list[[i]] = p_obj$plot
	}
  
	al_obj = list(location = where_to_place, null_additive_fits = null_additive_fits, 
					null_ices = null_ices, frac_to_build_null = frac_to_build_null, plots = plot_list)
	class(al_obj) = "additivityLineup"  
	invisible(al_obj)
}
