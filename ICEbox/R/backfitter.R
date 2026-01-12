#' Backfitting for Additive Models
#'
#' Fits a model of the form \eqn{\hat{f}(x) = \hat{g}_{1}(x_S) + \hat{g}_{2}(x_C)} using backfitting.
#'
#' @param X The design matrix.
#' @param y The response vector.
#' @param predictor The name or index of the predictor of interest (\eqn{x_S}).
#' @param fitMethod A function that accepts \code{X} and \code{y} and returns a fitted model.
#' @param predictfcn A function that accepts \code{object} and \code{newdata} and returns predictions.
#' @param eps Convergence threshold.
#' @param iter.max Maximum number of iterations.
#' @param verbose If \code{TRUE}, prints progress messages.
#' @param ... Additional arguments passed to \code{fitMethod}.
#'
#' @return An object of class \code{backfitter}.
#'
#' @export
backfitter = function(X, y, predictor, fitMethod, predictfcn, eps = 0.01, iter.max = 10, verbose = TRUE, ...){

	assert_data_frame(X)
	assert_numeric(y, len = nrow(X))
	if (is.character(predictor)){
		assert_character(predictor, len = 1)
		assert_subset(predictor, names(X))
	} else {
		assert_integerish(predictor, lower = 1, upper = ncol(X))
	}
	assert_function(fitMethod)
	assert_function(predictfcn)
	assert_number(eps, lower = 0)
	assert_count(iter.max, positive = TRUE)
	assert_flag(verbose)

	# order by the predictor
	pred_vals = if (is.character(predictor)) X[[predictor]] else X[, predictor]
	xorder = order(pred_vals)
	X = X[xorder, ]
	y = y[xorder]

	if (!is.numeric(predictor)){
		which_col = which(names(X) == predictor)
		Xc = X[, -which_col]
	} else {
		Xc = X[, -predictor]
	}
	Xs = if (is.character(predictor)) X[[predictor]] else X[, predictor]

	# initialize
	N = nrow(X)
	g1_of_Xs = rep(0, N)
	g2_of_Xc = rep(0, N)
	current_g2 = NULL

	# supsmu will condense ties in x, so length(new_g1) = length(unique(Xs))
	times_to_repeat = as.numeric(table(Xs)) # remember Xs is sorted

	OneStep = function(){
		# do g2 first
		new_g2_mod = fitMethod(X = Xc, y = (y - g1_of_Xs), ...)
		new_g2 = predictfcn(object = new_g2_mod, Xc)
		new_g1 = supsmu(x = Xs, y = (y - new_g2))$y
		new_g1 = rep(new_g1, times_to_repeat) # matches length of new_g2 now.
		return(list(new_g1 = new_g1, new_g2 = new_g2, new_g2_mod = new_g2_mod))
	}

	delta = Inf
	iter = 0
	while(delta > eps && iter < iter.max){
		# one iteration
		nextStep = OneStep()
		
		# compute delta
		delta = sum((nextStep$new_g1 - g1_of_Xs)^2) / sum(g1_of_Xs^2)
		delta = delta + sum((nextStep$new_g2 - g2_of_Xc)^2) / sum(g2_of_Xc^2)

		# update
		current_g2 = nextStep$new_g2_mod
		g1_of_Xs = nextStep$new_g1
		g2_of_Xc = nextStep$new_g2
		iter = iter + 1
		
		# print message
		if(verbose){
			cat(paste("iter", iter, " delta: ", round(delta, 3), sep = ""))
			cat("\n")
		}
	}
	
	bf_obj = list(
		g1_of_Xs = g1_of_Xs,
		g2_of_Xc = g2_of_Xc,
		g2_mod = current_g2,
		X = X,
		y = y,
		predictor = predictor,
		iter = iter,
		delta = delta,
		fitMethod = fitMethod,
		predictfcn = predictfcn
	)
	class(bf_obj) = "backfitter"
	return(bf_obj)
}
