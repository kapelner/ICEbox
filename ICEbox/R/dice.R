#' Creates an object of class \code{dice}.
#'
#' Estimates the partial derivative function for each curve in an \code{ice} object.
#' See Goldstein et al (2013) for further details.
#'
#' @param ice_obj Object of class \code{ice}. This function generates partial derivative
#'   estimates for each row in \code{ice_obj$ice_curves}.
#' @param DerivEstimator Optional function with a single argument \code{y}. Returns the estimated
#'   partial derivative of a function sampled at the points (\code{ice_obj$gridpts},\code{y}).
#'   If NULL, the default uses a Savitzky-Golay filter to estimate the first derivative.
#' @param use_supsmu If \code{TRUE}, uses the old \code{supsmu} based derivative estimation logic.
#'   This is much slower than the default Savitzky-Golay filter.
#' @param verbose If \code{TRUE}, prints messages about the procedure's progress.
#' @param num_cores Integer number of cores to use for parallel derivative estimation. Defaults to 1.
#' @param sg_poly_order Polynomial order for Savitzky-Golay filter. Default is 2.
#' @param sg_window_size Window size for Savitzky-Golay filter. Default is 30\% of the grid.
#'
#' @return A list of class \code{dice} with the following elements. Most are passed directly through
#' from \code{ice_object} and exist to enable various plotting facilities.
#'
#'   \item{d_ice_curves}{Matrix of dimension \code{nrow(Xice)} by \code{length(gridpts)}.
#'   Each row corresponds to an observation's d-ICE curve, estimated at the values of \code{predictor} in \code{gridpts}.}
#'   \item{xj}{The actual values of \code{predictor} observed in the data in the order
#'   of \code{Xice}.}
#'   \item{actual_deriv}{Vector of length \code{nrow(Xice)} containing the estimated partial derivatives
#'   at the value of the \code{predictor} actually found in \code{Xice}.}
#'   \item{sd_deriv}{Vector of length \code{length(gridpts)} with the cross-observation sd of partial derivative
#'   estimates. For instance \code{sd_deriv[1]} equals \code{sd(d_ice_curves[,1])}.}
#'   \item{logodds}{Passed from \code{ice_object}. If \code{TRUE}, \code{d_ice_curves} are
#'   estimated derivatives of the centered log-odds.}
#'   \item{gridpts}{Passed from \code{ice_object}.}
#'   \item{predictor}{Passed from \code{ice_object}.}
#'   \item{xlab}{Passed from \code{ice_object}.}
#'   \item{nominal_axis}{Passed from \code{ice_object}.}
#'   \item{range_y}{Passed from \code{ice_object}.}
#'   \item{Xice}{Passed from \code{ice_object}.}
#'   \item{dpdp}{The estimated partial derivative of the PDP.}
#'
#' @references
#' Goldstein, A., Kapelner, A., Bleich, J., and Pitkin, E., Peeking
#'   Inside the Black Box: Visualizing Statistical Learning With Plots of
#'   Individual Conditional Expectation. (2014) Journal of Computational
#'   and Graphical Statistics, in press \cr
#'
#' @seealso \code{\link{ice}}, \code{\link{dice}}
#' @examples
#' \dontrun{
#' # same examples as for 'ice', but now create a derivative estimate as well.
#' require(ICEbox)
#' require(randomForest)
#' require(MASS) #has Boston Housing data, Pima
#'
#' ########  regression example
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
#' # make a dice object:
#' bhd.dice = dice(bhd.ice)
#'
#' #### classification example
#' data(Pima.te)  #Pima Indians diabetes classification
#' y = Pima.te$type
#' X = Pima.te
#' X$type = NULL
#'
#' ## build a RF:
#' pima_rf = randomForest(x = X, y = y)
#'
#' ## Create an 'ice' object for the predictor "skin":
#' # For classification we plot the centered log-odds. If we pass a predict
#' # function that returns fitted probabilities, setting logodds = TRUE instructs
#' # the function to set each ice curve to the centered log-odds of the fitted
#' # probability.
#' pima.ice = ice(object = pima_rf, X = X, predictor = "skin", logodds = TRUE,
#'                     predictfcn = function(object, newdata){
#'                          predict(object, newdata, type = "prob")[, 2]
#'                     }
#'               )
#'
#' # make a dice object:
#' pima.dice = dice(pima.ice)
#' }
#' @export
dice = function(ice_obj, DerivEstimator = NULL, use_supsmu = FALSE, verbose = TRUE, num_cores = 1, sg_poly_order = 2, sg_window_size = NULL){

	assert_class(ice_obj, "ice")
	if (!is.null(DerivEstimator)){
		assert_function(DerivEstimator)
	}
	assert_flag(use_supsmu)
	assert_flag(verbose)
	assert_integerish(num_cores, len = 1, lower = 1)
	assert_integerish(sg_poly_order, len = 1, lower = 1)
	if (!is.null(sg_window_size)){
		assert_integerish(sg_window_size, len = 1, lower = 3)
	}

	#error checking:
	if(!inherits(ice_obj, "ice")){
		stop("ice_obj is not a valid ice object.")
	}

	gridpts = ice_obj$gridpts
	dice_obj = ice_obj

	if (is.null(DerivEstimator)){
		if (use_supsmu){
			if (verbose) cat("Estimating derivatives using supsmu + centered differences\n")
			EstimatorWrapper = function(y){
				as.vector(derivative_cpp(matrix(supsmu(x = gridpts, y = y)$y, nrow = 1), gridpts))
			}
		} else {
			# Use Savitzky-Golay filter (Rcpp implementation)
			if (is.null(sg_window_size)){
				n_grid = length(gridpts)
				# heuristic: 30% of grid, min 5, ensure odd
				sg_window_size = max(5, round(n_grid * 0.3))
				if (sg_window_size %% 2 == 0) sg_window_size = sg_window_size + 1
			}
			
			if (verbose) cat("Estimating derivatives using Savitzky-Golay filter (window =", sg_window_size, ", order =", sg_poly_order, ")\n")

			# (1) Smooth the curves using SG (deriv=0)
			# (2) Compute the derivative using centered differences on the actual grid (derivative_cpp)
			# This correctly handles non-uniform grids, mirroring the old supsmu + D1tr logic.
			# Polynomial order 2 is generally more stable for derivative estimation.
			smoothed_curves = sg_smooth_cpp(ice_obj$ice_curves, window_size = sg_window_size, order = sg_poly_order, deriv = 0, n_cores = num_cores)
			dice_obj$d_ice_curves = derivative_cpp(smoothed_curves, gridpts, n_cores = num_cores)
			
			# do it for the pdp as well
			smoothed_pdp = as.vector(sg_smooth_cpp(matrix(ice_obj$pdp, nrow = 1), window_size = sg_window_size, order = sg_poly_order, deriv = 0, n_cores = num_cores))
			dice_obj$dpdp = as.vector(derivative_cpp(matrix(smoothed_pdp, nrow = 1), gridpts, n_cores = num_cores))
			dice_obj$ice_curves = NULL
		}
	} 
	
	if (!is.null(DerivEstimator) || (is.null(DerivEstimator) && use_supsmu)) {
		if (is.null(DerivEstimator) && use_supsmu){
			# EstimatorWrapper already defined above
		} else {
			#argument checking???
			EstimatorWrapper = function(y){
				DerivEstimator(y=y,x=gridpts)		
			}
		}
		
		#compute derivatives
		ice_curves = ice_obj$ice_curves
		n_curves = nrow(ice_curves)

		if (num_cores > 1 && n_curves > 1){
			num_workers = min(num_cores, n_curves)
			estimate_row = function(row_idx){
				EstimatorWrapper(ice_curves[row_idx, ])
			}
			if (.Platform$OS.type == "windows"){
				cl = parallel::makeCluster(num_workers)
				on.exit(parallel::stopCluster(cl), add = TRUE)
				loaded_pkgs = loadedNamespaces()
				parallel::clusterExport(cl, varlist = c("loaded_pkgs", "ice_curves", "EstimatorWrapper"), envir = environment())
				parallel::clusterEvalQ(cl, {
					invisible(lapply(loaded_pkgs, require, character.only = TRUE))
				})
				deriv_list = parallel::parLapply(cl, seq_len(n_curves), estimate_row)
			} else {
				deriv_list = parallel::mclapply(seq_len(n_curves), estimate_row, mc.cores = num_workers)
			}
			dice_obj$d_ice_curves = do.call(rbind, deriv_list)
		} else {
			dice_obj$d_ice_curves = t(apply(ice_obj$ice_curves, 1, FUN = EstimatorWrapper))
		}
		dice_obj$ice_curves = NULL

		#do it for the pdp as well.
		dice_obj$dpdp = EstimatorWrapper(ice_obj$pdp)
	}


	#figure out point on each curve that corresponds to observed X.
    col_idx_of_actual = c(1, 1 + cumsum(diff(ice_obj$xj)>0))
	row_idx = 1:nrow(dice_obj$d_ice_curves)
	actual_deriv_idx = cbind(row_idx, col_idx_of_actual)
	dice_obj$actual_deriv = dice_obj$d_ice_curves[actual_deriv_idx]

	#compute the sd of the derivatives at each gridpt.
	dice_obj$sd_deriv = colSds_cpp(dice_obj$d_ice_curves, n_cores = num_cores)
	
	#clean up, make it of class 'dice'	
	dice_obj$actual_prediction = NULL
	dice_obj$predictfcn = NULL
	dice_obj$frac_to_build = NULL
	dice_obj$indices_to_build = NULL
	class(dice_obj) = "dice"
	
	invisible(dice_obj)
}
