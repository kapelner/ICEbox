#' Creates an object of class \code{ice}.
#'
#' Creates an \code{ice} object with individual conditional expectation curves
#' for the passed model object, \code{X} matrix, predictor, and response. See
#' Goldstein et al (2013) for further details.
#'
#' @param object The fitted model to estimate ICE curves for.
#' @param X The design matrix we wish to estimate ICE curves for. Rows are observations, columns are
#'   predictors. Typically this is taken to be \code{object}'s training data, but this is not
#'   strictly necessary.
#' @param y Optional vector of the response values \code{object} was trained on. It is used
#'   to compute y-axis ranges that are useful for plotting. If not passed, the range
#'   of predicted values is used and a warning is printed.
#' @param predictor The column number or variable name in \code{X} of the predictor of interest,
#'   (\eqn{x_S = X[, j]}{x_S= X[, j]}).
#' @param predictfcn Optional function that accepts two arguments, \code{object} and \code{newdata}, and
#'   returns an \code{N} vector of \code{object}'s predicted response for data \code{newdata}.
#'   If this argument is not passed, the procedure attempts to find a generic \code{predict}
#'   function corresponding to \code{class(object)}.
#' @param verbose If \code{TRUE}, prints messages about the procedure's progress.
#' @param frac_to_build Number between 0 and 1, with 1 as default. For large \code{X} matrices or fitted models
#'   that are slow to make predictions, specifying \code{frac_to_build} less than 1 will choose
#'   a subset of the observations to build curves for. The subset is chosen such that the remaining
#'   observations' values of \code{predictor} are evenly spaced throughout the quantiles of the
#'   full \code{X[,predictor]} vector.
#' @param indices_to_build Vector of indices, \eqn{\subset \{1, \ldots, nrow(X)\}}{each element between \code{1} and \code{nrow(X)}} specifying which observations to build ICE curves for. As this is an alternative to setting \code{frac_to_build}, both
#'   cannot be specified.
#' @param num_grid_pts Optional number of values in the range of \code{predictor} at which to estimate each curve.
#'   If missing, the curves are estimated at each unique value of \code{predictor}
#'   in the \code{X} observations we estimate ICE curves for.
#' @param logodds If \code{TRUE}, for classification creates PDPs by plotting the centered log-odds implied by the fitted
#'   probabilities. We assume that the generic or passed predict function
#'   returns probabilities, and so the flag tells us to transform these to centered logits after
#'   the predictions are generated. Note: \code{probit} cannot be \code{TRUE}.
#' @param probit If \code{TRUE}, for classification creates PDPs by plotting the probit implied by the fitted
#'   probabilities. We assume that the generic or passed predict function
#'   returns probabilities, and so the flag tells us to transform these to probits after
#'   the predictions are generated. Note: \code{logodds} cannot be \code{TRUE}.
#' @param num_cores Integer number of cores to use for parallel prediction. Defaults to 1.
#' @param ... Other arguments to be passed to \code{object}'s generic predict function.
#'
#' @return A list of class \code{ice} with the following elements:
#'   \item{gridpts}{Sorted values of \code{predictor} at which each curve is estimated. Duplicates
#'   are removed -- by definition, elements of \code{gridpts} are unique.}
#'   \item{ice_curves}{Matrix of dimension \code{nrow(X)} by \code{length(gridpts)}.
#'   Each row corresponds to an observation's ICE curve, estimated at the values of \code{predictor} in
#'   \code{gridpts}.}
#'   \item{xj}{The actual values of \code{predictor} observed in the data in the order
#'   of \code{Xice}.}
#'   \item{actual_prediction}{Vector of length \code{nrow(X)} containing the model's
#'   predictions at the actual value of the predictors in the order of \code{Xice}.}
#'   \item{xlab}{String with the predictor name corresponding to \code{predictor}. If \code{predictor}
#'   is a column number, \code{xlab} is set to \code{colnames(X)[, predictor]}.}
#'   \item{nominal_axis}{If \code{TRUE}, \code{length(gridpts)} is 5 or fewer; otherwise \code{FALSE}.
#'   When \code{TRUE} the \code{plot} function treats the x-axis as if x is nominal.}
#'   \item{range_y}{If \code{y} was passed, the range of the response. Otherwise it defaults to be
#'   \code{max(ice_curves)} - \code{min(ice_curves)} and a message is printed to the console.}
#'   \item{sd_y}{If \code{y} was passed, the standard deviation of the response. Otherwise it is defaults to
#'   \code{sd(actual_prediction)} and a message is printed to the console.}
#'   \item{Xice}{A matrix containing the subset of \code{X} for which ICE curves are estimated.
#'   Observations are ordered to be increasing in \code{predictor}. This ordering is the same one
#'   as in \code{ice_curves}, \code{xj} and \code{actual_prediction}, meaning for all these objects
#'   the \code{i}-th element refers to the same observation in \code{X}.}
#'   \item{pdp}{A vector of size \code{length(gridpts)} which is a numerical approximation to the partial
#'   dependence function (PDP) corresponding to the estimated ICE curves. See Goldstein et al (2013) for a discussion
#'   of how the PDP is a form of post-processing. See Friedman (2001) for a description of PDPs.}
#'   \item{predictor}{Same as the argument, see argument description.}
#'   \item{logodds}{Same as the argument, see argument description.}
#'   \item{indices_to_build}{Same as the argument, see argument description.}
#'   \item{frac_to_build}{Same as the argument, see argument description.}
#'   \item{predictfcn}{Same as the argument, see argument description.}
#'
#' @references
#' Jerome Friedman. Greedy Function Approximation: A Gradient Boosting Machine. The Annals of Statistics,
#'   29(5): 1189-1232, 2001.
#'
#' Goldstein, A., Kapelner, A., Bleich, J., and Pitkin, E., Peeking
#'   Inside the Black Box: Visualizing Statistical Learning With Plots of
#'   Individual Conditional Expectation. (2014) Journal of Computational
#'   and Graphical Statistics, in press
#'
#' @seealso \code{\link{plot.ice}}, \code{\link{print.ice}}, \code{\link{summary.ice}}
#' @examples
#' \dontrun{
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
#' }
#' @export
ice = function(object, X, y,
		predictor, predictfcn, 
		verbose = TRUE, frac_to_build = 1, indices_to_build = NULL, 
		num_grid_pts, logodds = FALSE, probit = FALSE, num_cores = 1, ...){

	assert_true(!missing(object))
	assert_true(!missing(X))
	assert_true(is.data.frame(X) || is.matrix(X))
	if (!missing(y)){
		assert_numeric(y, any.missing = TRUE)
	}
	assert_true(!missing(predictor))
	assert_true(is.numeric(predictor) || is.character(predictor))
	if (is.character(predictor)){
		assert_character(predictor, len = 1)
		if (!(predictor %in% colnames(X))) {
			stop("predictor '", predictor, "' not found in colnames(X)")
		}
	} else {
		assert_integerish(predictor, len = 1, lower = 1, upper = ncol(X))
	}
	if (!missing(predictfcn)){
		assert_function(predictfcn)
	}
	assert_flag(verbose)
	assert_number(frac_to_build)
	assert_integerish(indices_to_build, min.len = 1, null.ok = TRUE)
	if (!missing(num_grid_pts)){
		assert_integerish(num_grid_pts, len = 1, lower = 1)
	}
	assert_flag(logodds)
	assert_flag(probit)
	assert_list(list(...))
	assert_integerish(num_cores, len = 1, lower = 1)

	MAX_NUM_UNIQUE_PTS_NOMINAL = 5

	# Convert X to data.table for faster subsetting and assignment
	# We use copy() to avoid modifying the user's object in place if it's already a data.table
	X = as.data.table(copy(X))
	
	# Ensure predictor is a name for data.table operations
	pred_name = if(is.numeric(predictor)) names(X)[predictor] else predictor

	#check for factor
	if (inherits(X[[pred_name]], "factor") || inherits(X[[pred_name]], "character")){
		stop("ICE does not support factor attributes")
	}
	
	if(!is.numeric(frac_to_build) || frac_to_build > 1 || frac_to_build < 0 ){
		stop("frac_to_build must be in (0, 1]")
	}
  
  if(!missing(y) && inherits(y, "factor")){
    stop("Do not pass y when it is categorical variable.")
  }
  
  if (logodds && probit){
	  stop("You must employ either logodds OR probit but not both.")
  }
  

	######## (1) check inputs
	# (a) check for valid prediction routine...
	if (!missing(predictfcn)){
		fcn_args = names(formals(predictfcn))
		if (!("object" %in% fcn_args && "newdata" %in% fcn_args)){
			stop("predictfcn must have 'object' and 'newdata' arguments")
		} else {
			use_generic = FALSE
		}
	} else { 
		#check for default prediction associated with class(object)
		#this is harder than it should be because some classes have class(object)
		#return a vector (randomForest, for instance).
		classvec = class(object) #may have multiple classes
		found_predict = FALSE
		i = 1
		for (i in 1 : length(classvec)){
			if (length(grep("predict", methods(class = classvec[i]))) > 0){
				found_predict = TRUE
				break
			}
		}
		if (!found_predict){
			stop("No generic predict method found for this object.")
		} else {
			use_generic = TRUE
		}
	}

	if (missing(predictfcn)){
		predictfcn = NULL
	}

	
		######## (2)
		N = nrow(X)
		# grid points
		
			# Convert X to data.table for faster subsetting and assignment. 
			# We force a copy to prevent side effects on the user's data.
			X = data.table::copy(data.table::as.data.table(X))		
		# Ensure predictor is a name for data.table operations
		pred_name = if(is.numeric(predictor)) names(X)[predictor] else predictor
	
		#now create xj-to-predict values
		xj = X[[pred_name]]  #fix so predictor can be given as a name. 
		grid_pts = sort(X[[pred_name]])
	
		X[, .original_indices := .I]
		# there are 3 cases: frac_to_build specificed, indices specified, or nothing specified.
	if (frac_to_build < 1 && !missing(indices_to_build)){
		stop("\"frac_to_build\" and \"indices_to_build\" cannot both be specified simultaneously")
	}

	# 1: check fraction to build
	if (frac_to_build < 1){
		# we don't sample randomly -- we ensure uniform sampling across
		# quantiles of xj so as to not leave out portions of the dist'n of x.
		setorderv(X, pred_name)
		nskip = round(1 / frac_to_build)
		xj_sorted_indices_to_build = seq(1, N, by = nskip)
		X = X[xj_sorted_indices_to_build]
		xj = X[[pred_name]]
		grid_pts = sort(xj)
		indices_to_build = X$.original_indices
	} else { 
		
	  #2: indices specified:
	  if (!missing(indices_to_build)){
	  	#extract the indicies the user asks for first, THEN order the remaining
	  	#X matrix by the xj's left in the sub-sample defined by the user.
		X = X[indices_to_build]
		setorderv(X, pred_name)
		xj = X[[pred_name]]
		grid_pts = sort(xj)
	  } else { #3: nothing specified, so just re-order by xj
		setorderv(X, pred_name)
		xj = X[[pred_name]]
	  }	
	}
	X[, .original_indices := NULL]
	grid_pts = unique(grid_pts)
	num_unique_pts = length(grid_pts)
	
	
	#now handle less grid pts
	if (!missing(num_grid_pts)){
		if (num_grid_pts > num_unique_pts){
			warning(paste("the number of grid points specified,", num_grid_pts, "is larger than the number of unique values,", num_unique_pts, "(defaulting to number of unique values)"))	
		} else {
			grid_pts = grid_pts[round(seq(from = 1, to = num_unique_pts, length.out = num_grid_pts))]	
		}		
	}
	
	# generate partials
	predict_args = list(...)
	if (use_generic){
		actual_prediction = do.call("predict", c(list(object, X), predict_args))
	} else {
		actual_prediction = predictfcn(object = object, newdata = X)
	}
	
	if (inherits(actual_prediction, "factor")){
		stop("The predict function must return probabilities (not levels of a factor).")
	}
	if (logodds || probit){ 	
		min_pred = min(actual_prediction)
		max_pred = max(actual_prediction)
		
		#do some p_hat \in [0, 1] checking
		if (min_pred < 0){ 
			stop("the logodds option is on but predict returns values less than 0 (these should be probabilities!)")
		} else if (max_pred > 1){
			stop("the logodds option is on but predict returns values greater than 1 (these should be probabilities!)")
		}
		if (min_pred == 0){
			second_lowest = min(actual_prediction[actual_prediction > 0])
			if (is.na(second_lowest)){ 
				second_lowest = .0001
			}
			lowest_practical_prob = mean(c(second_lowest, 0))
			actual_prediction[actual_prediction == 0] = lowest_practical_prob
			#warning("At least one probability was predicted to be 0, ICEbox is using ", lowest_practical_prob, " for the value(s) instead.")
		}
		if (max_pred == 1){
			second_highest = max(actual_prediction[actual_prediction < 1])
			if (is.na(second_highest)){ 
				second_highest = .9999
			}
			highest_practical_prob = mean(c(second_highest, 1)) #arbitrarily, halfway between 1 and second_highest
			actual_prediction[actual_prediction == 1] = highest_practical_prob
			#warning("At least one probability was predicted to be 1, ICEbox is using ", highest_practical_prob, " for the value(s) instead.")
		}
		
		if (logodds){
			# centered logit formula: log(p) - 0.5 * (log(p) + log(1-p)) = 0.5 * log(p/(1-p))
			actual_prediction = 0.5 * qlogis(actual_prediction)
		} else if (probit){
			actual_prediction = qnorm(actual_prediction)
		}
		
	}
	
	ice_curves = matrix(NA, nrow = nrow(X), ncol = length(grid_pts))
	colnames(ice_curves) = grid_pts

	#Compute actual pdp. Note that this is averaged over the observations
	#we sample, so this might be different from the 'true' pdp if frac_to_build < 0.
	xvec_temp = X[[pred_name]]  #actual setting of this column after ordering.
	n_obs = nrow(X)
	grid_len = length(grid_pts)
	max_rows_per_predict = getOption("ICEbox.max_rows_per_predict", 2e5)
	chunk_size = max(1L, floor(max_rows_per_predict / n_obs))
	chunk_size = min(chunk_size, grid_len)
	chunk_starts = seq(1L, grid_len, by = chunk_size)

	predict_chunk = function(start_idx){
		end_idx = min(grid_len, start_idx + chunk_size - 1L)
		chunk_grid = grid_pts[start_idx : end_idx]
		X_chunk = X[rep(seq_len(n_obs), times = length(chunk_grid))]
		set(X_chunk, j = pred_name, value = rep(chunk_grid, each = n_obs))
		if (use_generic){
			preds = do.call("predict", c(list(object, X_chunk), predict_args))
		} else {
			preds = predictfcn(object = object, newdata = X_chunk)
		}
		list(start_idx = start_idx, end_idx = end_idx, preds = preds)
	}

	if (num_cores > 1 && length(chunk_starts) > 1){
		num_workers = min(num_cores, length(chunk_starts))
		if (.Platform$OS.type == "windows"){
			cl = parallel::makeCluster(num_workers)
			on.exit(parallel::stopCluster(cl), add = TRUE)
			loaded_pkgs = loadedNamespaces()
			parallel::clusterExport(cl, varlist = c("loaded_pkgs", "X", "pred_name", "use_generic", "object", "predictfcn",
				"predict_args", "grid_pts", "grid_len", "n_obs", "chunk_size"), envir = environment())
			parallel::clusterEvalQ(cl, {
				invisible(lapply(loaded_pkgs, require, character.only = TRUE))
			})
			chunk_results = parallel::parLapply(cl, chunk_starts, predict_chunk)
		} else {
			chunk_results = parallel::mclapply(chunk_starts, predict_chunk, mc.cores = num_workers)
		}
		for (chunk in chunk_results){
			chunk_len = chunk$end_idx - chunk$start_idx + 1L
			ice_curves[, chunk$start_idx : chunk$end_idx] = matrix(chunk$preds, nrow = n_obs, ncol = chunk_len)
		}
		if (verbose){cat(rep(".", grid_len), sep = "")}
	} else {
		for (start_idx in chunk_starts){
			end_idx = min(grid_len, start_idx + chunk_size - 1L)
			chunk_grid = grid_pts[start_idx : end_idx]
			X_chunk = X[rep(seq_len(n_obs), times = length(chunk_grid))]
			set(X_chunk, j = pred_name, value = rep(chunk_grid, each = n_obs))
			if (use_generic){
				preds = do.call("predict", c(list(object, X_chunk), predict_args))
			}
			else{
				preds = predictfcn(object = object, newdata = X_chunk)
			}
			ice_curves[, start_idx : end_idx] = matrix(preds, nrow = n_obs, ncol = length(chunk_grid))
			if (verbose){cat(rep(".", length(chunk_grid)), sep = "")}
		}
	}
	#return X to its original state.
	set(X, j = pred_name, value = xvec_temp)

	#do logit if necessary
	if (logodds || probit){
		#prevent log(0) error
		min_val = min(ice_curves)
		max_val = max(ice_curves)
		if (min_val < 0){
			stop("logodds is TRUE but predict returns negative values (these should be probabilities!)")
		} 
		if (min_val == 0){
			second_lowest = min(ice_curves[ice_curves > 0])
			if (is.na(second_lowest)){ 
				second_lowest = .0001 #arbitrary epsilon value
			}
			lowest_practical_prob = mean(c(second_lowest, 0)) #arbitrarily, halfway between 0 and second_lowest
			ice_curves[(ice_curves == 0)] = lowest_practical_prob
			if (verbose) warning("At least one probability was predicted to be 0, ICEbox is using ", lowest_practical_prob, " for the value(s) instead.")
		} 
		if (max_val == 1){
			second_highest = max(ice_curves[ice_curves < 1])
			if (is.na(second_highest)){ 
				second_highest = .9999 #arbitrary 1 - epsilon value
			} 
			highest_practical_prob = mean(c(second_highest, 1)) #arbitrarily, halfway between 1 and second_highest
			ice_curves[(ice_curves == 1)] = highest_practical_prob
			if (verbose) warning("At least one probability was predicted to be 1, ICEbox is using ", highest_practical_prob, " for the value(s) instead.")
		}
		
		# Efficient in-place transformation via Rcpp
		method = ifelse(logodds, 1L, 2L)
		ice_curves = transform_ice_curves_cpp(ice_curves, method, n_cores = num_cores)
		colnames(ice_curves) = grid_pts # restore colnames lost during transformation
	}
	if (verbose){cat("\n")}
	
	if (!inherits(predictor, "character")){
		xlab = paste("x", predictor, sep = "_")  #x_1, x_2 etc.
	} else {
		xlab = predictor #the actual name of the feature.
	}
	
	if (num_unique_pts <= MAX_NUM_UNIQUE_PTS_NOMINAL){
		nominal_axis = TRUE
	} else {
		nominal_axis = FALSE
	}	

	range_y = NULL
	sd_y = NULL
	if(!missing(y)){
		range_y = max(y) - min(y)
		sd_y = sd(y)
	} else if (!logodds && !probit){
		range_y = (max(ice_curves) - min(ice_curves))
		sd_y = sd(actual_prediction)
		if (verbose) cat("y not passed, so range_y is range of ice curves and sd_y is sd of predictions on real observations\n")
	}

	#compute friedman's pdp:
	pdp = colMeans(ice_curves) # pdp = average over the columns
  		if(missing(predictfcn)){
	    predictfcn=NULL
	}
	
	ice_obj = list(ice_curves = ice_curves, gridpts = grid_pts, predictor = predictor, xj = xj, actual_prediction = actual_prediction, 
				logodds = logodds, probit = probit, xlab = xlab, nominal_axis = nominal_axis, range_y = range_y, sd_y = sd_y, Xice = X, pdp = pdp,
				indices_to_build = indices_to_build, frac_to_build = frac_to_build, predictfcn = predictfcn) 
	class(ice_obj) = "ice"
		
	invisible(ice_obj)
}