#' Print method for \code{ice} objects.
#'
#' Prints a summary of an \code{ice} object.
#'
#' @param x Object of class \code{ice}.
#' @param ... Ignored for now.
#'
#' @export
print.ice = function(x, ...){
	assert_class(x, "ice")
	assert_list(list(...))
	cat("ice object generated on data with n = ", nrow(x$ice_curves), " for predictor \"", x$predictor, "\"\n", sep = "")
	cat("predictor considered ", ifelse(x$nominal_axis, "discrete", "continuous"), ", logodds ", ifelse(x$logodds, "on", "off"), "\n", sep = "")
}
