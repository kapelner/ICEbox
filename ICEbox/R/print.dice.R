#' Print method for \code{dice} objects.
#'
#' Prints a summary of a \code{dice} object.
#'
#' @param x Object of class \code{dice}.
#' @param ... Ignored for now.
#'
#' @export
print.dice = function(x, ...){
	assert_class(x, "dice")
	assert_list(list(...))
	cat("dice object generated on data with n = ", nrow(x$d_ice_curves), " for predictor \"", x$predictor, "\"\n", sep = "")
	cat("predictor considered ", ifelse(x$nominal_axis, "discrete", "continuous"), ", logodds ", ifelse(x$logodds, "on", "off"), "\n", sep = "")
}
