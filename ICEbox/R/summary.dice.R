#' Summary function for \code{ice} objects.
#'
#' Alias of \code{print} method.
#'
#' @param object Object of class \code{ice}.
#' @param ... Ignored for now.
#'
#' @export
summary.ice = function(object, ...){
	assert_class(object, "ice")
	assert_list(list(...))
	print(object,...)
}
