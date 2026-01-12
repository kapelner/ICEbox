#' Summary function for \code{dice} objects.
#'
#' Alias of \code{print} method.
#'
#' @param object Object of class \code{dice}.
#' @param ... Ignored for now.
#'
#' @export
summary.dice = function(object, ...){
	assert_class(object, "dice")
	assert_list(list(...))
	print(object, ...)
}
