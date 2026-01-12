#' @import data.table
#' @import checkmate
#' @import ggplot2
#' @importFrom grDevices rgb
#' @importFrom stats ecdf kmeans predict qnorm qlogis quantile rbinom runif sd setNames supsmu
#' @importFrom utils methods
#' @importFrom Rcpp evalCpp
#' @useDynLib ICEbox, .registration = TRUE
utils::globalVariables(c(
	"cluster_id",
	"color",
	"curve_id",
	"label",
	"line_size",
	"x",
	"y",
	".original_indices",
	".I",
	":="
))
