# ICEbox 1.2.0.9000

* Added multi-platform GitHub Actions checks and automated pkgdown deployment.
* Added a pkgdown website with a curated function reference and installation
  documentation.
* Added Citation File Format and CodeMeta records, and expanded the package
  citation to cover both the software and its methodological article.
* Added regression tests for ICE and d-ICE construction, input validation, and
  compiled helpers.
* Removed compiled build artifacts from version control, replaced the malformed
  repository license copy with canonical GPL-2 and GPL-3 texts, and modernized
  package-development and discovery metadata.

# ICEbox 1.2

* Migrated plotting to ggplot2.
* Added `additivityLineup()` and its supporting `backfitter()` function.
* Improved performance with data.table and optional parallelization in
  `ice()`, `dice()`, and `clusterICE()`.
* Added an Rcpp implementation of the Savitzky-Golay filter for numerical
  derivative calculation and moved other performance-sensitive operations to
  Rcpp.
* Added roxygen2 documentation, argument validation with checkmate, tests, and
  benchmarks against the previous version.

# ICEbox 1.1.5

* `plot.dice()` uses grid points by default.

# ICEbox 1.1.4

* `plot.ice()` uses grid points by default.
* Made ICEbox more flexible following a suggestion from Dave Armstrong.

# ICEbox 1.1.3

* Fixed a bug in `plot.ice()` when `plot_pdp = FALSE`.
* Fixed a bug in `plot.ice()` when specifying `pts_preds_size`.
* The PDP in an ICE object remains based on all ICE curves unless
  `plot_points_indices` is specified.

# ICEbox 1.1.2

* Made `colorvec` consistent.

# ICEbox 1.1.1

* Fixed several bugs.

# ICEbox 1.1

* Fixed a segmentation fault when using GBM that originated in R's methods
  dispatch.
* Added a console legend for line colors when using `color_by`.
* Added handling for predicted classification probabilities of zero or one.
* Added an error when a prediction function returns factor levels instead of
  probabilities.
* Added probit plots for classification models.
* Added coloring by arbitrary data vectors, including the response.
* Added point labels and controls for selecting which points to plot.
* Fixed various other bugs.

# ICEbox 1.0

* Improved the color scheme, rug tick customization, and axis scaling.
* Allowed additional graphical parameters to be passed to plotting methods.
* Added citation information for the ICE publication.

# ICEbox 0.2

* Initial release.
