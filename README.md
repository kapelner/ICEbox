ICEbox
======

This is the official repository for the R code for the R package ICEbox --- Individual Conditional Expectation (ICE) plots, a tool for visualizing the model estimated by any supervised learning algorithm. Classical partial dependence plots (PDPs) help visualize the average partial relationship between the predicted response and one or more features. In the presence of substantial interaction effects, the partial response relationship can be heterogeneous and an average curve, such as the PDP, can obfuscate the complexity of the modeled relationship. Accordingly, ICE plots refine the partial dependence plot by graphing the functional relationship between the predicted response and the feature for individual observations. Specifically, ICE plots highlight the variation in the fitted values across the range of a covariate, suggesting where and to what extent they may exist. In addition to providing a plotting suite for exploratory analysis, we include a visual test for additive structure in the data generating model. Through simulated examples and real data sets, we demonstrate how ICE plots can shed light on estimated models in ways PDPs cannot.

For a full writeup of the method with examples, see [Goldstein et al (2015)](https://www.tandfonline.com/doi/abs/10.1080/10618600.2014.907095) or the [free arxiv version](https://arxiv.org/abs/1309.6392) which is largely the same. 

Recent News
---------

January, 2026 - version 1.2 released

* upgraded to ggplot2
* rotated the additivityLineup test into the package (was previously an experimental script) as well as its backfitter function
* speedups via data.table
* speedups via optional parallelization for ice, dice, clusterICE
* speedup via Savitzky-Golay Filter for numerical derivative calculation implemented in Rcpp
* other speedups in Rcpp
* added Roxygen documentation
* checkmate checks for all functions' arguments
* tests
* benchmark vs previous version

