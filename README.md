# ICEbox: Individual Conditional Expectation Plots <img src="ICEbox/man/figures/logo.png" align="right" height="139" alt="ICEbox hex logo" />

[![CRAN status](https://www.r-pkg.org/badges/version/ICEbox)](https://CRAN.R-project.org/package=ICEbox)
[![R-CMD-check](https://github.com/kapelner/ICEbox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kapelner/ICEbox/actions/workflows/R-CMD-check.yaml)
[![R-universe version](https://kapelner.r-universe.dev/ICEbox/badges/version)](https://kapelner.r-universe.dev/ICEbox)
[![R-universe checks](https://kapelner.r-universe.dev/ICEbox/badges/checks)](https://kapelner.r-universe.dev/ICEbox)
[![License: GPL-2 or GPL-3](https://img.shields.io/badge/license-GPL--2%20or%20GPL--3-blue.svg)](LICENSE)

This is the official repository for the R code for the R package ICEbox --- Individual Conditional Expectation (ICE) plots, a tool for visualizing the model estimated by any supervised learning algorithm. Classical partial dependence plots (PDPs) help visualize the average partial relationship between the predicted response and one or more features. In the presence of substantial interaction effects, the partial response relationship can be heterogeneous and an average curve, such as the PDP, can obfuscate the complexity of the modeled relationship. Accordingly, ICE plots refine the partial dependence plot by graphing the functional relationship between the predicted response and the feature for individual observations. Specifically, ICE plots highlight the variation in the fitted values across the range of a covariate, suggesting where and to what extent they may exist. In addition to providing a plotting suite for exploratory analysis, we include a visual test for additive structure in the data generating model. Through simulated examples and real data sets, we demonstrate how ICE plots can shed light on estimated models in ways PDPs cannot.

For a full writeup of the method with examples, see [Goldstein et al (2015)](https://www.tandfonline.com/doi/abs/10.1080/10618600.2014.907095) or the [free arxiv version](https://arxiv.org/abs/1309.6392) which is largely the same. 

Installation
------------

Install the CRAN release with:

```r
install.packages("ICEbox")
```

Or install the latest R-universe build with:

```r
install.packages(
  "ICEbox",
  repos = c(
    kapelner = "https://kapelner.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
```

The package website is at <https://kapelner.github.io/ICEbox/>.

Citation
--------

To cite the software and the methodological article, run:

```r
citation("ICEbox")
```

The methodological article is Goldstein, Kapelner, Bleich, and Pitkin
(2015), “Peeking Inside the Black Box: Visualizing Statistical Learning With
Plots of Individual Conditional Expectation,”
<https://doi.org/10.1080/10618600.2014.907095>. The permanent CRAN package DOI
is <https://doi.org/10.32614/CRAN.package.ICEbox>.

License
-------

ICEbox is available under either the GNU General Public License version 2 or
version 3, at your option. See [LICENSE](LICENSE) for details.

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
