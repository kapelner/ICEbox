# ICEbox

ICEbox provides tools for constructing and visualizing Individual Conditional
Expectation (ICE) plots. ICE plots show the fitted relationship between a
feature and the predicted response for individual observations, revealing
heterogeneity and interactions that an averaged partial dependence plot can
hide.

For the method and worked examples, see Goldstein, Kapelner, Bleich, and Pitkin
(2015), [“Peeking Inside the Black Box: Visualizing Statistical Learning With
Plots of Individual Conditional
Expectation”](https://doi.org/10.1080/10618600.2014.907095), or the freely
available [arXiv version](https://arxiv.org/abs/1309.6392).

## Installation

Install the CRAN release with:

```r
install.packages("ICEbox")
```

Or install the latest build from Adam Kapelner's R-universe:

```r
install.packages(
  "ICEbox",
  repos = c(
    kapelner = "https://kapelner.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
```

## Basic use

```r
library(ICEbox)

fit <- lm(mpg ~ wt + hp + disp, data = mtcars)
X <- mtcars[c("wt", "hp", "disp")]

ice_fit <- ice(fit, X = X, y = mtcars$mpg, predictor = "wt", verbose = FALSE)
plot(ice_fit)
```
