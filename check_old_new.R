library(ICEbox)
library(sfsmisc)
library(randomForest)

set.seed(123)
N <- 1000
P <- 5
X <- as.data.frame(matrix(rnorm(N * P), ncol = P))
colnames(X) <- paste0("x", 1:P)
y <- X$x1 + X$x2^2 + rnorm(N)
rf_mod <- randomForest(X, y, ntree=10)
ice_obj <- ice(object = rf_mod, X = X, y = y, predictor = "x1", verbose = FALSE)

# Old method
gridpts <- ice_obj$gridpts
old_derivs <- t(apply(ice_obj$ice_curves, 1, function(y) D1tr(supsmu(gridpts, y)$y, gridpts)))
cat("Old average derivative:", mean(old_derivs), "\n")

# New method
avg_spacing <- (max(gridpts) - min(gridpts)) / (length(gridpts) - 1)
dice_obj <- dice(ice_obj, sg_window_size = 101) # Will use order 2
# Note: current dice.R uses derivative_cpp. I want to test sg_smooth_cpp(deriv=1)
d_sg_direct <- ICEbox::sg_smooth_cpp(ice_obj$ice_curves, window_size = 21, order = 2, deriv = 1) / avg_spacing
cat("New average derivative (SG direct, window=21, order=2):", mean(d_sg_direct), "\n")

