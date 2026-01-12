library(ICEbox)
library(sfsmisc)

set.seed(123)
# Non-equidistant grid
x <- sort(runif(100, 0, 10))
y <- sin(x) + rnorm(100, 0, 0.1)

# Old method
y_smooth_old <- supsmu(x, y)$y
deriv_old <- D1tr(y_smooth_old, x)

# New SG method

Y <- matrix(y, nrow = 1)

window_size <- 21

deriv_sg <- as.vector(ICEbox::sg_smooth_cpp(Y, window_size = window_size, order = 3, deriv = 1))





# Scaling fix

avg_spacing <- mean(diff(x))

deriv_sg_scaled <- deriv_sg / avg_spacing



# Compare

cat("Correlation between old and new derivative:", cor(deriv_old, deriv_sg_scaled), "\n")

cat("Mean absolute difference (scaled):", mean(abs(deriv_old - deriv_sg_scaled)), "\n")

cat("SD of old deriv:", sd(deriv_old), "\n")

cat("SD of SG deriv (scaled):", sd(deriv_sg_scaled), "\n")
