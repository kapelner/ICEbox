library(ICEbox)
library(sfsmisc)

set.seed(123)
x <- sort(runif(100, 0, 10))
y <- sin(x) + rnorm(100, 0, 0.1)

# Old D1tr
d_old <- D1tr(y, x)

# New derivative_cpp
d_new <- as.vector(ICEbox:::derivative_cpp(matrix(y, nrow=1), x))

cat("Correlation:", cor(d_old, d_new), "\n")
cat("Mean diff:", mean(abs(d_old - d_new)), "\n")
cat("SD old:", sd(d_old), "\n")
cat("SD new:", sd(d_new), "\n")
