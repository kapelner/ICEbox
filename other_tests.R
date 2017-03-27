library(ICEbox)
library(gbm)

n = 250
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
y <- x1*x2*x3 + rnorm(n, 0.1)
X = data.frame(x1, x2, x3)


object <- gbm.fit(X, y, distribution = "gaussian")
oice <- ice(object, X, y,predictor = "x1", n.trees = 100)
plot(oice, color_by = "x2", frac_to_plot = 0.99)


library(randomForest)
yb = factor(ifelse(y > median(y), 1, 0))
object <- randomForest(X, yb)
oicep <- ice(object, X, predictor = "x1", probit = TRUE, predictfcn = function(object, newdata){predict(object, newdata, type = "prob")[, 2]})
plot(oicep, color_by = "x2")
plot(oicep, color_by = "x2", point_labels = 1 : n)
plot(oicep, color_by = yb, point_labels = 1 : n)
plot(oicep, color_by = y, point_labels = 1 : n ,frac_to_plot = 0.1)
plot(oicep, color_by = yb, point_labels = 1 : n , frac_to_plot = 0.02)

par(mfrow = c(2, 2))
idx = which(x3 < quantile(x3, .25))
p1=plot(oicep, color_by = yb, point_labels = 1 : n , plot_points_indices = idx, point_labels_size = 1, centered = TRUE)
idx = which(x3 >= quantile(x3, .25) & x3 < quantile(x3, .5))
p2=plot(oicep, color_by = yb, point_labels = 1 : n , plot_points_indices = idx, point_labels_size = 1, centered = TRUE)
idx = which(x3 >= quantile(x3, .5) & x3 < quantile(x3, .75))
p3=plot(oicep, color_by = yb, point_labels = 1 : n , plot_points_indices = idx, point_labels_size = 1, centered = TRUE)
idx = which(x3 >= quantile(x3, .75))
p4=plot(oicep, color_by = yb, point_labels = 1 : n , plot_points_indices = idx, point_labels_size = 1, centered = TRUE)

all.equal(p1$pdp, p4$pdp)
windows()
oicel <- ice(object, X, predictor = "x1", logodds = TRUE, predictfcn = function(object, newdata){predict(object, newdata, type = "prob")[, 2]})
plot(oicel, color_by = "x2")