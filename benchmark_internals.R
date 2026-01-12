args <- commandArgs(trailingOnly = TRUE)
lib_path <- args[1]
result_csv <- args[2]
plot_pdf <- args[3]
n_iter <- as.integer(args[4])
version_label <- args[5]
task <- args[6]

# Ensure we use the custom library
.libPaths(c(lib_path, .libPaths()))

if (!require("ICEbox", quietly = TRUE)) {
  stop("ICEbox not found in ", lib_path)
}
if (!require("randomForest", quietly = TRUE)) {
  stop("randomForest not found")
}
if (!require("MASS", quietly = TRUE)) {
    # MASS is needed for Boston data
    warning("MASS package not found, Boston tasks might fail")
}

# --- Setup Data and Model ---

setup_wine <- function() {
    data("WhiteWine", package = "ICEbox")
    # Subset for speed but enough to be meaningful
    set.seed(42)
    df <- WhiteWine[sample(nrow(WhiteWine), 500), ]
    rf <- randomForest(quality ~ ., data = df, ntree = 50)
    list(data = df, model = rf, predictor = "alcohol", y = df$quality)
}

setup_boston <- function() {
    # Boston is in MASS
    data("Boston", package = "MASS")
    set.seed(42)
    # Boston is small (506 rows), use all or subset. Use all.
    df <- Boston
    rf <- randomForest(medv ~ ., data = df, ntree = 50)
    list(data = df, model = rf, predictor = "rm", y = df$medv)
}

# --- Run Task ---

runtimes <- numeric(n_iter)
pdf(plot_pdf)

if (task == "wine_ice") {
    setup <- setup_wine()
    # Benchmark ICE
    for (i in 1:n_iter) {
        start <- Sys.time()
        ice_obj <- ice(object = setup$model, X = setup$data, y = setup$y, predictor = setup$predictor, frac_to_build = 0.5, verbose = FALSE)
        end <- Sys.time()
        runtimes[i] <- as.numeric(difftime(end, start, units = "secs"))
        
        if (i == 1) {
            # Plot Uncentered
            plot(ice_obj, main = paste(version_label, "- Wine ICE Uncentered"), centered = FALSE)
            # Plot Centered
            plot(ice_obj, main = paste(version_label, "- Wine ICE Centered"), centered = TRUE)
        }
    }
} else if (task == "wine_dice") {
    setup <- setup_wine()
    # Pre-compute ICE for DICE benchmark
    ice_obj <- ice(object = setup$model, X = setup$data, y = setup$y, predictor = setup$predictor, frac_to_build = 0.5, verbose = FALSE)
    
    # Benchmark DICE
    for (i in 1:n_iter) {
        start <- Sys.time()
        if (grepl("New", version_label)) {
            dice_obj <- dice(ice_obj, verbose = FALSE, use_supsmu = TRUE)
        } else {
            dice_obj <- dice(ice_obj)
        }
        end <- Sys.time()
        runtimes[i] <- as.numeric(difftime(end, start, units = "secs"))
        
        if (i == 1) {
             # Plot Uncentered
            plot(dice_obj, main = paste(version_label, "- Wine DICE Uncentered"), centered = FALSE)
            # Plot Centered
            plot(dice_obj, main = paste(version_label, "- Wine DICE Centered"), centered = TRUE)
        }
    }

} else if (task == "wine_cluster") {
    setup <- setup_wine()
    # Pre-compute ICE
    ice_obj <- ice(object = setup$model, X = setup$data, y = setup$y, predictor = setup$predictor, frac_to_build = 0.5, verbose = FALSE)
    
    for (i in 1:n_iter) {
        start <- Sys.time()
        # clusterICE plots by default
        if (i == 1) {
            # Capture the plots
             # Page 1: Uncentered
             if (grepl("New", version_label)) {
                clusterICE(ice_obj, nClusters = 3, plot = TRUE, centered = FALSE, main = paste(version_label, "- Wine ClusterICE (Uncentered)"))
             } else {
                clusterICE(ice_obj, nClusters = 3, plot = TRUE, centered = FALSE)
                title(paste(version_label, "- Wine ClusterICE (Uncentered)"))
             }
             
             # Page 2: Centered
             if (grepl("New", version_label)) {
                clusterICE(ice_obj, nClusters = 3, plot = TRUE, centered = TRUE, main = paste(version_label, "- Wine ClusterICE (Centered)"))
             } else {
                clusterICE(ice_obj, nClusters = 3, plot = TRUE, centered = TRUE)
                title(paste(version_label, "- Wine ClusterICE (Centered)"))
             }
        } else {
             # Don't plot for timing (approximate, since plotting is part of the function usually, 
             # but clusterICE has a plot argument. If we want to benchmark the calculation+plotting, we should keep plot=TRUE.
             # However, plotting to a PDF device 20 times is slow.
             # Let's benchmark just the clustering part if possible, or accept plotting overhead.
             # The user asked for runtimes.
             # clusterICE(..., plot=FALSE) does just the kmeans.
             # But the fix was in the plotting code. 
             # We should probably benchmark the plotting too if we want to catch regressions there, 
             # but usually benchmarks exclude rendering time.
             # Given the "size" fix is in the plotting code, we MUST run the plotting code to verify it doesn't crash/warn.
             # But we can suppress output to a NULL device for iterations > 1 if we want speed,
             # OR just run it with plot=FALSE for speed and rely on i=1 for verification.
             # The 'size' warning happens during plot construction.
             # Let's run plot=FALSE for the measured iterations to measure algorithm speed, 
             # and plot=TRUE for the first one for visual verification.
             # Wait, if we only plot once, we only check for warnings once. That's fine.
             clusterICE(ice_obj, nClusters = 3, plot = FALSE, centered = TRUE)
        }
        end <- Sys.time()
        runtimes[i] <- as.numeric(difftime(end, start, units = "secs"))
    }

} else if (task == "boston_ice") {
    setup <- setup_boston()
    for (i in 1:n_iter) {
        start <- Sys.time()
        ice_obj <- ice(object = setup$model, X = setup$data, y = setup$y, predictor = setup$predictor, frac_to_build = 0.5, verbose = FALSE)
        end <- Sys.time()
        runtimes[i] <- as.numeric(difftime(end, start, units = "secs"))
        
        if (i == 1) {
            plot(ice_obj, main = paste(version_label, "- Boston ICE Uncentered"), centered = FALSE)
            plot(ice_obj, main = paste(version_label, "- Boston ICE Centered"), centered = TRUE)
        }
    }
} else if (task == "boston_dice") {
    setup <- setup_boston()
    ice_obj <- ice(object = setup$model, X = setup$data, y = setup$y, predictor = setup$predictor, frac_to_build = 0.5, verbose = FALSE)
    
    for (i in 1:n_iter) {
        start <- Sys.time()
        if (grepl("New", version_label)) {
            dice_obj <- dice(ice_obj, verbose = FALSE, use_supsmu = TRUE)
        } else {
            dice_obj <- dice(ice_obj)
        }
        end <- Sys.time()
        runtimes[i] <- as.numeric(difftime(end, start, units = "secs"))
        
        if (i == 1) {
            plot(dice_obj, main = paste(version_label, "- Boston DICE Uncentered"), centered = FALSE)
            plot(dice_obj, main = paste(version_label, "- Boston DICE Centered"), centered = TRUE)
        }
    }
} else {
    stop("Unknown task: ", task)
}

invisible(dev.off())

# Save results
write.csv(data.frame(version = version_label, task = task, time = runtimes), file = result_csv, row.names = FALSE)
