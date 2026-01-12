repo_root <- getwd()
# Absolute paths
old_repo_path <- normalizePath(file.path(repo_root, "../ICEbox_old/ICEbox"), mustWork = TRUE)
new_repo_path <- normalizePath(file.path(repo_root, "ICEbox"), mustWork = TRUE)
tmp_dir <- file.path(repo_root, "benchmark_env")
lib_old <- file.path(tmp_dir, "lib_old")
lib_new <- file.path(tmp_dir, "lib_new")
worker_script <- file.path(repo_root, "benchmark_internals.R")

n_iter <- 10 # Reduced iterations per task to keep total time reasonable

# Create libs
if (!dir.exists(lib_old)) dir.create(lib_old, recursive = TRUE)
if (!dir.exists(lib_new)) dir.create(lib_new, recursive = TRUE)

# Install Old
cat("Installing Old Version...\n")
install.packages(old_repo_path, lib = lib_old, repos = NULL, type = "source", quiet = TRUE)

# Install New
cat("Installing New Version...\n")
install.packages(new_repo_path, lib = lib_new, repos = NULL, type = "source", quiet = TRUE)

tasks <- c("wine_ice", "wine_dice", "wine_cluster", "boston_ice", "boston_dice")
results_list <- list()
pdf_list <- c()

gs_cmd <- Sys.which("gs")

for (task in tasks) {
    cat(sprintf("\n=== Running Task: %s ===\n", task))
    
    # Define files
    old_csv <- file.path(tmp_dir, paste0("old_", task, ".csv"))
    new_csv <- file.path(tmp_dir, paste0("new_", task, ".csv"))
    old_pdf <- file.path(tmp_dir, paste0("old_", task, ".pdf"))
    new_pdf <- file.path(tmp_dir, paste0("new_", task, ".pdf"))
    
    # Run Old
    cat("  Benchmarking Old Version...\n")
    exit_code <- system2("Rscript", args = c(worker_script, lib_old, old_csv, old_pdf, n_iter, shQuote("Old Plot"), task))
    if (exit_code != 0) stop(paste("Old version failed on task:", task))
    
    # Run New
    cat("  Benchmarking New Version...\n")
    exit_code <- system2("Rscript", args = c(worker_script, lib_new, new_csv, new_pdf, n_iter, shQuote("New Plot"), task))
    if (exit_code != 0) stop(paste("New version failed on task:", task))
    
    # Interleave Pages if gs is available
    if (gs_cmd != "") {
        base_name <- file.path(tmp_dir, task)
        p1_old <- paste0(base_name, "_old_p1.pdf")
        p2_old <- paste0(base_name, "_old_p2.pdf")
        p1_new <- paste0(base_name, "_new_p1.pdf")
        p2_new <- paste0(base_name, "_new_p2.pdf")

        # Split Old (Page 1 and Page 2)
        system2(gs_cmd, args = c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", paste0("-sOutputFile=", p1_old), "-dFirstPage=1", "-dLastPage=1", old_pdf))
        system2(gs_cmd, args = c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", paste0("-sOutputFile=", p2_old), "-dFirstPage=2", "-dLastPage=2", old_pdf))

        # Split New (Page 1 and Page 2)
        system2(gs_cmd, args = c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", paste0("-sOutputFile=", p1_new), "-dFirstPage=1", "-dLastPage=1", new_pdf))
        system2(gs_cmd, args = c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", paste0("-sOutputFile=", p2_new), "-dFirstPage=2", "-dLastPage=2", new_pdf))

        # Add to list in alternating order: Old P1, New P1, Old P2, New P2
        pdf_list <- c(pdf_list, p1_old, p1_new, p2_old, p2_new)
    } else {
        # Fallback if no gs: just append files
        pdf_list <- c(pdf_list, old_pdf, new_pdf)
    }
    
    # Read Results
    old_res <- read.csv(old_csv)
    new_res <- read.csv(new_csv)
    
    # Combine
    all_res <- rbind(old_res, new_res)
    results_list[[task]] <- all_res
}

# --- Analysis & Reporting ---

comp_pdf <- file.path(repo_root, "benchmark_comparison.pdf")
pdf(comp_pdf)

par(mfrow = c(2, 2)) # 2x2 grid for 4 tasks

for (task in tasks) {
    data <- results_list[[task]]
    
    # T-test
    t_res <- t.test(time ~ version, data = data)
    
    mean_old <- mean(data$time[data$version == "Old Plot"])
    mean_new <- mean(data$time[data$version == "New Plot"])
    
    # Boxplot
    boxplot(time ~ version, data = data, 
            main = paste0(task, "\nOld: ", round(mean_old, 3), "s | New: ", round(mean_new, 3), "s"),
            ylab = "Time (s)",
            col = c("lightblue", "lightgreen"),
            sub = paste("p =", format.pval(t_res$p.value, digits = 3)))
}

invisible(dev.off())

cat("\nBenchmark suite complete.\n")
cat("Comparison Summary:", comp_pdf, "\n")

# Merge all PDFs
final_merged_pdf <- file.path(repo_root, "benchmark_full_report.pdf")

# Gather all PDFs: Summary first, then task specific ones
all_pdfs <- c(comp_pdf, pdf_list)

if (gs_cmd != "") {
    cat("Merging all PDFs into", final_merged_pdf, "...\n")
    system2(gs_cmd, args = c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", paste0("-sOutputFile=", final_merged_pdf), all_pdfs))
    cat("Merged PDF created.\n")
} else {
    cat("Ghostscript not found. PDFs remain separate.\n")
}

# Cleanup
cat("Cleaning up temporary benchmark environment...\n")
unlink(tmp_dir, recursive = TRUE)