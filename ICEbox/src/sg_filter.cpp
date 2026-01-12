#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Helper function to calculate factorial
double factorial(int n) {
    if (n <= 1) return 1.0;
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

// Savitzky-Golay coefficients calculator
NumericVector sg_coeffs(int m, int order, int deriv) {
    int window_size = 2 * m + 1;
    int n_params = order + 1;
    NumericMatrix A(window_size, n_params);
    
    for (int i = -m; i <= m; ++i) {
        for (int j = 0; j <= order; ++j) {
            A(i + m, j) = std::pow(i, j);
        }
    }
    
    NumericMatrix AtA(n_params, n_params);
    NumericMatrix At(n_params, window_size);
    
    for(int i=0; i<window_size; ++i) {
        for(int j=0; j<n_params; ++j) {
            At(j, i) = A(i, j);
        }
    }
    
    for(int i=0; i<n_params; ++i) {
        for(int j=0; j<n_params; ++j) {
            double sum = 0;
            for(int k=0; k<window_size; ++k) {
                sum += At(i, k) * A(k, j);
            }
            AtA(i, j) = sum;
        }
    }
    
    NumericMatrix Inv(n_params, n_params);
    for(int i=0; i<n_params; ++i) Inv(i, i) = 1.0;
    
    for(int i=0; i<n_params; ++i) {
        double pivot = AtA(i, i);
        for(int j=0; j<n_params; ++j) {
            AtA(i, j) /= pivot;
            Inv(i, j) /= pivot;
        }
        for(int k=0; k<n_params; ++k) {
            if(k != i) {
                double factor = AtA(k, i);
                for(int j=0; j<n_params; ++j) {
                    AtA(k, j) -= factor * AtA(i, j);
                    Inv(k, j) -= factor * Inv(i, j);
                }
            }
        }
    }
    
    NumericVector coeffs(window_size);
    for(int i=0; i<window_size; ++i) {
        double sum = 0;
        for(int j=0; j<n_params; ++j) {
            sum += Inv(deriv, j) * At(j, i);
        }
        coeffs[i] = sum * factorial(deriv);
    }
    
    return coeffs;
}


//' Savitzky-Golay Filter for Matrix (Row-wise)
//'
//' Smooths each row of a matrix using a Savitzky-Golay filter.
//' @param x Matrix to smooth row-wise
//' @param window_size Size of the filter window (must be odd)
//' @param order Polynomial order
//' @param deriv Derivative order (0=smooth, 1=first deriv, etc.)
//' @param n_cores Number of cores to use
//' @export
// [[Rcpp::export]]
NumericMatrix sg_smooth_cpp(NumericMatrix x, int window_size, int order, int deriv, int n_cores = 1) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    
    if (window_size % 2 == 0) window_size++;
    if (window_size > ncol) window_size = ncol | 1;
    
    int m = (window_size - 1) / 2;
    NumericVector coeffs = sg_coeffs(m, order, deriv);
    NumericMatrix result(nrow, ncol);
    
    #pragma omp parallel for num_threads(n_cores)
    for (int i = 0; i < nrow; ++i) {
        // Left boundary
        for (int j = 0; j < m; ++j) {
            double val = 0;
            for (int k = 0; k < window_size; ++k) {
                int idx = j + k - m;
                if (idx < 0) idx = 0;
                val += x(i, idx) * coeffs[k];
            }
            result(i, j) = val;
        }
        // Middle (no boundary checks needed)
        for (int j = m; j < ncol - m; ++j) {
            double val = 0;
            for (int k = 0; k < window_size; ++k) {
                val += x(i, j + k - m) * coeffs[k];
            }
            result(i, j) = val;
        }
        // Right boundary
        for (int j = ncol - m; j < ncol; ++j) {
            double val = 0;
            for (int k = 0; k < window_size; ++k) {
                int idx = j + k - m;
                if (idx >= ncol) idx = ncol - 1;
                val += x(i, idx) * coeffs[k];
            }
            result(i, j) = val;
        }
    }
    
    return result;
}

//' Efficient Column Standard Deviations
//' 
//' @param x Numeric Matrix
//' @param n_cores Number of cores to use
//' @export
// [[Rcpp::export]]
NumericVector colSds_cpp(NumericMatrix x, int n_cores = 1) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector res(ncol);
    
    if (nrow < 2) {
        std::fill(res.begin(), res.end(), NA_REAL);
        return res;
    }

    #pragma omp parallel for num_threads(n_cores)
    for (int j = 0; j < ncol; ++j) {
        double mean = 0;
        for (int i = 0; i < nrow; ++i) {
            mean += x(i, j);
        }
        mean /= nrow;
        
        double variance = 0;
        for (int i = 0; i < nrow; ++i) {
            double diff = x(i, j) - mean;
            variance += diff * diff;
        }
        res[j] = std::sqrt(variance / (nrow - 1));
    }
    
    return res;
}

//' Row-wise Centering
//' 
//' Centers each row of a matrix by subtracting the row mean.
//' @param x Numeric Matrix
//' @param n_cores Number of cores to use
//' @export
// [[Rcpp::export]]
NumericMatrix rowCenter_cpp(NumericMatrix x, int n_cores = 1) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericMatrix res(nrow, ncol);
    
    #pragma omp parallel for num_threads(n_cores)
    for (int i = 0; i < nrow; ++i) {
        double sum = 0;
        for (int j = 0; j < ncol; ++j) {
            sum += x(i, j);
        }
        double mean = sum / ncol;
        for (int j = 0; j < ncol; ++j) {
            res(i, j) = x(i, j) - mean;
        }
    }
    return res;
}

//' Probability Transformation
//' 
//' Efficiently applies logodds or probit transformation to a matrix.
//' @param x Numeric Matrix (probabilities)
//' @param method 1 for centered logodds, 2 for probit
//' @param n_cores Number of cores to use
//' @export
// [[Rcpp::export]]
NumericMatrix transform_ice_curves_cpp(NumericMatrix x, int method, int n_cores = 1) {
    int n = x.length();
    NumericMatrix res(x.nrow(), x.ncol());
    
    #pragma omp parallel for num_threads(n_cores)
    for (int i = 0; i < n; ++i) {
        double p = x[i];
        if (method == 1) { // Centered Logodds: 0.5 * log(p / (1-p))
            res[i] = 0.5 * std::log(p / (1.0 - p));
        } else { // Probit
            res[i] = R::qnorm(p, 0.0, 1.0, 1, 0);
        }
    }
    return res;
}

//' Melt Matrix to Long Format Vector
//' 
//' Efficiently converts a matrix to a long-format vector (row-major order)
//' for plotting.
//' @param x Numeric Matrix
//' @param n_cores Number of cores to use
//' @export
// [[Rcpp::export]]
NumericVector melt_ice_curves_cpp(NumericMatrix x, int n_cores = 1) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector res(nrow * ncol);
    
    #pragma omp parallel for num_threads(n_cores)
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            res[i * ncol + j] = x(i, j);
        }
    }
    return res;
}

//' Efficient Numerical Derivative for Matrix (Row-wise)
//' 
//' Computes the first derivative using centered differences, mirroring sfsmisc::D1tr.
//' @param x Numeric Matrix (smoothed values)
//' @param gridpts Grid points corresponding to columns of x
//' @param n_cores Number of cores to use
//' @export
// [[Rcpp::export]]
NumericMatrix derivative_cpp(NumericMatrix x, NumericVector gridpts, int n_cores = 1) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericMatrix res(nrow, ncol);
    
    if (ncol < 2) return res;

    #pragma omp parallel for num_threads(n_cores)
    for (int i = 0; i < nrow; ++i) {
        // Forward difference for first point
        res(i, 0) = (x(i, 1) - x(i, 0)) / (gridpts[1] - gridpts[0]);
        
        // Centered difference for middle points: (y[j+1] - y[j-1]) / (x[j+1] - x[j-1])
        for (int j = 1; j < ncol - 1; ++j) {
            res(i, j) = (x(i, j + 1) - x(i, j - 1)) / (gridpts[j + 1] - gridpts[j - 1]);
        }
        
        // Backward difference for last point
        res(i, ncol - 1) = (x(i, ncol - 1) - x(i, ncol - 2)) / (gridpts[ncol - 1] - gridpts[ncol - 2]);
    }
    return res;
}
