#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

// Structure to represent a sparse matrix in CSR format
struct CSRMatrix {
    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    std::vector<double> values;
    int n; // number of rows/columns
};

// Function to compute PageRank using parallel Gauss-Seidel with blocking
void pagerank_gauss_seidel_blocking(const CSRMatrix& A, 
                                   std::vector<double>& x, 
                                   double alpha, 
                                   double tol, 
                                   int max_iter,
                                   int block_size) {
    int n = A.n;
    x.assign(n, 1.0 / n); // Initialize with uniform distribution
    
    std::vector<double> x_old(n);
    double norm_diff;
    int iter = 0;
    
    // Calculate number of blocks
    int num_blocks = (n + block_size - 1) / block_size;
    
    do {
        x_old = x;
        norm_diff = 0.0;
        
        // Process each block
        #pragma omp parallel for reduction(+:norm_diff)
        for (int block = 0; block < num_blocks; ++block) {
            int start = block * block_size;
            int end = std::min((block + 1) * block_size, n);
            
            // Process each row in the block
            for (int i = start; i < end; ++i) {
                double sum = 0.0;
                
                // Sum contributions from incoming links
                for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                    int col = A.col_idx[j];
                    sum += A.values[j] * x[col];
                }
                
                // Update PageRank value using Gauss-Seidel formula
                x[i] = alpha * sum + (1.0 - alpha) / n;
                
                // Update norm difference
                norm_diff += std::abs(x[i] - x_old[i]);
            }
        }
        
        iter++;
    } while (norm_diff > tol && iter < max_iter);
    
    std::cout << "Converged in " << iter << " iterations with norm_diff = " << norm_diff << std::endl;
}

int main() {
    // Example usage
    CSRMatrix A;
    A.n = 4;
    A.row_ptr = {0, 2, 3, 4, 5};
    A.col_idx = {1, 2, 0, 1, 2};
    A.values = {0.5, 0.5, 1.0, 1.0, 1.0};
    
    std::vector<double> x;
    double alpha = 0.85;
    double tol = 1e-6;
    int max_iter = 1000;
    int block_size = 32; // Adjust based on cache size
    
    pagerank_gauss_seidel_blocking(A, x, alpha, tol, max_iter, block_size);
    
    // Print results
    std::cout << "PageRank values:" << std::endl;
    for (int i = 0; i < A.n; ++i) {
        std::cout << "Node " << i << ": " << x[i] << std::endl;
    }
    
    return 0;
}
