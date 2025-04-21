#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <algorithm>
#include <cstring>

// Structure to represent a sparse matrix in CSR format
struct CSRMatrix {
    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    std::vector<double> values;
    int n; // number of rows/columns
};

// Normalize the web graph to create the Markov transition matrix
void normalizeCSRMatrix(CSRMatrix& matrix) {
    std::vector<double> col_sums(matrix.n, 0.0);
    
    // Calculate out-degree (sum of outgoing links) for each node
    #pragma omp parallel for
    for (int i = 0; i < matrix.n; ++i) {
        for (int j = matrix.row_ptr[i]; j < matrix.row_ptr[i + 1]; ++j) {
            col_sums[i] += matrix.values[j];
        }
    }
    
    // Normalize each column to sum to 1
    #pragma omp parallel for
    for (int i = 0; i < matrix.n; ++i) {
        if (col_sums[i] > 0) {
            for (int j = matrix.row_ptr[i]; j < matrix.row_ptr[i + 1]; ++j) {
                matrix.values[j] /= col_sums[i];
            }
        }
    }
}

// Identify dangling nodes (nodes with no outgoing links)
std::vector<bool> identifyDanglingNodes(const CSRMatrix& matrix) {
    std::vector<bool> is_dangling(matrix.n, true);
    
    #pragma omp parallel for
    for (int i = 0; i < matrix.n; ++i) {
        if (matrix.row_ptr[i] < matrix.row_ptr[i + 1]) {
            is_dangling[i] = false;
        }
    }
    
    return is_dangling;
}

// PageRank with Gauss-Seidel iterations
std::vector<double> pageRankGaussSeidel(
    const CSRMatrix& matrix,
    double alpha = 0.85,
    double tolerance = 1e-8,
    int max_iterations = 100
) {
    int n = matrix.n;
    std::vector<double> rank(n, 1.0 / n); // Initialize ranks evenly
    std::vector<double> old_rank(n);
    std::vector<bool> dangling_nodes = identifyDanglingNodes(matrix);
    
    double dangling_factor = alpha / n;
    double teleport_factor = (1.0 - alpha) / n;
    double error = 1.0;
    int iterations = 0;
    
    while (error > tolerance && iterations < max_iterations) {
        // Copy current ranks to old_rank
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            old_rank[i] = rank[i];
        }
        
        // Calculate dangling node contribution
        double dangling_sum = 0.0;
        #pragma omp parallel for reduction(+:dangling_sum)
        for (int i = 0; i < n; ++i) {
            if (dangling_nodes[i]) {
                dangling_sum += old_rank[i];
            }
        }
        dangling_sum = dangling_factor * dangling_sum;
        
        // Gauss-Seidel iteration - NO PARALLELISM HERE due to dependencies
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            
            // Calculate contribution from incoming links
            for (int j = 0; j < n; ++j) {
                for (int k = matrix.row_ptr[j]; k < matrix.row_ptr[j + 1]; ++k) {
                    if (matrix.col_idx[k] == i) {
                        // For Gauss-Seidel, use updated values when available
                        if (j < i) {
                            sum += alpha * matrix.values[k] * rank[j];
                        } else {
                            sum += alpha * matrix.values[k] * old_rank[j];
                        }
                    }
                }
            }
            
            // Add teleportation and dangling nodes contribution
            rank[i] = sum + teleport_factor + dangling_sum;
        }
        
        // Calculate error (can be parallelized)
        error = 0.0;
        #pragma omp parallel for reduction(max:error)
        for (int i = 0; i < n; ++i) {
            double node_error = std::fabs(rank[i] - old_rank[i]);
            error = std::max(error, node_error);
        }
        
        // Normalize to ensure sum is 1 (parallelized in two steps)
        double sum = 0.0;
        #pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < n; ++i) {
            sum += rank[i];
        }
        
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            rank[i] /= sum;
        }
        
        iterations++;
        std::cout << "Iteration " << iterations << ", error: " << error << std::endl;
    }
    
    std::cout << "Converged after " << iterations << " iterations." << std::endl;
    return rank;
}

// Helper function to create a CSR matrix from an adjacency list
CSRMatrix createCSRFromAdjacencyList(const std::vector<std::vector<std::pair<int, double>>>& adj_list) {
    CSRMatrix matrix;
    int n = adj_list.size();
    matrix.n = n;
    
    matrix.row_ptr.resize(n + 1, 0);
    
    // Count number of non-zero elements per row and set row_ptr
    int nnz = 0;
    for (int i = 0; i < n; ++i) {
        matrix.row_ptr[i] = nnz;
        nnz += adj_list[i].size();
    }
    matrix.row_ptr[n] = nnz;
    
    matrix.col_idx.resize(nnz);
    matrix.values.resize(nnz);
    
    // Fill col_idx and values
    int idx = 0;
    for (int i = 0; i < n; ++i) {
        for (const auto& edge : adj_list[i]) {
            matrix.col_idx[idx] = edge.first;
            matrix.values[idx] = edge.second;
            idx++;
        }
    }
    
    return matrix;
}

int main(int argc, char** argv) {
    // Check input arguments
    if (argc < 5) {
        std::cout << "Error in arguments! Required: graph_filename N threshold damping_factor [num_threads]" << std::endl;
        std::cout << "Example: ./pagerank web-Google.txt 916428 0.00001 0.85 8" << std::endl;
        return 1;
    }

    // Parse command line arguments
    char filename[256];
    strcpy(filename, argv[1]);
    int N = atoi(argv[2]);               // Number of nodes
    double threshold = atof(argv[3]);    // Convergence threshold
    double alpha = atof(argv[4]);        // Damping factor (called d in the reference code)
    int num_threads = 1;                 // Default to single thread
    
    if (argc > 5) {
        num_threads = atoi(argv[5]);
    }
    
    // Set number of OpenMP threads
    omp_set_num_threads(num_threads);
    
    std::cout << "Graph file: " << filename << std::endl;
    std::cout << "Number of nodes: " << N << std::endl;
    std::cout << "Convergence threshold: " << threshold << std::endl;
    std::cout << "Damping factor: " << alpha << std::endl;
    std::cout << "Number of threads: " << num_threads << std::endl;
    
    // Read the graph data and build CSR matrix
    std::cout << "Reading graph data..." << std::endl;
    
    // Initialize empty adjacency list
    std::vector<std::vector<std::pair<int, double>>> adj_list(N);
    
    // Read the input file
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }
    
    std::string line;
    int from_idx, to_idx;
    int edges_count = 0;
    
    // Skip comment lines that start with #
    while (std::getline(infile, line)) {
        if (line[0] == '#') continue;
        
        std::stringstream ss(line);
        if (ss >> from_idx >> to_idx) {
            // Add edge with weight 1.0
            if (from_idx < N && to_idx < N) {
                adj_list[from_idx].push_back({to_idx, 1.0});
                edges_count++;
            }
        }
    }
    
    std::cout << "Loaded " << edges_count << " edges." << std::endl;
    
    // Create CSR matrix from adjacency list
    std::cout << "Building CSR matrix..." << std::endl;
    CSRMatrix web_graph = createCSRFromAdjacencyList(adj_list);
    
    // Normalize the matrix
    std::cout << "Normalizing transition matrix..." << std::endl;
    normalizeCSRMatrix(web_graph);
    
    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Run PageRank with Gauss-Seidel
    std::cout << "Running PageRank Gauss-Seidel algorithm..." << std::endl;
    std::vector<double> page_ranks = pageRankGaussSeidel(web_graph, alpha, threshold, 1000);
    
    // End timer
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    
    std::cout << "Total time: " << elapsed.count() << " seconds" << std::endl;
    
    // Print top 10 results if not too many nodes
    if (N <= 100) {
        std::cout << "PageRank results for all nodes:" << std::endl;
        for (int i = 0; i < web_graph.n; ++i) {
            std::cout << "Page " << i << ": " << page_ranks[i] << std::endl;
        }
    } else {
        // Create a vector of pairs (node index, rank value)
        std::vector<std::pair<int, double>> ranked_nodes;
        for (int i = 0; i < web_graph.n; ++i) {
            ranked_nodes.push_back({i, page_ranks[i]});
        }
        
        // Sort by rank value in descending order
        std::sort(ranked_nodes.begin(), ranked_nodes.end(), 
                 [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                     return a.second > b.second;
                 });
        
        // Print top 10 nodes
        std::cout << "Top 10 pages by PageRank:" << std::endl;
        for (int i = 0; i < std::min(10, web_graph.n); ++i) {
            std::cout << "Page " << ranked_nodes[i].first << ": " 
                      << ranked_nodes[i].second << std::endl;
        }
    }
    
    return 0;
}