#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <sys/time.h>
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/sequence.h"
#include "parlaylib/include/parlay/primitives.h"
#include "pagerank_utils.hpp"

// Delta-push PageRank implementation
parlay::sequence<double> pagerank_delta_push(const Graph &G,
                                            double epsilon,
                                            double damping,
                                            int max_iters) {
  size_t n = G.n;
  double base_rank = (1.0 - damping) / n;
  
  // Initialize ranks with equal values 1/n
  auto ranks = parlay::sequence<double>(n, 1.0/n);
  
  // Delta values for each vertex
  auto deltas = parlay::sequence<double>(n, 0.0);
  
  // Initially, all vertices need to push their entire rank
  parlay::parallel_for(0, n, [&](size_t i) {
    deltas[i] = ranks[i];
  });
  
  // Keep track of max error for convergence
  double max_delta = 1.0;
  int iter = 0;
  
  while (max_delta > epsilon && iter < max_iters) {
    // Next iteration's deltas
    auto next_deltas = parlay::sequence<double>(n, 0.0);
    
    // Process all vertices with non-negligible deltas
    parlay::parallel_for(0, n, [&](size_t i) {
      if (std::fabs(deltas[i]) > epsilon * 0.01) {  // Using a threshold for active vertices
        // Update the rank with the pending delta
        ranks[i] += deltas[i];
        
        // Calculate outgoing delta to distribute to neighbors
        double outgoing_delta = damping * deltas[i];
        
        // Skip if no outgoing edges
        if (G.out_degrees[i] > 0) {
          double delta_per_edge = outgoing_delta / G.out_degrees[i];
          
          // Push delta to out-neighbors
          for (size_t j = G.out_offsets[i]; j < G.out_offsets[i+1]; ++j) {
            size_t neighbor = G.out_edges[j];
            // Use atomic add to safely update next_deltas in parallel
            __atomic_add_fetch((long*)&next_deltas[neighbor], 
                              *(long*)&delta_per_edge, 
                              __ATOMIC_RELAXED);
          }
        }
        
        // Reset this vertex's delta
        deltas[i] = 0.0;
      }
    });
    
    // Calculate convergence metric
    auto delta_abs = parlay::tabulate<double>(n, [&](size_t i) {
      return std::fabs(next_deltas[i]);
    });
    max_delta = parlay::reduce(delta_abs, parlay::maxm<double>());
    
    // Swap deltas for next iteration
    std::swap(deltas, next_deltas);
    
    iter++;
    if (iter % 10 == 0) {
      std::cout << "Iteration " << iter << ", max_delta = " << max_delta << std::endl;
    }
  }
  
  std::cout << "Converged in " << iter << " iterations; final error = " << max_delta << "\n";
  
  // Apply base rank and normalize
  double sum = 0.0;
  parlay::parallel_for(0, n, [&](size_t i) {
    // Apply any remaining deltas
    ranks[i] += deltas[i];
    // Add base rank component
    ranks[i] = base_rank + damping * ranks[i];
    // Atomic add to sum for normalization
    __atomic_add_fetch((long*)&sum, *(long*)&ranks[i], __ATOMIC_RELAXED);
  });
  
  // Normalize so sum = 1
  parlay::parallel_for(0, n, [&](size_t i) {
    ranks[i] /= sum;
  });
  
  return ranks;
}

// Original power iteration implementation for comparison
parlay::sequence<double> pagerank_power_iter(const Graph &G,
                                             double epsilon,
                                             double damping,
                                             int max_iters) {
  size_t n = G.n;
  double base_rank = (1.0 - damping) / n;

  auto old_ranks = parlay::sequence<double>(n, 1.0/n);
  auto new_ranks = parlay::sequence<double>(n, 0.0);

  double error = 1.0;
  int iter = 0;
  while (error > epsilon && iter < max_iters) {
    // — gather from in-neighbors —
    parlay::parallel_for(0, n, [&](size_t i) {
      double sum = 0.0;
      for (size_t k = G.in_offsets[i]; k < G.in_offsets[i+1]; ++k) {
        size_t j = G.in_edges[k];
        sum += old_ranks[j] / G.out_degrees[j];
      }
      new_ranks[i] = base_rank + damping * sum;
    });

    // — compute max-delta for convergence —
    auto deltas = parlay::tabulate<double>(n, [&](size_t i){
      return std::fabs(new_ranks[i] - old_ranks[i]);
    });
    error = parlay::reduce(deltas, parlay::maxm<double>());

    // swap buffers
    std::swap(old_ranks, new_ranks);
    iter++;
  }

  std::cout << "Converged in " << iter << " iterations; final error = " << error << "\n";

  // normalize so sum = 1
  double sum = parlay::reduce(old_ranks, parlay::addm<double>());
  parlay::parallel_for(0, n, [&](size_t i){
    old_ranks[i] /= sum;
  });

  return old_ranks;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <edge_list.txt> [num_vertices] [epsilon] [damping] [max_iters] [algorithm]\n";
    std::cerr << "  algorithm: 0 = power iteration (default), 1 = delta push\n";
    return 1;
  }
  std::string file = argv[1];
  size_t num_v = (argc>2 ? std::stoull(argv[2]) : 0);
  double epsilon = (argc>3 ? std::stod(argv[3]) : DEFAULT_EPSILON);
  double damping = (argc>4 ? std::stod(argv[4]) : DEFAULT_DAMPING);
  int max_iters = (argc>5 ? std::stoi(argv[5]) : DEFAULT_MAX_ITERS);
  int algorithm = (argc>6 ? std::stoi(argv[6]) : 0);

  std::cout << "Loading graph...\n";
  Graph G = Graph::load_and_build(file, num_v);
  std::cout << "Vertices: " << G.n
            << "   Edges: " << G.out_edges.size() << "\n";

  struct timeval t0, t1;
  gettimeofday(&t0, nullptr);
  
  parlay::sequence<double> ranks;
  if (algorithm == 0) {
    std::cout << "Running PageRank with power iteration...\n";
    ranks = pagerank_power_iter(G, epsilon, damping, max_iters);
  } else {
    std::cout << "Running PageRank with delta-push approach...\n";
    ranks = pagerank_delta_push(G, epsilon, damping, max_iters);
  }
  
  gettimeofday(&t1, nullptr);
  double secs = get_elapsed_time(t0, t1);
  std::cout << "Total PageRank time: " << secs << " seconds\n";

  print_top_pages(ranks, G.n);
  return 0;
}