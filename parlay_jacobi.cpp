#include "pagerank_utils.hpp"
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/primitives.h"
#include "parlaylib/include/parlay/sequence.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

double DAMPING = 0.85; // damping factor
double EPSILON = 1e-6; // convergence threshold
int MAX_ITERS = 100;   // maximum iterations

// pure gather-based PageRank (power iteration)
parlay::sequence<double> pagerank_power_iter(const Graph &G, double epsilon,
                                             double damping, int max_iters) {
  size_t n = G.n;
  double base_rank = (1.0 - damping) / n;

  auto old_ranks = parlay::sequence<double>(n, 1.0 / n);
  auto new_ranks = parlay::sequence<double>(n, 0.0);

  double error = 1.0;
  int iter = 0;
  while (error > epsilon && iter < max_iters) {
    // — gather from in-neighbors —
    parlay::parallel_for(0, n, [&](size_t i) {
      double sum = 0.0;
      for (size_t k = G.in_offsets[i]; k < G.in_offsets[i + 1]; ++k) {
        size_t j = G.in_edges[k];
        sum += old_ranks[j] / G.out_degrees[j];
      }
      new_ranks[i] = base_rank + damping * sum;
    });

    // — compute max-delta for convergence —
    auto deltas = parlay::tabulate<double>(
        n, [&](size_t i) { return std::fabs(new_ranks[i] - old_ranks[i]); });
    error = parlay::reduce(deltas, parlay::maxm<double>());

    // swap buffers
    std::swap(old_ranks, new_ranks);
    iter++;
  }

  std::cout << "Converged in " << iter << " iterations; final error = " << error
            << "\n";

  // normalize so sum = 1
  double sum = parlay::reduce(old_ranks, parlay::addm<double>());
  parlay::parallel_for(0, n, [&](size_t i) { old_ranks[i] /= sum; });

  return old_ranks;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr
        << "Usage: " << argv[0]
        << " <edge_list.txt> [num_vertices] [epsilon] [damping] [max_iters]\n";
    return 1;
  }
  std::string file = argv[1];
  size_t num_v = (argc > 2 ? std::stoull(argv[2]) : 0);
  double epsilon = (argc > 3 ? std::stod(argv[3]) : DEFAULT_EPSILON);
  double damping = (argc > 4 ? std::stod(argv[4]) : DEFAULT_DAMPING);
  int max_iters = (argc > 5 ? std::stoi(argv[5]) : DEFAULT_MAX_ITERS);

  std::cout << "Loading graph...\n";
  Graph G = Graph::load_and_build(file, num_v);
  std::cout << "Vertices: " << G.n << "   Edges: " << G.out_edges.size()
            << "\n";

  struct timeval t0, t1;
  gettimeofday(&t0, nullptr);
  auto ranks = pagerank_power_iter(G, epsilon, damping, max_iters);
  gettimeofday(&t1, nullptr);

  double secs = get_elapsed_time(t0, t1);
  std::cout << "Total PageRank time: " << secs << " seconds\n";

  ////////////////////////////////
  // Open the file in append mode
  std::ofstream output_file("pagerank_mc.txt", std::ios::app);
  if (!output_file.is_open()) {
    std::cerr << "Error opening output file for writing\n";
    return EXIT_FAILURE;
  }
  output_file << secs << "\n";
  output_file.close();
  ////////////////////////////////

  print_top_pages(ranks, G);
  return 0;
}
