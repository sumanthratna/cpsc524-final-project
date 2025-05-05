// monte_carlo_pagerank.cpp
// Monte Carlo (random‐walk) PageRank in ParlayLib + OpenMP,
// adapted to your input/output style and timing conventions.
//
// Compile with:
//   g++ -std=c++17 -O3 -fopenmp \
//       -I parlaylib/include \
//       monte_carlo_pagerank.cpp -o pagerank_mc \
//       -pthread
//
// Run with:
//   ./pagerank_mc <edge_list.txt> [num_vertices]
//                  [num_walkers] [walk_length]
//                  [damping]     [num_threads]

#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

#include "pagerank_utils.hpp"
#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/primitives.h"
#include "parlaylib/include/parlay/sequence.h"

// Monte Carlo random‐walk PageRank
parlay::sequence<double> monte_carlo_pagerank(const Graph &G,
                                              size_t num_walkers,
                                              size_t walk_length,
                                              double damping, int num_threads) {
  size_t n = G.n;
  // per-thread local visit counts
  parlay::sequence<parlay::sequence<size_t>> local_counts(
      num_threads, parlay::sequence<size_t>(n, 0));

#pragma omp parallel num_threads(num_threads)
  {
    int tid = omp_get_thread_num();
    auto &counts = local_counts[tid];
    std::mt19937_64 rng((uint64_t)time(nullptr) ^ ((uint64_t)tid << 32));
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // round‐robin assign walkers
    for (size_t w = tid; w < num_walkers; w += num_threads) {
      size_t node = rng() % n; // random start
      for (size_t step = 0; step < walk_length; ++step) {
        if (dist(rng) < damping && G.out_degrees[node] > 0) {
          size_t off = G.out_offsets[node];
          size_t deg = G.out_offsets[node + 1] - off;
          node = G.out_edges[off + (rng() % deg)];
        } else {
          node = rng() % n; // teleport
        }
      }
      counts[node]++;
    }
  }

  // sum up into a global count vector
  parlay::sequence<size_t> counts(n);
  parlay::parallel_for(0, n, [&](size_t i) {
    size_t s = 0;
    for (int t = 0; t < num_threads; ++t)
      s += local_counts[t][i];
    counts[i] = s;
  });

  // normalize to get ranks
  parlay::sequence<double> ranks(n);
  double total = (double)num_walkers;
  parlay::parallel_for(0, n,
                       [&](size_t i) { ranks[i] = (double)counts[i] / total; });

  return ranks;
}

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0]
              << " <graph_filename> <N> <threshold> <damping> <num_threads>\n";
    return 1;
  }

  std::string file = argv[1];
  size_t num_v = std::stoull(argv[2]);
  double threshold = std::stod(argv[3]);
  double damping = std::stod(argv[4]);
  int num_threads = std::stoi(argv[5]);

  // Set reasonable defaults for Monte Carlo parameters
  size_t num_walkers = 1000000; // 1 million walkers
  size_t walk_length = 50;      // 50 steps per walk

  std::cout << "Loading graph...\n";
  Graph G = Graph::load_and_build(file, num_v);
  std::cout << "Vertices: " << G.n << "   Edges: " << G.out_edges.size()
            << "\n";

  struct timeval t0, t1;
  gettimeofday(&t0, nullptr);
  auto ranks =
      monte_carlo_pagerank(G, num_walkers, walk_length, damping, num_threads);
  gettimeofday(&t1, nullptr);

  double secs = get_elapsed_time(t0, t1);
  std::cout << "Total PageRank time: " << secs << " seconds\n";

  print_top_pages(ranks, G);
  return 0;
}
