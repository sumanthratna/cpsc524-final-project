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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

#include "parlaylib/include/parlay/parallel.h"
#include "parlaylib/include/parlay/sequence.h"
#include "parlaylib/include/parlay/primitives.h"



struct Graph {
    size_t n;  
    parlay::sequence<size_t> out_degrees;
    parlay::sequence<size_t> out_offsets, out_edges;
    parlay::sequence<size_t> in_offsets,  in_edges;
  
    static Graph load_and_build(const std::string &filename, size_t num_vertices = 0) {
      std::ifstream file(filename);
      if (!file.is_open()) throw std::runtime_error("Cannot open " + filename);
  
      // --- 1) read edge list into std::vector ---
      size_t src, dst, max_id = 0;
      std::vector<std::pair<size_t,size_t>> edges;
      std::string line;
      while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream in(line);
        if (in >> src >> dst) {
          edges.emplace_back(src,dst);
          max_id = std::max(max_id, std::max(src,dst));
        }
      }
      file.close();
  
      size_t n = num_vertices>0 ? num_vertices : (max_id+1);
      Graph G;  
      G.n = n;
  
      // --- 2) build outgoing CSR (deg, offsets, edges) ---
      std::vector<size_t> deg(n, 0);
      for (auto &e : edges) deg[e.first]++;
  
      std::vector<size_t> out_offsets(n+1);
      out_offsets[0] = 0;
      for (size_t i = 0; i < n; i++) 
        out_offsets[i+1] = out_offsets[i] + deg[i];
  
      std::vector<size_t> out_edges(edges.size());
      { // fill
        auto cursor = out_offsets;
        for (auto &e : edges) {
          size_t u = e.first, v = e.second;
          out_edges[cursor[u]++] = v;
        }
      }
  
      // --- 3) build reverse CSR (in-degrees, offsets, edges) ---
      std::vector<size_t> in_deg(n, 0);
      for (auto &e : edges) 
        in_deg[e.second]++;
  
      std::vector<size_t> in_offsets(n+1);
      in_offsets[0] = 0;
      for (size_t i = 0; i < n; i++) 
        in_offsets[i+1] = in_offsets[i] + in_deg[i];
  
      std::vector<size_t> in_edges(edges.size());
      { // fill
        auto cursor = in_offsets;
        for (auto &e : edges) {
          size_t u = e.first, v = e.second;
          in_edges[cursor[v]++] = u;
        }
      }
  
      // --- 4) move into Parlay sequences ---
      G.out_degrees = parlay::sequence<size_t>(deg.begin(), deg.end());
      G.out_offsets = parlay::sequence<size_t>(out_offsets.begin(), out_offsets.end());
      G.out_edges   = parlay::sequence<size_t>(out_edges.begin(),   out_edges.end());
      G.in_offsets  = parlay::sequence<size_t>(in_offsets.begin(),  in_offsets.end());
      G.in_edges    = parlay::sequence<size_t>(in_edges.begin(),    in_edges.end());
  
      return G;
    }
  };

// Monte Carlo random‐walk PageRank
parlay::sequence<double> monte_carlo_pagerank(const Graph &G,
                                              size_t num_walkers,
                                              size_t walk_length,
                                              double damping,
                                              int num_threads) {
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
      size_t node = rng() % n;           // random start
      for (size_t step = 0; step < walk_length; ++step) {
        if (dist(rng) < damping && G.out_degrees[node] > 0) {
          size_t off = G.out_offsets[node];
          size_t deg = G.out_offsets[node+1] - off;
          node = G.out_edges[off + (rng() % deg)];
        } else {
          node = rng() % n;               // teleport
        }
      }
      counts[node]++;
    }
  }

  // sum up into a global count vector
  parlay::sequence<size_t> counts(n);
  parlay::parallel_for(0, n, [&](size_t i) {
    size_t s = 0;
    for (int t = 0; t < num_threads; ++t) s += local_counts[t][i];
    counts[i] = s;
  });

  // normalize to get ranks
  parlay::sequence<double> ranks(n);
  double total = (double)num_walkers;
  parlay::parallel_for(0, n, [&](size_t i) {
    ranks[i] = (double)counts[i] / total;
  });

  return ranks;
}

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0]
              << " <graph_filename> <N> <threshold> <damping> <num_threads>\n";
    return 1;
  }

  std::string file       = argv[1];
  size_t      num_v      = std::stoull(argv[2]);
  double      threshold  = std::stod(argv[3]);
  double      damping    = std::stod(argv[4]);
  int         num_threads = std::stoi(argv[5]);

  // Set reasonable defaults for Monte Carlo parameters
  size_t num_walkers = 1000000;  // 1 million walkers
  size_t walk_length = 50;       // 50 steps per walk

  std::cout << "Loading graph...\n";
  Graph G = Graph::load_and_build(file, num_v);
  std::cout << "Vertices: " << G.n
            << "   Edges: " << G.out_edges.size() << "\n";

  struct timeval t0, t1;
  gettimeofday(&t0, nullptr);
  auto ranks = monte_carlo_pagerank(G,
                                    num_walkers,
                                    walk_length,
                                    damping,
                                    num_threads);
  gettimeofday(&t1, nullptr);

  double secs = (t1.tv_sec - t0.tv_sec)
              + (t1.tv_usec - t0.tv_usec)*1e-6;
  std::cout << "Total PageRank time: " << secs << " seconds\n";

  // print top-10 just like power-iter version
  std::vector<std::pair<double,size_t>> top;
  top.reserve(G.n);
  for (size_t i = 0; i < G.n; ++i)
    top.emplace_back(ranks[i], i);

  std::sort(top.begin(), top.end(),
            [&](auto &a, auto &b){ return a.first > b.first; });

  std::cout << "Top 10 pages by rank:\n";
  for (int i = 0; i < 10 && i < (int)top.size(); ++i) {
    std::cout << (i+1) << ". Page " << top[i].second
              << " (rank=" << top[i].first << ")\n";
  }

  return 0;
}
