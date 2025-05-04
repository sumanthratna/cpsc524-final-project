#pragma once

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

// Common constants
constexpr double DEFAULT_DAMPING = 0.85; // damping factor
constexpr double DEFAULT_EPSILON = 1e-6; // convergence threshold
constexpr int DEFAULT_MAX_ITERS = 100;   // maximum iterations

// Graph structure for PageRank
struct Graph {
  size_t n;
  parlay::sequence<size_t> out_degrees;
  parlay::sequence<size_t> out_offsets, out_edges;
  parlay::sequence<size_t> in_offsets, in_edges;

  static Graph load_and_build(const std::string &filename,
                              size_t num_vertices = 0) {
    std::ifstream file(filename);
    if (!file.is_open())
      throw std::runtime_error("Cannot open " + filename);

    // --- 1) read edge list into std::vector ---
    size_t src, dst, max_id = 0;
    std::vector<std::pair<size_t, size_t>> edges;
    std::string line;
    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream in(line);
      if (in >> src >> dst) {
        edges.emplace_back(src, dst);
        max_id = std::max(max_id, std::max(src, dst));
      }
    }
    file.close();

    size_t n = num_vertices > 0 ? num_vertices : (max_id + 1);
    Graph G;
    G.n = n;

    // --- 2) build outgoing CSR (deg, offsets, edges) ---
    std::vector<size_t> deg(n, 0);
    for (auto &e : edges)
      deg[e.first]++;

    std::vector<size_t> out_offsets(n + 1);
    out_offsets[0] = 0;
    for (size_t i = 0; i < n; i++)
      out_offsets[i + 1] = out_offsets[i] + deg[i];

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

    std::vector<size_t> in_offsets(n + 1);
    in_offsets[0] = 0;
    for (size_t i = 0; i < n; i++)
      in_offsets[i + 1] = in_offsets[i] + in_deg[i];

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
    G.out_offsets =
        parlay::sequence<size_t>(out_offsets.begin(), out_offsets.end());
    G.out_edges = parlay::sequence<size_t>(out_edges.begin(), out_edges.end());
    G.in_offsets =
        parlay::sequence<size_t>(in_offsets.begin(), in_offsets.end());
    G.in_edges = parlay::sequence<size_t>(in_edges.begin(), in_edges.end());

    return G;
  }
};

// Utility function to print top pages
inline void print_top_pages(const parlay::sequence<double> &ranks, size_t n,
                            size_t num_to_print = 10) {
  std::vector<std::pair<double, size_t>> top;
  top.reserve(n);
  for (size_t i = 0; i < n; i++)
    top.emplace_back(ranks[i], i);
  std::sort(top.begin(), top.end(),
            [&](auto &a, auto &b) { return a.first > b.first; });
  std::cout << "Top " << num_to_print << " pages by rank:\n";
  for (int i = 0; i < num_to_print && i < (int)top.size(); i++) {
    std::cout << (i + 1) << ". Page " << top[i].second
              << " (rank=" << top[i].first << ")\n";
  }
}

// Utility function to get elapsed time
inline double get_elapsed_time(const timeval &t0, const timeval &t1) {
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) * 1e-6;
}