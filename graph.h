#pragma once
#include <string>
#include <vector>
#include <parlay/sequence.h>

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
