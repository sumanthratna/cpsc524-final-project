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


double DAMPING   = 0.85;    // damping factor
double EPSILON   = 1e-6;    // convergence threshold
int    MAX_ITERS = 100;     // maximum iterations

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


parlay::sequence<double> pagerank_gmres(const Graph &G,
                                        double epsilon,
                                        double damping,
                                        int max_iters) {
  size_t n = G.n;
  double base_rank = (1.0 - damping) / n;

  // Right-hand side vector b = (1 - c) * uniform
  auto b = parlay::sequence<double>(n, base_rank);
  auto x = parlay::sequence<double>(n, 1.0 / n);  // initial guess

  // Matrix-vector product A(x) = x - c * P^T x
  auto apply_A = [&](const parlay::sequence<double>& x_in) {
    parlay::sequence<double> out(n, 0.0);
    parlay::parallel_for(0, n, [&](size_t i) {
      double sum = 0.0;
      for (size_t k = G.in_offsets[i]; k < G.in_offsets[i+1]; ++k) {
        size_t j = G.in_edges[k];
        sum += x_in[j] / G.out_degrees[j];
      }
      out[i] = x_in[i] - damping * sum;
    });
    return out;
  };

  int restart = 20;
  int total_iters = 0;

  while (total_iters < max_iters) {
    auto r0 = b;
    auto Ax0 = apply_A(x);
    parlay::parallel_for(0, n, [&](size_t i) {
      r0[i] -= Ax0[i];
    });

    double beta = std::sqrt(parlay::reduce(parlay::tabulate(n, [&](size_t i) {
      return r0[i] * r0[i];
    }), parlay::addm<double>()));

    if (beta < epsilon) break;

    std::vector<parlay::sequence<double>> V;
    V.push_back(parlay::sequence<double>(n));
    parlay::parallel_for(0, n, [&](size_t i) {
      V[0][i] = r0[i] / beta;
    });

    std::vector<std::vector<double>> H(restart + 1, std::vector<double>(restart, 0.0));
    std::vector<double> cs(restart), sn(restart), e1(restart + 1, 0.0);
    e1[0] = beta;

    int actual_k = 0;
    for (int j = 0; j < restart && total_iters < max_iters; ++j, ++total_iters) {
      auto w = apply_A(V[j]);

      for (int i = 0; i <= j; ++i) {
        H[i][j] = parlay::reduce(parlay::tabulate(n, [&](size_t k) {
          return V[i][k] * w[k];
        }), parlay::addm<double>());
        parlay::parallel_for(0, n, [&](size_t k) {
          w[k] -= H[i][j] * V[i][k];
        });
      }

      H[j+1][j] = std::sqrt(parlay::reduce(parlay::tabulate(n, [&](size_t i) {
        return w[i] * w[i];
      }), parlay::addm<double>()));
      if (H[j+1][j] < 1e-12) break;

      auto v_next = parlay::sequence<double>(n);
      parlay::parallel_for(0, n, [&](size_t i) {
        v_next[i] = w[i] / H[j+1][j];
      });
      V.push_back(v_next);
      actual_k++;

      for (int i = 0; i < j; ++i) {
        double temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
        H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
        H[i][j] = temp;
      }

      double r = std::hypot(H[j][j], H[j+1][j]);
      cs[j] = H[j][j] / r;
      sn[j] = H[j+1][j] / r;
      H[j][j] = r;
      H[j+1][j] = 0.0;
      e1[j+1] = -sn[j] * e1[j];
      e1[j] = cs[j] * e1[j];

      if (std::abs(e1[j+1]) < epsilon) break;
    }

    std::vector<double> y(actual_k, 0.0);
    for (int i = actual_k - 1; i >= 0; --i) {
      y[i] = e1[i];
      for (int j = i + 1; j < actual_k; ++j)
        y[i] -= H[i][j] * y[j];
      y[i] /= H[i][i];
    }

    for (int j = 0; j < actual_k; ++j) {
      parlay::parallel_for(0, n, [&](size_t i) {
        x[i] += y[j] * V[j][i];
      });
    }
  }

  double sum = parlay::reduce(x, parlay::addm<double>());
  parlay::parallel_for(0, n, [&](size_t i) {
    x[i] /= sum;
  });

  std::cout << "GMRES finished in " << total_iters << " iterations\n";
  return x;
}


int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <edge_list.txt> [num_vertices] [epsilon] [damping] [max_iters]\n";
    return 1;
  }
  std::string file = argv[1];
  size_t num_v = (argc>2 ? std::stoull(argv[2]) : 0);
  if (argc>3) EPSILON   = std::stod(argv[3]);
  if (argc>4) DAMPING   = std::stod(argv[4]);
  if (argc>5) MAX_ITERS = std::stoi(argv[5]);

  std::cout << "Loading graph...\n";
  Graph G = Graph::load_and_build(file, num_v);
  std::cout << "Vertices: " << G.n
            << "   Edges: " << G.out_edges.size() << "\n";

  struct timeval t0, t1;
  gettimeofday(&t0, nullptr);
  auto ranks = pagerank_gmres(G, EPSILON, DAMPING, MAX_ITERS);
  gettimeofday(&t1, nullptr);

  double secs = (t1.tv_sec - t0.tv_sec)
              + (t1.tv_usec - t0.tv_usec)*1e-6;
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

  // print top 10
  std::vector<std::pair<double,size_t>> top;
  top.reserve(G.n);
  for (size_t i = 0; i < G.n; i++) 
    top.emplace_back(ranks[i], i);
  std::sort(top.begin(), top.end(),
            [&](auto &a, auto &b){ return a.first > b.first; });
  std::cout << "Top 10 pages by rank:\n";
  for (int i = 0; i < 10 && i < (int)top.size(); i++) {
    std::cout << (i+1) << ". Page " << top[i].second
              << " (rank=" << top[i].first << ")\n";
  }

  return 0;
}
