// Independent exhaustive checker for the q=13 Paley one-vertex scan.
// It deliberately uses a reduction, rather than the producer's atomic pruning.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <omp.h>

int main() {
  constexpr int q = 13, n = 14, states = 1 << 13;
  int chi[q] = {};
  for (int a = 1; a < q; ++a) chi[(a * a) % q] = 1;
  for (int a = 1; a < q; ++a) if (chi[a] == 0) chi[a] = -1;
  int C[n][n] = {};
  for (int j = 1; j < n; ++j) C[0][j] = C[j][0] = 1;
  for (int i = 0; i < q; ++i) for (int j = i + 1; j < q; ++j)
    C[i + 1][j + 1] = C[j + 1][i + 1] = chi[(j - i) % q];
  std::vector<int16_t> values(states);
  std::vector<int8_t> spins(states * n);
  int source = 0;
  for (int u = 0; u < states; ++u) {
    int8_t* x = &spins[u * n]; x[0] = 1;
    for (int i = 1; i < n; ++i) x[i] = ((u >> (i - 1)) & 1) ? 1 : -1;
    int e = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) e += C[i][j] * x[i] * x[j];
    values[u] = e; source = std::max(source, std::abs(e));
  }
  int global = 1 << 30;
  #pragma omp parallel for reduction(min:global) schedule(dynamic, 8)
  for (int u = 0; u < states; ++u) {
    int b[n]; b[0] = 1;
    for (int i = 1; i < n; ++i) b[i] = ((u >> (i - 1)) & 1) ? 1 : -1;
    int maximum = 0;
    for (int v = 0; v < states; ++v) {
      int dot = 0;
      for (int i = 0; i < n; ++i) dot += b[i] * spins[v*n + i];
      maximum = std::max(maximum, std::abs((int)values[v]) + std::abs(dot));
    }
    global = std::min(global, maximum);
  }
  std::cout << "{\"status\":\"PASS\",\"source_phi\":" << source
            << ",\"extension_phi\":" << global << ",\"workers\":"
            << omp_get_max_threads() << "}\n";
}
