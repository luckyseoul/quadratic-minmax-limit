// Exhaustive one-vertex extension of the q=13 Paley conference signing.
// Finite diagnostic only: it evaluates min_b max_x (|Q_C(x)|+|b.x|).
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <omp.h>

int main() {
  constexpr int q = 13, n = q + 1, states = 1 << (n - 1);
  int C[n][n] = {};
  int chi[q] = {};
  for (int a = 1; a < q; ++a) chi[(a * a) % q] = 1;
  for (int a = 1; a < q; ++a) if (!chi[a]) chi[a] = -1;
  for (int j = 1; j < n; ++j) C[0][j] = C[j][0] = 1;
  for (int i = 0; i < q; ++i) for (int j = i + 1; j < q; ++j)
    C[i + 1][j + 1] = C[j + 1][i + 1] = chi[(j - i) % q];
  std::vector<int8_t> spin(states * n);
  std::vector<int16_t> energy(states);
  int phi = 0;
  for (int mask = 0; mask < states; ++mask) {
    auto x = &spin[mask * n]; x[0] = 1;
    for (int i = 1; i < n; ++i) x[i] = ((mask >> (i - 1)) & 1) ? 1 : -1;
    int v = 0;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) v += C[i][j] * x[i] * x[j];
    energy[mask] = v; phi = std::max(phi, std::abs(v));
  }
  std::atomic<int> best(1 << 30), best_mask(0);
  #pragma omp parallel for schedule(static)
  for (int mask = 0; mask < states; ++mask) {
    int b[n]; b[0] = 1;
    for (int i = 1; i < n; ++i) b[i] = ((mask >> (i - 1)) & 1) ? 1 : -1;
    int local = 0;
    for (int xmask = 0; xmask < states; ++xmask) {
      const auto x = &spin[xmask * n]; int dot = 0;
      for (int i = 0; i < n; ++i) dot += b[i] * x[i];
      local = std::max(local, std::abs((int)energy[xmask]) + std::abs(dot));
      if (local >= best.load(std::memory_order_relaxed)) break;
    }
    int old = best.load();
    while (local < old && !best.compare_exchange_weak(old, local)) {}
    if (local == best.load() && mask < best_mask.load()) best_mask.store(mask);
  }
  const int winner = best_mask.load();
  std::cout << "{\"n\":" << n << ",\"source_phi\":" << phi
            << ",\"extension_phi\":" << best.load() << ",\"b_mask\":" << winner
            << ",\"workers\":" << omp_get_max_threads() << "}\n";
}
