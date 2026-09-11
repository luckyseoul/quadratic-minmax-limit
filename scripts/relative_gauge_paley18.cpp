// Exact labeled relative-gauge occupancy for the balanced 9+9 Paley C18 split.
//
// This is a finite diagnostic for the good-fiber target, not an asymptotic
// theorem or a claim that either 9-vertex principal block is optimal.
// It builds shell transforms over G=P_9 x P_9 x {+/-1}, applies the exact
// Fourier product formula, and reports occupancies below the first odd energy
// strictly above 2^(3/2) m_9 = 24 sqrt(2), namely 35.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <omp.h>

namespace {
constexpr int P = 17;
constexpr int N = 18;
constexpr int HALF = 9;
constexpr int SPIN_BITS = HALF - 1;
constexpr int GROUP_BITS = 2 * HALF - 1;
constexpr int GROUP = 1 << GROUP_BITS;

using Vec = std::vector<std::int64_t>;

int legendre(int a) {
  a %= P;
  if (a < 0) a += P;
  if (a == 0) return 0;
  int value = 1;
  int base = a;
  int exponent = (P - 1) / 2;
  while (exponent) {
    if (exponent & 1) value = (value * base) % P;
    base = (base * base) % P;
    exponent >>= 1;
  }
  return value == 1 ? 1 : -1;
}

std::array<std::array<int, N>, N> conference() {
  std::array<std::array<int, N>, N> c{};
  for (int j = 1; j < N; ++j) c[0][j] = c[j][0] = 1;
  for (int i = 1; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) c[i][j] = c[j][i] = legendre((i - 1) - (j - 1));
  }
  return c;
}

template <std::size_t R, std::size_t C>
int quadratic(const std::array<std::array<int, C>, R>& a, int mask) {
  std::array<int, R> x{};
  x[0] = 1;
  for (std::size_t i = 1; i < R; ++i) x[i] = ((mask >> (i - 1)) & 1) ? -1 : 1;
  int total = 0;
  for (std::size_t i = 0; i < R; ++i) for (std::size_t j = i + 1; j < R; ++j) total += a[i][j] * x[i] * x[j];
  return total;
}

int rectangular(const std::array<std::array<int, HALF>, HALF>& b, int left, int right) {
  std::array<int, HALF> x{}, y{};
  x[0] = y[0] = 1;
  for (int i = 1; i < HALF; ++i) {
    x[i] = ((left >> (i - 1)) & 1) ? -1 : 1;
    y[i] = ((right >> (i - 1)) & 1) ? -1 : 1;
  }
  int total = 0;
  for (int i = 0; i < HALF; ++i) for (int j = 0; j < HALF; ++j) total += b[i][j] * x[i] * y[j];
  return total;
}

void fwht(Vec& a) {
  for (int width = 1; width < GROUP; width <<= 1) {
    for (int base = 0; base < GROUP; base += 2 * width) {
      for (int j = 0; j < width; ++j) {
        const auto u = a[base + j];
        const auto v = a[base + width + j];
        a[base + j] = u + v;
        a[base + width + j] = u - v;
      }
    }
  }
}

using Shells = std::map<int, Vec>;
Vec& bin(Shells& shells, int deficit) {
  auto [it, inserted] = shells.try_emplace(deficit);
  if (inserted) it->second.assign(GROUP, 0);
  return it->second;
}
void transform(Shells& shells) { for (auto& [_, v] : shells) fwht(v); }
std::int64_t mass(const Shells& shells) {
  std::int64_t answer = 0;
  for (const auto& [_, v] : shells) answer += std::accumulate(v.begin(), v.end(), std::int64_t{0});
  return answer;
}
}  // namespace

int main() {
  const auto c = conference();
  std::array<std::array<int, HALF>, HALF> left{}, right{}, cross{};
  for (int i = 0; i < HALF; ++i) for (int j = 0; j < HALF; ++j) {
    left[i][j] = c[i][j]; right[i][j] = c[i + HALF][j + HALF]; cross[i][j] = c[i][j + HALF];
  }
  constexpr int spins = 1 << SPIN_BITS;
  int ml = 0, mr = 0, mc = 0;
  std::vector<int> le(spins), re(spins);
  for (int x = 0; x < spins; ++x) {
    le[x] = quadratic(left, x); re[x] = quadratic(right, x);
    ml = std::max(ml, std::abs(le[x])); mr = std::max(mr, std::abs(re[x]));
  }
  std::vector<int> ce(spins * spins);
  for (int x = 0; x < spins; ++x) for (int y = 0; y < spins; ++y) {
    ce[x * spins + y] = rectangular(cross, x, y);
    mc = std::max(mc, std::abs(ce[x * spins + y]));
  }
  Shells ls, rs, cs;
  for (int x = 0; x < spins; ++x) for (int sign : {-1, 1}) {
    bin(ls, ml - sign * le[x])[x | ((sign < 0) << (GROUP_BITS - 1))]++;
    bin(rs, mr - sign * re[x])[(x << SPIN_BITS) | ((sign < 0) << (GROUP_BITS - 1))]++;
  }
  for (int x = 0; x < spins; ++x) for (int y = 0; y < spins; ++y)
    bin(cs, mc - std::abs(ce[x * spins + y]))[x | (y << SPIN_BITS)]++;
  const int ceiling = ml + mr + mc;
  const int target = 35;  // first odd integer strictly greater than 24 sqrt(2)
  const int cutoff = ceiling - target;
  const auto scalar_before = mass(ls) * mass(rs) * mass(cs);
  if (scalar_before != (std::int64_t{1} << (GROUP_BITS * 2))) throw std::runtime_error("state mass");
  transform(ls); transform(rs); transform(cs);
  Vec out(GROUP, 0);
  std::int64_t terms = 0;
  std::vector<std::array<const Vec*, 3>> triples;
  for (const auto& [dl, vl] : ls) for (const auto& [dr, vr] : rs) for (const auto& [dc, vc] : cs)
    if (dl + dr + dc <= cutoff) triples.push_back({&vl, &vr, &vc});
  terms = static_cast<std::int64_t>(triples.size());
  #pragma omp parallel for schedule(static)
  for (int index = 0; index < GROUP; ++index) {
    std::int64_t total = 0;
    for (const auto& t : triples) total += (*t[0])[index] * (*t[1])[index] * (*t[2])[index];
    out[index] = total;
  }
  const std::int64_t scalar = out[0];
  fwht(out);
  std::int64_t zeros = 0, minimum = INT64_MAX, maximum = 0;
  bool integral = true;
  for (auto value : out) {
    integral = integral && value % GROUP == 0;
    const auto occupancy = value / GROUP;
    zeros += occupancy == 0;
    if (occupancy) minimum = std::min(minimum, occupancy);
    maximum = std::max(maximum, occupancy);
  }
  std::cout << "{\"status\":\"PASS\",\"classification\":\"exact_fixed_order_relative_gauge_diagnostic_not_all_orders\""
            << ",\"order\":18,\"split\":[9,9],\"maxima\":[" << ml << ',' << mr << ',' << mc << ']'
            << ",\"target_energy\":35,\"independent_ceiling\":" << ceiling << ",\"deficit_cutoff\":" << cutoff
            << ",\"group_size\":" << GROUP << ",\"included_shell_triples\":" << terms
            << ",\"scalar_subthreshold_product_triples\":" << scalar << ",\"integral_inverse\":" << (integral ? "true" : "false")
            << ",\"empty_fibers\":" << zeros << ",\"minimum_positive_occupancy\":" << minimum
            << ",\"maximum_occupancy\":" << maximum << ",\"workers\":" << omp_get_max_threads() << "}" << std::endl;
  return integral ? 0 : 1;
}
