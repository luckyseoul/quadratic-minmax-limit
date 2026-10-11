/*
 * Tight classes of a Paley conference matrix: frame identity and block-extension bounds.
 *
 * For a prime q = 1 (mod 4) let C be the symmetric Paley conference matrix of order
 * n = q+1 (C_00 = 0, C_0j = C_j0 = 1, C_{i+1,j+1} = chi(i-j) with chi the Legendre
 * symbol, zero diagonal).  Q_C(x) = sum_{i<j} C_ij x_i x_j and Phi(C) = max_x |Q_C(x)|.
 * A class is a pair {x,-x}; it is tight if |Q_C(x)| = Phi(C).  Let S be the set of
 * tight classes, N0 = |S|, and G = sum_{x in S} x x^T.
 *
 * The program prints, by exhaustive Gray-code enumeration of all 2^(n-1) classes:
 *   Phi, N0, whether G = N0 * I exactly, the extreme eigenvalues of G, and the resulting
 *   one-row bound sqrt(n * lambda_min(G) / N0) (equal to sqrt(n) for a frame),
 *   m1 = min over rows c in {+-1}^n of max_{x in S} |x.c|          (one-row cost on S),
 *   S1 = min over rows c of sum_{x in S} |x.c|, and S1/(N0 sqrt n)  (beta of the note).
 *
 * With a second argument "core" the first vertex is deleted before the enumeration, so the
 * order is n = q and the matrix is the Seidel matrix of the Paley graph (for q = 29 the
 * restriction of the order-30 matrix to 29 vertices).
 *
 * Build:  gcc -O3 -march=native -fopenmp -o one_vertex_tight_frame_paley \
 *             one_vertex_tight_frame_paley.c -lm        (OpenMP optional)
 * Run:    ./one_vertex_tight_frame_paley 13     (n=14, about 1 s)
 *         ./one_vertex_tight_frame_paley 17     (n=18, about 1 s)
 *         ./one_vertex_tight_frame_paley 29     (n=30, about 20 s on 8 cores)
 *         ./one_vertex_tight_frame_paley 29 core   (n=29, order-30 matrix minus one vertex)
 * Only the C standard library is used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAXN 32
#define MAXK 4096
static int n, q, C[MAXN][MAXN];
static int8_t S[MAXK][MAXN];
static int K = 0;

/* Eigenvalue range of the symmetric matrix G (cyclic Jacobi rotations). */
static void eig_range(long Gl[MAXN][MAXN], double *lo, double *hi) {
  double a[MAXN][MAXN];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) a[i][j] = (double)Gl[i][j];
  for (int sweep = 0; sweep < 100; sweep++) {
    double off = 0;
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
    if (off < 1e-18) break;
    for (int p = 0; p < n; p++)
      for (int r = p + 1; r < n; r++) {
        if (fabs(a[p][r]) < 1e-300) continue;
        double th = (a[r][r] - a[p][p]) / (2 * a[p][r]);
        double t = (th >= 0 ? 1.0 : -1.0) / (fabs(th) + sqrt(th * th + 1));
        double c = 1 / sqrt(t * t + 1), sn = t * c;
        for (int k = 0; k < n; k++) {
          double akp = a[k][p], akr = a[k][r];
          a[k][p] = c * akp - sn * akr;
          a[k][r] = sn * akp + c * akr;
        }
        for (int k = 0; k < n; k++) {
          double apk = a[p][k], ark = a[r][k];
          a[p][k] = c * apk - sn * ark;
          a[r][k] = sn * apk + c * ark;
        }
      }
  }
  *lo = 1e300;
  *hi = -1e300;
  for (int i = 0; i < n; i++) {
    if (a[i][i] < *lo) *lo = a[i][i];
    if (a[i][i] > *hi) *hi = a[i][i];
  }
}

static void build(void) {
  int qr[256];
  memset(qr, 0, sizeof qr);
  for (int i = 1; i < q; i++) qr[(i * i) % q] = 1;
  memset(C, 0, sizeof C);
  for (int j = 1; j < n; j++) C[0][j] = C[j][0] = 1;
  for (int i = 0; i < q; i++)
    for (int j = 0; j < q; j++) {
      int d = ((i - j) % q + q) % q;
      C[i + 1][j + 1] = (d == 0) ? 0 : (qr[d] ? 1 : -1);
    }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      int s = 0;
      for (int k = 0; k < n; k++) s += C[i][k] * C[k][j];
      if (s != (i == j ? q : 0)) { fprintf(stderr, "C^2 != qI\n"); exit(1); }
    }
}

/* Enumerate classes with x_0 = +1: the top PRE free bits are fixed by the prefix,
 * the low bits are traversed in Gray order.  mode 0: track max |Q|.  mode 1: collect
 * states with |Q| == target. */
static void enumerate(int mode, long target, long *maxq) {
  const int pre = (n - 1 > 8) ? 6 : 0;
  const int fb = n - 1 - pre;
  long gmax = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(max : gmax)
  for (int pfx = 0; pfx < (1 << pre); pfx++) {
    int x[MAXN], h[MAXN];
    long Q = 0;
    x[0] = 1;
    for (int i = 1; i <= fb; i++) x[i] = 1;
    for (int b = 0; b < pre; b++) x[fb + 1 + b] = ((pfx >> b) & 1) ? -1 : 1;
    for (int i = 0; i < n; i++) {
      h[i] = 0;
      for (int j = 0; j < n; j++) h[i] += C[i][j] * x[j];
    }
    for (int i = 0; i < n; i++) Q += x[i] * h[i];
    Q /= 2;
    uint64_t tot = 1ULL << fb;
    for (uint64_t g = 0; g < tot; g++) {
      if (g) {
        int i = 1 + __builtin_ctzll(g);
        Q -= 2L * x[i] * h[i];
        int xo = x[i];
        x[i] = -xo;
        for (int j = 0; j < n; j++) h[j] -= 2 * C[j][i] * xo;
      }
      long a = labs(Q);
      if (mode == 0) {
        if (a > gmax) gmax = a;
      } else if (a == target) {
#pragma omp critical
        {
          if (K >= MAXK) { fprintf(stderr, "too many tight classes\n"); exit(1); }
          for (int i = 0; i < n; i++) S[K][i] = (int8_t)x[i];
          K++;
        }
      }
    }
  }
  *maxq = gmax;
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) { fprintf(stderr, "usage: %s q [core]\n", argv[0]); return 1; }
  q = atoi(argv[1]);
  n = q + 1;
  if (n > MAXN || q % 4 != 1) { fprintf(stderr, "need prime q = 1 mod 4, q+1 <= %d\n", MAXN); return 1; }
  build();
  if (argc == 3 && strcmp(argv[2], "core") == 0) {
    for (int i = 0; i < q; i++)
      for (int j = 0; j < q; j++) C[i][j] = C[i + 1][j + 1];
    n = q;
  }
  long phi;
  enumerate(0, 0, &phi);
  long dummy;
  enumerate(1, phi, &dummy);
  printf("q=%d n=%d%s Phi=%ld N0=%d\n", q, n, n == q ? " (core)" : "", phi, K);
  static long G[MAXN][MAXN];
  for (int k = 0; k < K; k++)
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) G[i][j] += S[k][i] * S[k][j];
  int frame = 1;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (G[i][j] != (i == j ? K : 0)) frame = 0;
  printf("Gram of tight classes equals N0*I exactly: %s\n", frame ? "yes" : "no");
  double lo, hi;
  eig_range(G, &lo, &hi);
  printf("eigenvalues of Gram: min %.4f max %.4f (N0 = %d);  sqrt(n*min/N0) = %.4f\n", lo, hi, K, sqrt(n * lo / K));

  /* min over rows c (c_0 = +1) of max_x |x.c| and of sum_x |x.c| */
  static int16_t ST[MAXN][MAXK];
  for (int k = 0; k < K; k++)
    for (int i = 0; i < n; i++) ST[i][k] = S[k][i];
  const int pre = (n - 1 > 8) ? 6 : 0;
  const int fb = n - 1 - pre;
  long best1 = 1L << 60, bestinf = 1L << 60;
#pragma omp parallel for schedule(dynamic, 1) reduction(min : best1, bestinf)
  for (int pfx = 0; pfx < (1 << pre); pfx++) {
    int c[MAXN];
    int16_t ip[MAXK] __attribute__((aligned(64)));
    c[0] = 1;
    for (int i = 1; i <= fb; i++) c[i] = 1;
    for (int b = 0; b < pre; b++) c[fb + 1 + b] = ((pfx >> b) & 1) ? -1 : 1;
    for (int k = 0; k < K; k++) {
      int t = 0;
      for (int i = 0; i < n; i++) t += S[k][i] * c[i];
      ip[k] = (int16_t)t;
    }
    uint64_t tot = 1ULL << fb;
    for (uint64_t g = 0; g < tot; g++) {
      if (g) {
        int i = 1 + __builtin_ctzll(g);
        c[i] = -c[i];
        int16_t d = (int16_t)(2 * c[i]);
        for (int k = 0; k < K; k++) ip[k] += d * ST[i][k];
      }
      long s = 0, m = 0;
      for (int k = 0; k < K; k++) {
        long a = ip[k] < 0 ? -ip[k] : ip[k];
        s += a;
        if (a > m) m = a;
      }
      if (s < best1) best1 = s;
      if (m < bestinf) bestinf = m;
    }
  }
  printf("min_c max_x |x.c| = %ld   (one-row cost over the tight classes)\n", bestinf);
  printf("min_c sum_x |x.c| = %ld   avg = %.5f   avg/sqrt(n) = %.4f\n", best1, (double)best1 / K,
         (double)best1 / K / sqrt((double)n));
  return 0;
}
