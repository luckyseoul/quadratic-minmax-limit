/* Signed covariance of the tight set of a Paley conference matrix of prime order q+1 (q = 1 mod 4 prime).
 *   build:  gcc -O2 -fopenmp -o paley_signed_tight_covariance paley_signed_tight_covariance.c
 *   run:    ./paley_signed_tight_covariance q          (q = 5 13 17 29 37; q=37 takes ~6 min on 8 cores)
 * Enumerates all classes {x,-x} (x_0 = +1) by Gray code, finds Phi = max |Q_C(x)| and accumulates over the tight
 * classes (|Q|=Phi): N0, the unsigned Gram G = sum x x^T and the signed matrix M = sum sigma(x) x x^T, sigma = sign Q(x).
 * Prints whether G = N0 I and whether M = (N0 Phi / N) C, N = n(n-1)/2, exactly (integer arithmetic). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#define MAXN 48
int n, q;
int8_t A[MAXN][MAXN];
static int legendre(int a, int p) { a %= p; if (a < 0) a += p; if (a == 0) return 0;
  int r = 1, e = (p - 1) / 2, b = a; while (e) { if (e & 1) r = r * b % p; b = b * b % p; e >>= 1; } return r == 1 ? 1 : -1; }
int main(int argc, char **argv) {
  q = atoi(argv[1]); n = q + 1;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
    if (i == j) A[i][j] = 0; else if (i == q || j == q) A[i][j] = 1; else A[i][j] = legendre(i - j, q); }
  int pre = (n >= 22) ? 6 : 0, fb = n - 1 - pre;
  long best = 0, N0 = 0; static long G[MAXN][MAXN], M[MAXN][MAXN];
  #pragma omp parallel
  {
    long lbest = 0, lN0 = 0; static __thread long lG[MAXN][MAXN], lM[MAXN][MAXN];
    memset(lG, 0, sizeof lG); memset(lM, 0, sizeof lM);
    #pragma omp for schedule(dynamic, 1)
    for (int pfx = 0; pfx < (1 << pre); pfx++) {
      int8_t x[MAXN]; int fld[MAXN];
      x[0] = 1; for (int i = 1; i <= fb; i++) x[i] = 1;
      for (int b = 0; b < pre; b++) x[fb + 1 + b] = ((pfx >> b) & 1) ? -1 : 1;
      for (int i = 0; i < n; i++) { int t = 0; for (int j = 0; j < n; j++) t += A[i][j] * x[j]; fld[i] = t; }
      long Q = 0; for (int i = 0; i < n; i++) Q += x[i] * fld[i]; Q /= 2;
      uint64_t total = 1ULL << fb;
      for (uint64_t g = 0; g < total; g++) {
        if (g) { int i = 1 + __builtin_ctzll(g); Q -= 2 * x[i] * fld[i]; x[i] = -x[i]; int d = 2 * x[i];
                 for (int j = 0; j < n; j++) fld[j] += d * A[j][i]; }
        long a = labs(Q);
        if (a >= lbest) {
          if (a > lbest) { lbest = a; lN0 = 0; memset(lG, 0, sizeof lG); memset(lM, 0, sizeof lM); }
          lN0++; int s = Q > 0 ? 1 : -1;
          for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { long v = (long)x[i] * x[j]; lG[i][j] += v; lM[i][j] += s * v; }
        }
      }
    }
    #pragma omp critical
    { if (lbest > best) { best = lbest; N0 = 0; memset(G, 0, sizeof G); memset(M, 0, sizeof M); }
      if (lbest == best) { N0 += lN0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { G[i][j] += lG[i][j]; M[i][j] += lM[i][j]; } } }
  }
  long N = (long)n * (n - 1) / 2; int okG = 1, okM = 1;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
    if (G[i][j] != (i == j ? N0 : 0)) okG = 0;
    /* M = (N0*Phi/N) C  <=>  M*N == N0*Phi*C */
    if (M[i][j] * N != N0 * best * (long)A[i][j]) okM = 0; }
  printf("q=%d n=%d Phi=%ld N0=%ld  Gram==N0*I: %s   signed M==(N0*Phi/N)*C: %s\n", q, n, best, N0, okG ? "yes" : "no", okM ? "yes" : "no");
  return 0;
}
