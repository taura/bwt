#include <stdio.h>
#include "wavelet_matrix.h"

int test_wavelet_matrix(bwt::idx_t n, bwt::idx_t m, 
			bwt::alpha_t alpha0, bwt::alpha_t alpha1,
			bwt::idx_t segment_sz) {
  unsigned short rg[3] = { 1, 2, 3 };
  /* generate random alphabets */
  bwt::alpha_t * s = new bwt::alpha_t[n];
  bwt::alpha_t * w = new bwt::alpha_t[n];
  for (bwt::idx_t i = 0; i < n; i++) {
    bwt::alpha_t a = alpha0 + (nrand48(rg) % (alpha1 - alpha0 + 1));
    s[i] = a;
  }
  bwt::wavelet_matrix wm;

  printf("building a wavelet matrix of %lu elements ...\n", n);
  bwt::tsc_t c0 = bwt::get_tsc();
  wm.init(s, n, w, alpha0, alpha1, segment_sz);
  bwt::tsc_t c1 = bwt::get_tsc();
  printf("took %llu cycles\n", c1 - c0);
  for (bwt::idx_t i = 0; i < m; i++) {
    bwt::alpha_t c = alpha0 + nrand48(rg) % (alpha1 - alpha0);
    bwt::idx_t p = (bwt::idx_t)nrand48(rg) % n;
    bwt::tsc_t c0 = bwt::get_tsc();
    bwt::idx_t x = wm.rank(c, p);
    bwt::tsc_t c1 = bwt::get_tsc();
    bwt::idx_t y = wm.rank_slow(c, p);
    if (x == y) {
      printf("[OK] %llu cycles to find rank(%d, %lu) = %lu\n", 
	     c1 - c0, c, p, x);
    } else {
      printf("[NG] %llu cycles to find rank(%d, %lu) = %lu != %lu\n", 
	     c1 - c0, c, p, x, y);
      return 0;
    }
  }
  delete s;
  printf("OK\n");
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 1000000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  size_t segment_sz = (argc > 3 ? atol(argv[3]) : 100);
  if (test_wavelet_matrix(n, m, 'a', 'z', segment_sz)) {
    return 0;
  } else {
    return 1;
  }
}
