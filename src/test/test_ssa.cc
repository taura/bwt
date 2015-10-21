#include <stdio.h>
#include <string.h>
#include <vector>

#include "util.h"
#include "bwt.h"

#define bwt_check(exp) bwt::check_(exp, #exp, __FILE__, __LINE__, __func__)

/* test sampled SA. see if sampled sa returns
   the same result with the real SA */
int test_ssa(bwt::idx_t n, bwt::idx_t m,
	     bwt::alpha_t alpha0, bwt::alpha_t alpha1,
	     uint64_t seed) {
  /* generate random alphabets */
  assert(alpha0 > 0);
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };
  bwt::alpha_t * T = new bwt::alpha_t[n];
  for (bwt::idx_t i = 0; i + 1 < n; i++) {
    T[i] = alpha0 + (nrand48(rg) % (alpha1 - alpha0 + 1));
  }
  T[n - 1] = 0;
  /* workspace */
  //bwt::alpha_t * W = new bwt::alpha_t[n];
  /* bwt of SA[a:b] into SA */
  bwt::idx_t * SA = new bwt::idx_t[n];
  /* bwt of T[a:b] into L[a:b] */
  bwt::alpha_t * L = new bwt::alpha_t[n];

  bwt::bwt_opt opt;
  opt.set_defaults(0, n);
  bwt::mallocator mem(n, opt);

  for (bwt::idx_t j = 0; j < m; j++) {
    /* always test the full string case */
    bwt::idx_t a = (j == 0 ? 0 : nrand48(rg) % n);
    bwt::idx_t b = (j == 0 ? n : a + nrand48(rg) % (n + 1 - a));
    if (b - a == 0) b = a + 1;
    assert(a < n);
    assert(a < b);
    assert(b <= n);
    bwt::random_init(L, n, rg);
    bwt::sa_range(T, n, a, b, SA, mem);
    bwt::bwt_opt opt;
    opt.set_defaults(T, n);

    bwt::bwt bwt = bwt_leaf(T, n, a, b, L, mem, opt);
    /* compare bs computes SA[r] for all r */
    for (bwt::idx_t r = a; r < b; r++) {
      if (!bwt_check(bwt.sa(r) == SA[r])) {
	printf("%lu NG\n", j);
	return 0;
      }
    }
    bwt.fini(mem, opt);
    printf("%lu OK\n", j);
  }
  delete[] T;
  delete[] SA;
  delete[] L;
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 10000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  bwt::alpha_t alpha0 = (argc > 3 ? (bwt::alpha_t)argv[3][0] : 'a');
  bwt::alpha_t alpha1 = (argc > 4 ? (bwt::alpha_t)argv[4][0] : 'z');
  uint64_t seed = (argc > 5 ? atol(argv[5]) : 20);
  if (test_ssa(n, m, alpha0, alpha1, seed)) {
    return 0;
  } else {
    return 1;
  }
}
