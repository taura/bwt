#include <stdio.h>
#include <string.h>
#include <vector>

#include "util.h"
#include "bwt.h"

int test_bwt(bwt::idx_t n, bwt::idx_t m,
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
  
  /* bwt of T[a:b] into L */
  bwt::alpha_t * L = new bwt::alpha_t[n];
  /* inverse BWT */
  bwt::alpha_t * I = new bwt::alpha_t[n];
  /* workspace */
  bwt::alpha_t * W = new bwt::alpha_t[n];
  
  bwt::bwt_opt opt;
  opt.set_defaults();
  bwt::mallocator mem(opt);

  for (bwt::idx_t j = 0; j < m; j++) {
    /* always test the full string case */
    bwt::idx_t a = (j == 0 ? 0 : nrand48(rg) % n);
    bwt::idx_t b = (j == 0 ? n : a + nrand48(rg) % (n - a));
    if (a >= b) b = a + 1;
    assert(a < n);
    assert(a < b);
    assert(b <= n);
    bwt::random_init(L, n, rg);
    /* here we give ns == 0; this is wrong in other
       situations; here it's safe because we know 
       we do not use sa */
    /* use I for workspace */

    mem.reset(0, opt);
    bwt::bwt bwt = bwt_leaf(T, n, a, b, L, mem, opt);
    bwt::random_init(I, n, rg);
    bwt.init_extra(W, mem, opt);
    bwt.ibwt(I);
    /* check if we got the identical string back */
    if (bwt::check_equal(T, I, a, b)) {
      printf("%3lu OK\n", j);
    } else {
      printf("%3lu NG\n", j);
      return 0;
    }
    bwt.fini(mem, opt);
  }
  delete[] T;
  delete[] L;
  delete[] I;
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 10000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  bwt::alpha_t alpha0 = (argc > 3 ? (bwt::alpha_t)argv[3][0] : 'a');
  bwt::alpha_t alpha1 = (argc > 4 ? (bwt::alpha_t)argv[4][0] : 'z');
  uint64_t seed = (argc > 5 ? atol(argv[5]) : 20);
  if (test_bwt(n, m, alpha0, alpha1, seed)) {
    return 0;
  } else {
    return 1;
  }
}
