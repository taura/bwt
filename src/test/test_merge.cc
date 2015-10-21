#include <stdio.h>
#include <string.h>
#include <vector>

#include "util.h"
#include "bwt.h"

#define bwt_check(exp) bwt::check_(exp, #exp, __FILE__, __LINE__, __func__)

int test_merge(bwt::idx_t n, bwt::idx_t m,
	       bwt::alpha_t alpha0, bwt::alpha_t alpha1,
	       uint64_t seed) {
  /* generate random alphabets */
  assert(alpha0 > 0);
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };
  /* since this is a merge example, we must
     have at least 2 elements (1 element from
     each bwt to merge) */
  bwt_check(n >= 2);
  bwt::alpha_t * T = new bwt::alpha_t[n];
  for (bwt::idx_t i = 0; i + 1 < n; i++) {
    T[i] = alpha0 + (nrand48(rg) % (alpha1 - alpha0 + 1));
  }
  T[n - 1] = 0;
  
  /* bwt of T[a:b] into L[a:b] */
  bwt::alpha_t * L = new bwt::alpha_t[n];
  /* bwt of T[a:c] into M[a:c] and T[c:b] into M[c:b] */
  bwt::alpha_t * M = new bwt::alpha_t[n];
  /* space to get the original string back from bwt */
  bwt::alpha_t * I = new bwt::alpha_t[n];
  
  for (bwt::idx_t j = 0; j < m; j++) {
    /* never forget to test the full string case */
    bwt::idx_t a = (j == 0 ? 0 : nrand48(rg) % (n - 1));
    assert(a < n - 1);
    bwt::idx_t b = (j == 0 ? n : a + 2 + nrand48(rg) % (n - 1 - a));
    assert(b <= n);
    assert(b - a >= 2);
    bwt::random_init(L, n, rg);
    bwt::idx_t c = (a + b) / 2;
    assert(a < c);
    assert(c < b);
    bwt::bwt_opt opt;
    opt.set_defaults();
    bwt::mallocator mem(opt);

    /* workspace to merge M[a:c] and M[c:b] into M[a:b] */
    bwt::alpha_t * W = mem.new_<bwt::alpha_t>(n, bwt::mem_reason_workspace_to_merge);
    bwt::bwt t0_ = bwt_leaf(T, n, a, c, M, mem, opt);
    bwt::bwt t1_ = bwt_leaf(T, n, c, b, M, mem, opt);
    bwt::bwt t01 = bwt_merge(t0_, t1_, W, mem, opt);
    bwt::stat.print();
    bwt::random_init(I, n, rg);
    t01.init_extra(W, mem, opt);
    t01.ibwt(I);
    t01.fini(mem, opt);
    //bwt::delete_(W, n, "workspace to merge");
    mem.delete_(W, n, bwt::mem_reason_workspace_to_merge);
    mem.print();

    if (bwt::check_equal(T, I, a, b)) {
      printf("%lu OK\n", j);
    } else {
      printf("%lu NG\n", j);
      return 0;
    }
  }
  delete[] T;
  delete[] L;
  delete[] M;
  delete[] I;
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 10000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  bwt::alpha_t alpha0 = (argc > 3 ? (bwt::alpha_t)argv[3][0] : 1);
  bwt::alpha_t alpha1 = (argc > 4 ? (bwt::alpha_t)argv[4][0] : 255);
  uint64_t seed = (argc > 5 ? atol(argv[5]) : 20);
  if (test_merge(n, m, alpha0, alpha1, seed)) {
    return 0;
  } else {
    return 1;
  }
}
