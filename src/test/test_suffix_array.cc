#include <stdio.h>
#include <string.h>
#include <vector>

#include "util.h"
#include "suffix_array.h"

#define bwt_check(exp) bwt::check_(exp, #exp, __FILE__, __LINE__, __func__)

int cmp_string(const bwt::alpha_t * a, const bwt::alpha_t * b, 
	       const bwt::alpha_t * e) {
  const bwt::alpha_t * p;
  const bwt::alpha_t * q;
  for (p = a, q = b; ; p++, q++) {
    bwt_check(p < e);
    bwt_check(q < e);
    if (*p < *q) return -1;
    if (*p > *q) return 1;
  }
  bwt_check(0);
}

/* generate a string of n alphabets, each between
   alpha0 (inclusive) and alpha1 (exclusive)
   and check the inequality of randomly chosen 
   m back-to-back pairs */
int test_sa(bwt::idx_t n, bwt::idx_t m, 
	    bwt::alpha_t alpha0, bwt::alpha_t alpha1) {
  unsigned short rg[3] = { 1, 2, 3 };
  bwt::alpha_t * T = new bwt::alpha_t[n];
  bwt::bwt_opt opt;
  opt.set_defaults(T, n);
  bwt::mallocator mem(n, opt);

  for (bwt::idx_t i = 0; i < n - 1; i++) {
    T[i] = alpha0 + (nrand48(rg) % (alpha1 - alpha0));
  }
  T[n - 1] = 0;
  bwt::idx_t * SA = new bwt::idx_t[n];
  if (!bwt::sa_range(T, n, 0, n, SA, mem)) return 0;
  for (bwt::idx_t i = 0; i < m; i++) {
    bwt::idx_t k = nrand48(rg) % (n - 1);
    if (!bwt_check(cmp_string(T + SA[k], T + SA[k + 1], T + n) < 0)) {
      return 0;
    }
  }
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 1000000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  if (test_sa(n, m, 'a', 'z' + 1)) {
    printf("OK\n");
    return 0;
  } else {
    printf("NG\n");
    return 1;
  }
}
