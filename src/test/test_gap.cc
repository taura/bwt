#include <math.h>
#include <stdio.h>
#include <string.h>
#include <map>

#include "typedefs.h"
#include "util.h"
#include "gap_array.h"

#define bwt_check(exp) bwt::check_(exp, #exp, __FILE__, __LINE__, __func__)

int test_gap_array(bwt::idx_t n, bwt::idx_t K, bwt::idx_t m, uint64_t seed) {
  /* generate random alphabets */
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };
  
  for (bwt::idx_t j = 0; j < m; j++) {
    bwt::idx_t a = (j == 0 ? 0 : nrand48(rg) % (n - 1));
    assert(a < n - 1);
    bwt::idx_t b = (j == 0 ? n : a + 2 + nrand48(rg) % (n - 1 - a));
    assert(b <= n);

    /* generate K distinct values */
    bwt::idx_t * distinct_values = new bwt::idx_t[K];
    for (bwt::idx_t i = 0; i < K; i++) {
      distinct_values[i] = a + nrand48(rg) % (b - a + 1);
    }

    bwt::gap_array gap;
    gap.init(a, b, 0);		/* granularity does not matter */
    std::map<bwt::idx_t,bwt::idx_t> counts;
    /* generate data and put it in the map */
    printf("%ld a=%ld, b=%ld, %ld elems\n", j, a, b, b - a + 1);
    bwt::idx_t * values = new bwt::idx_t [b - a + 1] - a;
    for (bwt::idx_t i = a; i <= b; i++) {
      bwt::idx_t k = nrand48(rg) % K;
      bwt::idx_t v = distinct_values[k];
      assert(a <= v);
      assert(     v <= b);
      values[i] = v;
      if (counts.count(v) == 0) 
	counts[v] = 0;
      counts[v]++;
    }
    /* insert all values to gap array concurrently */
#pragma omp parallel 
    {
#pragma omp barrier
#pragma omp for
      for (bwt::idx_t i = a; i <= b; i++) {
	gap.inc(values[i]);
      }
    }
    /* check */
    for (bwt::idx_t v = a; v <= b; v++) {
      if (counts.count(v) == 0) {
	if (!bwt_check(gap.get(v) == 0)) return 0;
      } else {
	if (!bwt_check(counts[v] == gap.get(v))) return 0;
      }
    }
    printf("%lu OK\n", j);
  }
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 100000);
  size_t K = (argc > 2 ? atol(argv[2]) : sqrt(n));
  size_t m = (argc > 3 ? atol(argv[3]) : 100);
  uint64_t seed = (argc > 3 ? atol(argv[4]) : 20);
  if (test_gap_array(n, K, m, seed)) {
    return 0;
  } else {
    return 1;
  }
}
